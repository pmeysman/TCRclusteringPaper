#!/usr/local/bin/python3

import pandas as pd
import numpy as np
import scipy as sp
from sklearn import metrics
import matplotlib.pyplot as plt

method = "DBSCAN"
from TCRclust import *

from pt_tcr_distances import *
from vj_distances import vdist,jdist,levenshteinDistance
from trm_distances import compareTrimer, trimercounts
from bm_distances import compareDimer, dimercounts
from profile_distances import profile_distance_allprop

datadir = "../data/"

#root = "Dash"
root = "VDJ"

calcdir = "../results/Distances/"
recalculate = False


file = root+"_known_dataset.txt"


#Reads in VDJ db file with headers Vgene, Jgene, CDR3 and Epitope
data = pd.read_table(datadir+file)

#Figure directory
figdir = "../results/Distances/"+root+"/"

def distanceTrimer (firstrec, secrec, cdr3 = "CDR3"):
    return compareTrimer(firstrec[cdr3],secrec[cdr3])

def distanceDimer (firstrec, secrec, cdr3 = "CDR3"):
    return compareDimer(firstrec[cdr3],secrec[cdr3])

def distancePT (firstrec, secrec, cdr3 = "CDR3"):
    return weighted_cdr3_distance( firstrec[cdr3], secrec[cdr3], default_distance_params)
    
def distanceProfile (firstrec, secrec, cdr3 = "CDR3"):
    return profile_distance_allprop(firstrec[cdr3],secrec[cdr3])
    
def compareLength (firstrec, secrec, cdr3 = "CDR3"):
    return abs(len(firstrec[cdr3]) - len(secrec[cdr3]))
    
def distanceVJ (firstrec, secrec, cdr3 = "CDR3"):
    return vdist(firstrec['V gene'],secrec['V gene'])+jdist(firstrec['J gene'],secrec['J gene'])
    
def distanceLvsh(firstrec, secrec, cdr3 = "CDR3"):
    return levenshteinDistance(firstrec[cdr3],secrec[cdr3])
    


colordict = {"Length":"black","Trimer": 'tab:blue',"Dimer": 'tab:cyan',"GapAlign": 'tab:orange',"Profile":'tab:green',"VJedit": 'tab:brown',"Levenshtein":'tab:red'}

methods = {"Length": compareLength,"Trimer": distanceTrimer,"Dimer": distanceDimer,"GapAlign": distancePT,"Profile":distanceProfile,"VJedit":distanceVJ, "Levenshtein": distanceLvsh}

resultsMethods = {}

res_columns = ("Method","meanEpi","meanOther","MUstat","pval")
results = pd.DataFrame(columns=res_columns)
c = 0

for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)

    
    if(recalculate):
        Clust.calc_dist(methods[methodName])
    else:
        Clust.read_dist(calcdir + root + "/" + methodName + "_distmat.txt")
    
    
    
    posArr = []
    negArr = []
    
    
    for epitope in Clust.tcrdata[Clust.epi].unique():
        print(epitope)
        matchLoc = Clust.tcrdata[Clust.epi] == epitope
        print(np.where(matchLoc)[0])
        if(len(np.where(matchLoc)[0])> 0):
            matchDist = Clust.dist.loc[np.where(matchLoc)[0]][np.where(matchLoc)[0]].values
            
            for i in range(len(np.where(matchLoc)[0])-1):
                for j in range(i+1,len(np.where(matchLoc)[0])):
                    posArr.append(np.negative(matchDist[i,j])) #Signs flipped because greater distance = worse
                    
            nomatchDist = Clust.dist.loc[np.where(matchLoc)[0]][np.where(matchLoc == False)[0]].values
            negArr.extend(np.negative(nomatchDist.flatten()))
    
    print("Found "+str(len(posArr))+" positives")
    print("Found "+str(len(negArr))+" negatives")
    
    stat,pval=sp.stats.mannwhitneyu(posArr,negArr)
    
    results.loc[c] = [methodName,np.mean(posArr),np.mean(negArr),stat,pval]
    c = c + 1
    
    
    fullArr = posArr + negArr 
    fullLabel = ['P']*len(posArr) + ['N']*len(negArr)
    
    print("Total: "+str(len(fullLabel)))
    
    fpr, tpr, thresholds = metrics.roc_curve(fullLabel,fullArr,pos_label='P')
    resultsMethods[methodName] = {"fpr":fpr,"tpr":tpr,"thr":thresholds}
    
for methodName in methods:
    plt.plot(resultsMethods[methodName]["fpr"],resultsMethods[methodName]["tpr"],label=methodName, color=colordict[methodName])
plt.plot([0, 1], [0, 1], 'k--')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(figdir+"roc_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()


results.to_csv(path_or_buf=figdir+"mannU.txt",sep="\t")
