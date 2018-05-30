#!/usr/local/bin/python3

import pandas as pd
import numpy as np

method = "DBSCAN"

from TCRclust import *

from pt_tcr_distances import *
from vj_distances import vdist,jdist,levenshteinDistance
from trm_distances import compareTrimer, trimercounts
from bm_distances import compareDimer, dimercounts
from profile_distances import profile_distance_allprop

datadir = "../data/"

calcdir = "../results/Distances/"
recalculate = False

#root = "Dash"
root = "VDJ"

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
    return vdist(firstrec['Vgene'],secrec['Vgene'])+jdist(firstrec['Jgene'],secrec['Jgene'])
    
def distanceLvsh(firstrec, secrec, cdr3 = "CDR3"):
    return levenshteinDistance(firstrec[cdr3],secrec[cdr3])



file = root+"_known_dataset.txt"


#Reads in VDJ db file with headers Vgene, Jgene, CDR3 and Epitope
data = pd.read_table(datadir+file)


#Figure directory
figdir = "../results/"+method+"/"+root+"/singlethreshold/"

methodsVal = {"Length": 0.5,"Trimer": 0.34,"Dimer": 0.23,"GapAlign": 18,"Profile":1,"VJedit":0.408,"Levenshtein":1}
methods = {"Length": compareLength,"Trimer": distanceTrimer,"Dimer": distanceDimer,"GapAlign": distancePT,"Profile":distanceProfile,"VJedit":distanceVJ, "Levenshtein": distanceLvsh}
colordict = {"Length":"black","Trimer": 'tab:blue',"Dimer": 'tab:cyan',"GapAlign": 'tab:orange',"Profile":'tab:green',"VJedit": 'tab:brown',"Levenshtein":'tab:red'}


resultsMethods = {}




for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)
    
    if(recalculate):
        print("Please precalculate the distance matrix for " + methodName + " as a speed-up")
        Clust.calc_dist(methods[methodName])
        print(Clust.dist.loc[1:10,1:10])
    else:
        Clust.read_dist(calcdir + root + "/" + methodName + "_distmat.txt")
        
    res_columns = ("cl","tot","chi-stat","log-chi-pval","acc","ent","prec","rnd_chi-stat","rnd_log-chi-pval","rnd_acc","rnd_ent","rnd_prec")
    results = pd.DataFrame(columns=res_columns)


    if(method == "DBSCAN"):
        Clust.cluster_dbscan(distance = methodsVal[methodName])
    else:
        print("Error, this script only works for DBscan")
    Clust.calc_confmat()
    chi = Clust.chisquared
    acc = Clust.acc
    ent = Clust.entropy
    prec = Clust.prec
    print(Clust.conf_mat.values.sum())
    print(acc)
    print(prec)
    rnd_chi = Clust.rnd_chisquared
    rnd_acc = Clust.rnd_acc
    rnd_ent = Clust.rnd_entropy
    rnd_prec = Clust.rnd_prec
    results.loc[0] = [Clust.conf_mat.shape[1],Clust.conf_mat.values.sum(),chi[0],-np.log(chi[1]),acc,ent,prec,rnd_chi[0],-np.log(rnd_chi[1]),rnd_acc,rnd_ent,rnd_prec]
    
    bestChiEpsfract = results["log-chi-pval"].idxmax()

    print(Clust.conf_mat)
    Clust.conf_mat.to_csv(path_or_buf=figdir+methodName+"_confmat.txt",sep="\t")
    resultsMethods[methodName] = results
    
    clstr = "Clusters"
    res_cdr3 = pd.DataFrame({Clust.cdr3:Clust.tcrdata[Clust.cdr3],Clust.epi:Clust.tcrdata[Clust.epi],clstr:pd.Series(Clust.clust.labels_)})
    res_cdr3.sort_values(clstr,inplace=True)
    res_cdr3.to_csv(path_or_buf=figdir+methodName+"_clustoutput.txt",sep="\t")
    
    
for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["acc"],'o',label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["rnd_acc"],'.',label = "Random_"+methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Accuracy")
plt.savefig(figdir+"acc_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["chi-stat"],'o',label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["rnd_chi-stat"],'.',label = "Random_"+methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Chi test statistic")
plt.savefig(figdir+"chi-stat_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["log-chi-pval"],'o',label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["rnd_log-chi-pval"],'.',label = "Random_"+methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Log Chi-squared P-value")
plt.savefig(figdir+"log-chi-pval_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["ent"],'o',label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["rnd_ent"],'.',label = "Random_"+methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Entropy")
plt.savefig(figdir+"ent_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["prec"],'o',label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["rnd_prec"],'.',label = "Random_"+methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.savefig(figdir+"prec_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["tot"],resultsMethods[methodName]["cl"],label=methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Clusters")
plt.savefig(figdir+"cl_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()