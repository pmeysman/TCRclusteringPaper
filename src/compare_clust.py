#!/usr/local/bin/python3

import pandas as pd
import numpy as np

#Choose clustering method (only DBSCAN supported)
method = "DBSCAN"

from TCRclust import *

from pt_tcr_distances import *
from vj_distances import vdist,jdist,levenshteinDistance
from trm_distances import compareTrimer, trimercounts
from bm_distances import compareDimer, dimercounts
from profile_distances import profile_distance_allprop, profile_distance_allprop_realign

datadir = "../data/"

calcdir = "../results/Distances/"
recalculate = False

#Choose dataset to use from ../data
#root is also used as the directory in the results folder to group all results from the same dataset
root = "Dash"
#root = "VDJ"
file = root+"_known_dataset.txt"

if(method == "DBSCAN"):
    #Clustering thresholds as a fraction of the score distribution
    vec = (0.0001,0.0005,0.001,0.002,0.003,0.004,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.3,0.4)

#Reads in VDJ db file with headers Vgene, Jgene, CDR3 and Epitope
data = pd.read_table(datadir+file)

#Figure directory
figdir = "../results/"+method+"/"+root+"/"

#Define distance functions
#Should always take the same variables and return a numerical distance
def distanceTrimer (firstrec, secrec, cdr3 = "CDR3"):
    return compareTrimer(firstrec[cdr3],secrec[cdr3])

def distanceDimer (firstrec, secrec, cdr3 = "CDR3"):
    return compareDimer(firstrec[cdr3],secrec[cdr3])

def distancePT (firstrec, secrec, cdr3 = "CDR3"):
    return weighted_cdr3_distance( firstrec[cdr3], secrec[cdr3], default_distance_params)
    
def distanceProfile (firstrec, secrec, cdr3 = "CDR3"):
    return profile_distance_allprop(firstrec[cdr3],secrec[cdr3])
    
def distanceVJ (firstrec, secrec, cdr3 = "CDR3"):
    return vdist(firstrec['V gene'],secrec['V gene'])+jdist(firstrec['J gene'],secrec['J gene'])
    
def compareLength (firstrec, secrec, cdr3 = "CDR3"):
    return abs(len(firstrec[cdr3]) - len(secrec[cdr3]))
    
def distanceLvsh(firstrec, secrec, cdr3 = "CDR3"):
    return levenshteinDistance(firstrec[cdr3],secrec[cdr3])
    

#Color definitions
#If you add a new method, add a new color here!
colordict = {"Length":"black","Trimer": 'tab:blue',"Dimer": 'tab:cyan',"GapAlign": 'tab:orange',"Profile":'tab:green',"VJedit": 'tab:brown',"Levenshtein":'tab:red'}

#Method names and function assignments
#Adding any new methods will automatically be used in the figure
methods = {"Trimer": distanceTrimer,"Dimer": distanceDimer,"GapAlign": distancePT,"Profile":distanceProfile,"VJedit": distanceVJ, "Levenshtein": distanceLvsh}




resultsMethods = {}



#Try clustering with a distance methods
for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)
    
    if(recalculate):
        Clust.calc_dist(methods[methodName])
    else:
        Clust.read_dist(calcdir + root + "/" + methodName + "_distmat.txt")
        
    res_columns = ("cl","rec","chi-stat","log-chi-pval","acc","ent","prec","rnd_chi-stat","rnd_log-chi-pval","rnd_acc","rnd_ent","rnd_prec")
    results = pd.DataFrame(columns=res_columns)

    for i in vec:
        if(method == "DBSCAN"):
            Clust.cluster_dbscan(i)
        else:
            Clust.cluster_aggcl(i)
        Clust.calc_confmat()
        chi = Clust.chisquared
        acc = Clust.acc
        ent = Clust.entropy
        prec = Clust.prec
        print(prec)
        rnd_chi = Clust.rnd_chisquared
        rnd_acc = Clust.rnd_acc
        rnd_ent = Clust.rnd_entropy
        rnd_prec = Clust.rnd_prec
        results.loc[i] = [Clust.conf_mat.shape[1],Clust.conf_mat.values.sum()/Clust.size,chi[0],-np.log(chi[1]),acc,ent,prec,rnd_chi[0],-np.log(rnd_chi[1]),rnd_acc,rnd_ent,rnd_prec]
    
    bestChiEpsfract = results["log-chi-pval"].idxmax()
    if(method == "DBSCAN"):
        Clust.cluster_dbscan(bestChiEpsfract)
    else:
        Clust.cluster_aggcl(bestChiEpsfract)
    Clust.calc_confmat()
    print(Clust.conf_mat)
    resultsMethods[methodName] = results
    

#Make tons of plots

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["acc"],label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["rnd_acc"],label = "Random_"+methodName,linestyle='--', color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Accuracy")
plt.savefig(figdir+"acc_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["chi-stat"],label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["rnd_chi-stat"],label = "Random_"+methodName,linestyle='--', color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Chi test statistic")
plt.savefig(figdir+"chi-stat_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["log-chi-pval"],label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["rnd_log-chi-pval"],label = "Random_"+methodName,linestyle='--', color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Log Chi-squared P-value")
plt.savefig(figdir+"log-chi-pval_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["ent"],label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["rnd_ent"],label = "Random_"+methodName,linestyle='--', color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Entropy")
plt.savefig(figdir+"ent_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["prec"],label=methodName, color=colordict[methodName])
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["rnd_prec"],label = "Random_"+methodName,linestyle='--', color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.savefig(figdir+"prec_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()

for methodName in methods:
    plt.plot(resultsMethods[methodName]["rec"],resultsMethods[methodName]["cl"],label=methodName, color=colordict[methodName])

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Recall")
plt.ylabel("Clusters")
plt.savefig(figdir+"cl_comparison.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()