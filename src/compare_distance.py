#!/usr/local/bin/python3

import pandas as pd
import numpy as np
import scipy as sp

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

recalculate = True


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




for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)
    
    if(recalculate):
        Clust.calc_dist(methods[methodName])
        print(Clust.dist.loc[1:10,1:10])
        Clust.dist.to_csv(path_or_buf=figdir+methodName+"_distmat.txt",sep="\t")
    else:
        Clust.read_dist(figdir+methodName+ "_distmat.txt")
    
    
    
    res_columns = ("Epi","meanEpi","meanOther","pval")
    results = pd.DataFrame(columns=res_columns)
    
    c = 0
    for epitope in Clust.tcrdata[Clust.epi].unique():
        print(epitope)
        matchLoc = Clust.tcrdata[Clust.epi] == epitope
        print(np.where(matchLoc)[0])
        if(len(np.where(matchLoc)[0])> 2):
            matchDist = Clust.dist.loc[np.where(matchLoc)[0],np.where(matchLoc)[0]].values
            matchDist = matchDist.flatten()
            matchDist = matchDist[np.nonzero(matchDist)]
            nomatchDist = Clust.dist.loc[np.where(matchLoc)[0],np.where(matchLoc == False)[0]].values
            nomatchDist = nomatchDist.flatten()
            nomatchDist = nomatchDist[np.nonzero(nomatchDist)]
            distWC = sp.stats.ranksums(nomatchDist,matchDist)
            print(distWC)
            results.loc[c] = [epitope,np.mean(matchDist),np.mean(nomatchDist),distWC[1]]
            c = c + 1
    
    
    totalWC = sp.stats.wilcoxon(results["meanEpi"],results["meanOther"])
    results.loc[c] = ["Total",np.mean(results["meanEpi"]),np.mean(results["meanOther"]),totalWC[1]]
    
    results.to_csv(path_or_buf=figdir+methodName+"_dist.txt",sep="\t")
