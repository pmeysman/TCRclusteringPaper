#!/usr/local/bin/python3

import pandas as pd
import numpy as np
import scipy as sp
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.decomposition import PCA

method = "DBSCAN"
from TCRclust import *

from pt_tcr_distances import *
from vj_distances import vdist,jdist, levenshteinDistance
from trm_distances import compareTrimer, trimercounts
from bm_distances import compareDimer, dimercounts
from profile_distances import profile_distance_allprop

datadir = "../data/"

#root = "Dash"
root = "VDJ"

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
    


methods = {"Length": compareLength,"Trimer": distanceTrimer,"Dimer": distanceDimer,"GapAlign": distancePT,"Profile":distanceProfile,"VJedit":distanceVJ, "Levenshtein": distanceLvsh}

resultsMethods = {}




for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)
    if(recalculate):
        Clust.calc_dist(methods[methodName])
    else:
        Clust.read_dist(figdir+methodName+ "_distmat.txt")
    
    distPCA = PCA(n_components=2).fit_transform(Clust.dist)
    
    colors = iter(cm.rainbow(np.linspace(0, 1, Clust.tcrdata[Clust.cdr3].str.len().max()-Clust.tcrdata[Clust.cdr3].str.len().min()+1)))
    for l in range(Clust.tcrdata[Clust.cdr3].str.len().max()-Clust.tcrdata[Clust.cdr3].str.len().min()+1):
        print(l + Clust.tcrdata[Clust.cdr3].str.len().min())
        matchLoc = Clust.tcrdata[Clust.cdr3].str.len() == l + Clust.tcrdata[Clust.cdr3].str.len().min()
        if(len(np.where(matchLoc)[0])> 0):
            plt.scatter(distPCA[np.where(matchLoc)[0], 0], distPCA[np.where(matchLoc)[0], 1], color=next(colors))
    
    plt.savefig(figdir+methodName+"_pca_length.pdf")
    plt.close()
    
    
    