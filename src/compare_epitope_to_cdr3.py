#!/usr/local/bin/python3

import pandas as pd
import numpy as np
import scipy as sp

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

file = root+"_known_dataset.txt"

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



#Reads in VDJ db file with headers Vgene, Jgene, CDR3 and Epitope
data = pd.read_table(datadir+file)

#Figure directory
figdir = "../results/Distances/"+root+"/"

methods = {"Length": compareLength,"Trimer": distanceTrimer,"Dimer": distanceDimer,"GapAlign": distancePT,"Profile":distanceProfile, "VJedit":distanceVJ, "Levenshtein":distanceLvsh}


for methodName in methods:
    print("Using metric "+methodName)
    Clust = TCRclust(data)

    if(recalculate):
        print("Please precalculate the distance matrix for " + methodName)
        Clust.calc_dist(methods[methodName])
    else:
        Clust.read_dist(calcdir + root + "/" + methodName + "_distmat.txt")
    
    res_columns = ("Epitopes","epiDist","cdr3Dist","intraDist","sameMHC","sameClass")
    results = pd.DataFrame(columns=res_columns)
    c = 0
    
    for epitopePrime in Clust.tcrdata[Clust.epi].unique():
        
        epitopePrimeIndex = Clust.tcrdata[Clust.epi] == epitopePrime
        mhcPrime = Clust.tcrdata.loc[epitopePrimeIndex].iloc[0][4]
        mhcPrimeClass = 'I'
        if(len(mhcPrime) > 5):
            mhcPrimeClass = 'II'
        
        #Filter out the MHC class II presented peptides
        #This works because class II's have longer names
        #if((len(Clust.tcrdata.loc[epitopePrimeIndex].iloc[0][4])) > 5):
        #    continue
        
        #Filtering out smaller epitope datasets has no difference on the results
        if(len(Clust.tcrdata.loc[epitopePrimeIndex].index) > 0):
            print(epitopePrime+" ("+str(len(Clust.tcrdata.loc[epitopePrimeIndex].index))+")")
        
            for epitopeSec in Clust.tcrdata[Clust.epi].unique():
                if(not epitopePrime == epitopeSec):
                    epitopeSecIndex = Clust.tcrdata[Clust.epi] == epitopeSec
                    mhcSec = Clust.tcrdata.loc[epitopeSecIndex].iloc[0][4]
                    mhcSecClass = 'I'
                    if(len(mhcSec) > 5):
                        mhcSecClass = 'II'
                    
                    if(len(Clust.tcrdata.loc[epitopeSecIndex].index) > 0):
                        
                        #print('...'+epitopeSec+" ("+str(len(Clust.tcrdata.loc[epitopeSecIndex].index))+")")
                
                        #distance = weighted_cdr3_distance(epitopePrime, epitopeSec, default_distance_params)
                        distance = len(Clust.tcrdata.loc[epitopePrimeIndex].index)
                        if(len(Clust.tcrdata.loc[epitopeSecIndex].index) > distance):
                            distance = len(Clust.tcrdata.loc[epitopeSecIndex].index)
                        #print(distance)
                
                        subDist = Clust.dist.loc[np.where(epitopePrimeIndex)[0]][np.where(epitopeSecIndex)[0]].values
                
                
                        subMean = np.median(subDist)
                        #print(subMean)
                        
                        subPrimeDist = Clust.dist.loc[np.where(epitopePrimeIndex)[0]][np.where(epitopePrimeIndex)[0]].values
                        subPrimeMean = np.median(subPrimeDist)
                        
                        colorMean = subPrimeMean
                        
                        subSecDist = Clust.dist.loc[np.where(epitopeSecIndex)[0]][np.where(epitopeSecIndex)[0]].values
                        subSecMean = np.median(subSecDist)
                        
                        if(colorMean < subSecMean):
                            colorMean = subSecMean
                        
                        sameMHC = 0
                        if(mhcPrime is mhcSec):
                            sameMHC = 1
                        sameClass = 0
                        if(mhcPrimeClass is mhcSecClass):
                            sameClass = 1
                        
                        results.loc[c] = [epitopePrime+"-"+epitopeSec,distance,subMean,colorMean,sameMHC,sameClass]
                        c = c+1
                
    plt.scatter(results["epiDist"], results["cdr3Dist"], c=results["intraDist"])
    #plt.xlabel('Epitope distance')
    plt.xlabel('Number of distinct TCRs binding the epitope')
    plt.ylabel('CDR3 distance')
    
    pr, pval = sp.stats.spearmanr(results["epiDist"], results["cdr3Dist"])
    
    plt.title('Correlation: '+str(round(pr,2))+' (P-val = '+str('{:0.3e}'.format(pval))+')')
    
    #plt.savefig(figdir+"epitope_allclasses_median_density_"+methodName+".pdf")
    plt.savefig(figdir+"records_allclasses_median_density_"+methodName+".pdf")
    plt.close()
    
    mhcstats, mhcpval = sp.stats.mannwhitneyu(results[results['sameMHC'] == 1]['cdr3Dist'],results[results['sameMHC'] == 0]['cdr3Dist'])
    print('Same MHC MWU P-value : '+str('{:0.3e}'.format(mhcpval)))
    plt.figure()
    results.boxplot(column="cdr3Dist",by="sameMHC")
    plt.title('MHC distance for'+methodName+' (P-val = '+str('{:0.3e}'.format(mhcpval))+')')
    plt.savefig(figdir+"mhc_distance_"+methodName+".pdf")
    plt.close()
    
    classstats, classpval = sp.stats.mannwhitneyu(results[results['sameClass'] == 1]['cdr3Dist'],results[results['sameClass'] == 0]['cdr3Dist'])
    print('Same MHC Class MWU P-value : '+str('{:0.3e}'.format(classpval)))
    plt.figure()
    results.boxplot(column="cdr3Dist",by="sameClass")
    plt.title('MHC distance for'+methodName+' (P-val = '+str('{:0.3e}'.format(classpval))+')')
    plt.savefig(figdir+"mhcclass_distance_"+methodName+".pdf")
    plt.close()
    
    
    results.to_csv(path_or_buf=figdir+methodName+"_meandist.txt",sep="\t")