#!/usr/local/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import pandas as pd
import math

from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.metrics import confusion_matrix


class TCRclust:
    #Class to test a distance metric on a TCR sequence data set
    def __init__(self, tcrdata, size = 0):
        #Initialize with
        # tcrdata, a pandas dataframe containing the V, CDR3, J and epitope binding info (in that order!)
        # size, the number of TCR sequences to include, mostly so that the coe can be run quickly on fewer records
        
        #Size defaults to the entire dataframe
        if(size == 0):
            size = len(tcrdata.index)
            
        self.size = size
        self.tcrdata = tcrdata.iloc[range(size)]
        
        #Save the names of the columns so that the code remains understandable
        self.v = tcrdata.columns[0]
        self.cdr3 = tcrdata.columns[1]
        self.j = tcrdata.columns[2]
        self.epi = tcrdata.columns[3]
        
        
    
    def calc_dist(self,distfunct):
         # distfunct, a distance function that calculates a numerical distance between two tcrdata records
        
        def calc_row(recOne):
            if(recOne.name%10 == 0):
                print(recOne.name,flush=True)
            return self.tcrdata.apply(lambda recTwo: distfunct(recOne, recTwo, cdr3 = self.cdr3), axis=1)
        
        self.dist = self.tcrdata.apply(lambda recOne: calc_row(recOne), axis=1)

                
        
    def read_dist(self,filename):
        # instead of recalculating the distance matrix, just load in a precalculated one
        self.dist = pd.read_csv(filename,sep="\t",index_col=0)
        
        
    def cluster_aggcl(self, n_clusters = 10):
        self.clust = AgglomerativeClustering(n_clusters = n_clusters, affinity="precomputed", linkage="complete")
        try:
            self.clust.fit(self.dist)
        except:
            print("Invalid distance matrix")
    
            
    def cluster_dbscan(self, epsfract = 0.01, distance = None, min_samples = 2):
        if(distance == None):
            sort = np.unique(self.dist.as_matrix())
            distance = sort[int(len(sort)*epsfract)]
            if(distance == 0):
                distance = 10**-10 #DBSCAN does not like zero distances
        print("Using distance: "+str(distance))
        self.clust = DBSCAN(eps = distance, metric="precomputed", min_samples=min_samples)
        try:
            self.clust.fit(self.dist)
        except:
            print("Invalid distance matrix")
            
    def calc_confmat(self,rnd = 50):
        clstr = "Clusters"
        cnt = "count"
        res = pd.DataFrame({self.epi:self.tcrdata[self.epi],clstr:pd.Series(self.clust.labels_)})
        res[cnt] = 1
        res = res[res[clstr] != -1] #Remove samples not in a cluster= 
        res = res.groupby(clstr).filter(lambda x: len(x) > 1) #Remove clusters with a single sample (DBSCAN artefact)
        self.conf_mat = pd.pivot_table(res,values=cnt,index=[self.epi],columns=[clstr],aggfunc=np.sum,fill_value =0)
        self.entropy = sum(sc.stats.entropy(self.conf_mat.T)) #Entropy doesn't work so well in unbalanced cases
        try:
             self.chisquared = sc.stats.chi2_contingency(self.conf_mat)
        except:
             self.chisquared = [math.inf,0]
        self.acc = self.appr_acc(self.conf_mat) 
        self.prec = self.appr_prec(self.conf_mat) 
        
        if (rnd > 0): #Run random iterations
            self.rnd_conf_mat = dict()
            self.rnd_entropy_mat = np.zeros(rnd)
            self.rnd_chistat_mat = np.zeros(rnd)
            self.rnd_chipval_mat = np.zeros(rnd)
            self.rnd_acc_mat = np.zeros(rnd)
            self.rnd_prec_mat = np.zeros(rnd)
            for i in range(rnd):
                self.rndepi = np.random.permutation(self.tcrdata[self.epi])
                rndres = pd.DataFrame({self.epi:self.rndepi,clstr:pd.Series(self.clust.labels_)})
                rndres[cnt] = 1
                rndres = rndres[rndres[clstr] != -1]
                self.rnd_conf_mat[i] = pd.pivot_table(rndres,values=cnt,index=[self.epi],columns=[clstr],aggfunc=np.sum,fill_value =0)
                self.rnd_entropy_mat[i] = sum(sc.stats.entropy(self.rnd_conf_mat[i].T)) #Entropy doesn't work so well in unbalanced cases
                try:
                    chi_res = sc.stats.chi2_contingency(self.rnd_conf_mat[i])
                    self.rnd_chistat_mat[i] = chi_res[0]
                    self.rnd_chipval_mat[i] = chi_res[1]
                except:
                    self.rnd_chistat_mat[i] = math.inf
                    self.rnd_chipval_mat[i] = 1
                self.rnd_acc_mat[i] = self.appr_acc(self.rnd_conf_mat[i]) 
                self.rnd_prec_mat[i] = self.appr_prec(self.rnd_conf_mat[i]) 
            self.rnd_entropy = np.median(self.rnd_entropy_mat)
            self.rnd_chisquared = np.zeros(2)
            self.rnd_chisquared[0] = np.median(self.rnd_chistat_mat)
            self.rnd_chisquared[1] = np.median(self.rnd_chipval_mat)
            self.rnd_acc = np.median(self.rnd_acc_mat)
            self.rnd_prec = np.median(self.rnd_prec_mat)
        
        return self.conf_mat
    
    def appr_acc(self, conf_mat = None):
        #Method that pretends that we solved a supervised problem where each cluster corresponds to a single epitope
        #Returns the accuracy of the best solution
        
        #Run by default on the confusion matrix provided in the object
        if conf_mat is None:
            conf_mat = self.conf_mat
        
        #Define recursive function that finds the best fit for the diagonal
        def rec_max(mat):
            high  = mat.max().max()
            col = mat.max().idxmax()
            row = mat[col].idxmax()
            
            if(len(mat.index) > 1 and len(mat.columns) > 1):
                high = high + rec_max(mat.drop(row,0).drop(col,1))
            
            return high
        
        return rec_max(conf_mat)/self.size
        
    def appr_prec(self, conf_mat = None):
        #Method that estimates the precision of the solution
        #We assigned each cluster to the most common epitope
        #All other epitopes in the same cluster are considered false positives
        
        #Run by default on the confusion matrix provided in the object
        if conf_mat is None:
            conf_mat = self.conf_mat
        
        hits = np.sum(conf_mat.apply(np.max,axis=0))
        
        return hits/np.sum(conf_mat.values,axis=None)