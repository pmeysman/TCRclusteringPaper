#!/usr/local/bin/python3

import numpy as np
import scipy as sc
import pandas as pd
import math

vjdatadir = '../data/imgt-reference/'

#use Levenshtein Distance (edit distance)
def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2 or c1 == '.' or c2 == '.':
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]
    

def readvjfile(filename):
    file = open(filename,'r')
    
    df = pd.DataFrame(columns=['acc','seq'])
    
    for cnt, line in enumerate(file):
        if(line[0] == '>'):
            header = line[1:-1]
            acc = header.split("|")[1]
            
            df = df.append({'acc':acc,'seq':''},ignore_index=True)

        else:
            df['seq'].iloc[-1] = df['seq'].iloc[-1] + line
        
    file.close()
    
    return df

def calcvjdist(df):
    
    def calc_row(recOne,df):
        if(recOne.name%10 == 0):
                print(recOne.name,flush=True)
        return df.apply(lambda recTwo: levenshteinDistance(recOne['seq'], recTwo['seq']), axis=1)
        
    return df.apply(lambda recOne: calc_row(recOne, df), axis=1)

vdatafile = vjdatadir+'tcrbv_jun2016.fq'
vdata = readvjfile(vdatafile)

vindex = {}
for index, row in vdata.iterrows():
    vindex[row['acc']] = index

vprecalc = pd.read_table('tcrbv_jun2016_dist.csv')

jdatafile = vjdatadir+'tcrbj_jun2016.fq'
jdata = readvjfile(jdatafile)

jindex = {}
for index, row in jdata.iterrows():
    jindex[row['acc']] = index

jprecalc = pd.read_table('tcrbj_jun2016_dist.csv')

def vdist(vFirst,vSecond):
    
    return vprecalc.values[vindex[vFirst],vindex[vSecond]]/vprecalc.values.max()

def jdist(jFirst,jSecond):
    try: 
        return jprecalc.values[jindex[jFirst],jindex[jSecond]]/jprecalc.values.max()
    except KeyError:
        jFirsts= jFirst.split(",")
        jSeconds= jSecond.split(",")
        
        minVal = math.inf;
        
        for jF in jFirsts:
            for jS in jSeconds:
                val = jprecalc.values[jindex[jF],jindex[jS]]
                if(val<minVal):
                    minVal = val
                
        return minVal/jprecalc.values.max()

if __name__ == '__main__':

    print(vdist("TRBV7-3*01","TRBV20-1*01"))

