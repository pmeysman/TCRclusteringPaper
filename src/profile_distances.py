#!/usr/bin/env python3

import numpy as np

from scipy.stats import entropy
from scipy.stats import zscore


    
#Z-normalization of dictionaries
def znorm_scipy(d):
    keys, vals = zip(*d.items())
    return dict(zip(keys, zscore(vals, ddof=1)))
    
#Weighted euclidean
def weightedL2(a,b,w):
    q = a-b
    return np.sqrt((w*q*q).sum())


# chemical properties
basicity = {'A': 206.4, 'B': 210.7, 'C': 206.2, 'D': 208.6, 'E': 215.6, 'F': 212.1, 'G': 202.7,
            'H': 223.7, 'I': 210.8, 'K': 221.8, 'L': 209.6, 'M': 213.3, 'N': 212.8, 'P': 214.4,
            'Q': 214.2, 'R': 237.0, 'S': 207.6, 'T': 211.7, 'V': 208.7, 'W': 216.1, 'X': 210.2,
            'Y': 213.1, 'Z': 214.9}

hydrophobicity = {'A': 0.16, 'B': -3.14, 'C': 2.50, 'D': -2.49, 'E': -1.50, 'F': 5.00, 'G': -3.31,
                      'H': -4.63, 'I': 4.41, 'K': -5.00, 'L': 4.76, 'M': 3.23, 'N': -3.79, 'P': -4.92,
                      'Q': -2.76, 'R': -2.77, 'S': -2.85, 'T': -1.08, 'V': 3.02, 'W': 4.88, 'X': 4.59,
                      'Y': 2.00, 'Z': -2.13}

helicity = {'A': 1.24, 'B': 0.92, 'C': 0.79, 'D': 0.89, 'E': 0.85, 'F': 1.26, 'G': 1.15, 'H': 0.97,
                'I': 1.29, 'K': 0.88, 'L': 1.28, 'M': 1.22, 'N': 0.94, 'P': 0.57, 'Q': 0.96, 'R': 0.95,
                'S': 1.00, 'T': 1.09, 'V': 1.27, 'W': 1.07, 'X': 1.29, 'Y': 1.11, 'Z': 0.91}

elektrochem = {'basicity': znorm_scipy(basicity), 'hydrophobicity': znorm_scipy(hydrophobicity), 'helicity': znorm_scipy(helicity)}

def make_profile(sequence, prop = 'basicity'):
    return [elektrochem[prop][x] for x in sequence]

def distribution_similarity_realign(distriA, distriB):
    
    """Calculate the similarity between two distributions.
    Similarity is based on the distance measure used in the calc_distance function.
    If the distributions are of unequal length, slide the shorter one over the longer one
    and return the best similarity match.
    Return the best similarity score and the start position of the shorter distribution
    on the longer one that resulted in the best similarity score.
    """
    
    # if distributions are of equal length, return similarity
    if len(distriA)==len(distriB):
        return calc_distance(distriA, distriB), 0
    
    # if distributions are of unequal length
    else:
        # determine which is the longer distribution
        if len(distriA) > len(distriB):
            long_distri = distriA
            short_distri= distriB
        else:
            long_distri = distriB
            short_distri= distriA
        
        # slide shorter distribution over longer one, calculating the similarity between them
        # initialize first similarity value
        #print('Calculating initial similarity for position 0:\n{}\n{}.'.format(short_distri, long_distri[0:len(short_distri)]))
        similarity = calc_distance(short_distri, long_distri[0:len(short_distri)])
        pos = 0
        for i in range(1, len(long_distri)-len(short_distri)+1):
            #print('Calculating similarity for position {}:\n{}\n{}.'.format(i,short_distri, long_distri[i:i+len(short_distri)]))
            new_similarity = calc_distance(short_distri, long_distri[i:i+len(short_distri)])
            if new_similarity < similarity:
                similarity = new_similarity
                pos = i
        return similarity, pos
        
def distribution_similarity(distriA, distriB):
    
    # if distributions are of equal length, return similarity
    if len(distriA)==len(distriB):
        return calc_distance(distriA, distriB), 0
    
    else:
        # determine which is the longer distribution
        if len(distriA) > len(distriB):
            long_distri = distriA
            short_distri= distriB
        else:
            long_distri = distriB
            short_distri= distriA
        
        diff = len(long_distri) - len(short_distri)
        
        if diff % 2 == 0:
            shift = int(diff/2)
            return calc_distance(short_distri, long_distri[shift:len(long_distri)-shift]), diff/2
        
        else:
            
            #print("Odd distance : ",diff)
            shift = int(diff/2)
            #print("Shift : ",shift)
            
            #print("Short length : ", len(short_distri))
            #print("Old length : ", len(long_distri))
            #print("New length : ", len(long_distri[shift+1:len(long_distri)-shift]))
            
            first = calc_distance(short_distri, long_distri[shift+1:len(long_distri)-shift])
            second = calc_distance(short_distri, long_distri[shift:len(long_distri)-shift-1])
            
            if(first > second):
                return first,int(diff/2)+1
            else:
                return second,int(diff/2)
        
    
def calc_distance(distriA, distriB):
    
    """ Take two distributions of equal length and return their similarity. 
    Similarity is calculated as the Euclidean distance
    between the points in the distribution. """
    
    w = distriA.copy() #This is quicker than making an empty matrix
    middle = len(distriA)/2
    
    for i in range(len(distriA)):
        w[i] = (middle - abs(middle-i))+1
        
    tot = w.sum()
        
    wn = [weight/tot for weight in w]
    
    return weightedL2(distriA, distriB,wn)

        
def profile_distance(seqA,seqB,prop='hydrophobicity',st = 0,l = 0):
    distance, position = distribution_similarity(np.asarray(make_profile(seqA[start:len(seqA)-l], prop)),
                                                 np.asarray(make_profile(seqB[start:len(seqB)-l], prop)))
    return distance
    
def profile_distance_allprop(seqA,seqB,st = 0,l = 0):
    distanceALL = 0
    for prop in elektrochem:
        distance, position = distribution_similarity(np.asarray(make_profile(seqA[st:len(seqA)-l], prop)),
                                                     np.asarray(make_profile(seqB[st:len(seqB)-l], prop)))
        distanceALL = distanceALL + distance
    return distanceALL
    
def profile_distance_allprop_realign(seqA,seqB,st = 0,l = 0):
    distanceALL = 0
    for prop in elektrochem:
        distance, position = distribution_similarity_realign(np.asarray(make_profile(seqA[st:len(seqA)-l], prop)),
                                                     np.asarray(make_profile(seqB[st:len(seqB)-l], prop)))
        distanceALL = distanceALL + distance
    return distanceALL
    
if __name__ == '__main__':
    
    # Calculate distance between two CDR3 sequences as follows:
    CDR3_A = 'CASSLWTGSHEQYF'
    CDR3_B = 'CASSLWTGSHEQYF'
    #CDR3_B = 'CSARDRTGNGYTF'
    distance = profile_distance_allprop(CDR3_A,CDR3_B)
    print(distance)
