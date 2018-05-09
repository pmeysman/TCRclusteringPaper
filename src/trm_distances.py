#!/usr/local/bin/python3

import pandas as pd

naivedatadir = '../data/naive/'

naivedatafile = 'naive_VZVstudy.txt'

naivedata = pd.read_table(naivedatadir+naivedatafile)

trimercounts = dict()

total = 0

for seq in pd.unique(naivedata['CDR3']):
    total = total + 1
    for i in range(len(seq)-2):
        trimer = seq[i:i+3]
        trimercounts[trimer] = trimercounts.get(trimer, 0) + 1

print("Total: "+str(total))

trimercounts.update((x, 1-y/total) for x, y in trimercounts.items())

for trimer in trimercounts:
    if(trimercounts[trimer] < 0.95):
        print(trimer+' '+str(trimercounts[trimer]))
        
def compareTrimer(firstseq,secseq):
    count = 0
    norm = 0
    for i in range(len(firstseq)-2):
        trimer = firstseq[i:i+3]
        norm = norm + trimercounts.get(trimer, 1)
        if trimer in secseq:
            count = count + trimercounts.get(trimer, 1)
    return 1 - count/norm

