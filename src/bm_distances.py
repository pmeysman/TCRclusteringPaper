#!/usr/local/bin/python3

import pandas as pd

naivedatadir = '../data/naive/'

naivedatafile = 'naive_VZVstudy.txt'

naivedata = pd.read_table(naivedatadir+naivedatafile)

dimercounts = dict()

total = 0

for seq in pd.unique(naivedata['CDR3']):
    total = total + 1
    for i in range(len(seq)-1):
        dimer = seq[i:i+2]
        dimercounts[dimer] = dimercounts.get(dimer, 0) + 1

#print("Total: "+str(total))

dimercounts.update((x, 1-y/total) for x, y in dimercounts.items())

#for dimer in dimercounts:
#    if(dimercounts[dimer] < 0.95):
#        print(dimer+' '+str(dimercounts[dimer]))
        
def compareDimer(firstseq,secseq):
    count = 0
    norm = 0
    for i in range(len(firstseq)-1):
        dimer = firstseq[i:i+2]
        norm = norm + dimercounts.get(dimer, 1)
        if dimer in secseq:
            count = count + dimercounts.get(dimer, 1)
    return 1 - count/norm

