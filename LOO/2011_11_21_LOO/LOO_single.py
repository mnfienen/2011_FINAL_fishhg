#!/usr/bin/python
import numpy as np
import sys
from LOO import LOO_Hg

####################################################################################
# Main Function
####################################################################################

Masterfile = 'NatfishFinalAllobs_20111116_MNF.csv'
Connectionsfile = 'COMPARABILITY_2011116.dat'
IDS_file = 'all_IDS.dat'
allIDS = np.genfromtxt(IDS_file,dtype=int)
ind_to_drop = int(sys.argv[1])
ID_to_drop = allIDS[ind_to_drop]
del allIDS

cHg,obsHg = LOO_Hg(ID_to_drop,Masterfile,Connectionsfile)

set1 = open('summaryRESULTS.dat','r').readlines()

ofp = open('data_{0:d}.dat'.format(ind_to_drop),'w')
ofp.write('Dropped_ID--> {0}\n'.format(ID_to_drop))
for line in set1:
    ofp.write(line)
ofp.write('%20s%20s%20s\n' %('dropped_ID','modHg','measuredHg'))
ofp.write('%20d%20.8f%20.8f\n' %(ID_to_drop,cHg[0],obsHg[0]))
ofp.close()
