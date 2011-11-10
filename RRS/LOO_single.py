#! /usr/local/bin/python
import numpy as np
import sys
from LOO import LOO_Hg

####################################################################################
# Main Function
####################################################################################

Masterfile = 'NatfishFinalAllobs_20110617_MNF.csv'
Connectionsfile = 'comparable_calcs_20110630.dat'
IDS_file = 'all_IDS.dat'
allIDS = np.genfromtxt(IDS_file,dtype=int)
ind_to_drop = int(sys.argv[1])
ID_to_drop = allIDS[ind_to_drop]
del allIDS

cHg,obsHg = LOO_Hg(ID_to_drop,Masterfile,Connectionsfile)

set1 = open('summaryRESULTS.dat','r').readlines()

ofp = open('data_{0:d}.dat'.format(ind_to_drop),'w')
for line in set1:
    ofp.write(line)
ofp.write('Dropped_ID--> {0}\n'.format(ID_to_drop))
ofp.write('calculated_Hg--> {0:f}\n'.format(cHg[0]))
ofp.write('Modeled_Hg--> {0:f}\n'.format(obsHg[0]))
ofp.close()
