import numpy as np
from forward_NDMMF import calc_Hg


# SpC parameters
SpCpars = np.loadtxt('BestSPs')
# Event parameters
Eventpars = np.loadtxt('BestEPs')

Masterfile = 'NatfishFinalAllobs_20110617_MNF.csv'

allDat = np.genfromtxt(Masterfile,dtype=None,delimiter = ',', names = True)

ID = allDat['ID']
LEN = allDat['length']
Event = allDat['Event']
SpC = allDat['SpC']

# now assemble the relevant indices
spcind_all  = []
for sp in CpC:
    spcind_all.append(np.nonzero(SpCpars[:,0]==sp)[0])
evind_all = []
for ev in Event:
    evind_all.append(np.nonzero(Eventpars[:,0]==ev)[0])
evind_all = np.squeeze(np.array(evind_all))
spcind_all = np.squeeze(np.array(spcind_all))

lnHg = []
for ii,csp in enumerate(SpC):
    lnHg.append(calc_Hg(SpCpars[]))