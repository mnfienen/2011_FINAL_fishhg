import numpy as np
from forward_NDMMF import Hg_all
from inv_mills_corrections import Inverse_Mills_Ratio_Correction

# get sigma
sigdat = np.loadtxt('summaryRESULTS.dat',skiprows=1)
sigmaIM = sigdat[0]



# SpC parameters
SpCpars = np.loadtxt('BestSPs')
# Event parameters
Eventpars = np.loadtxt('BestEPs')

Masterfile = 'NatfishFinalAllobs_20111116_MNF.csv'

allDat = np.genfromtxt(Masterfile,dtype=None,delimiter = ',', names = True)

ID = allDat['ID']
LEN = allDat['length']
Event = allDat['Event']
SpC = allDat['SpC']
DL = allDat['DL']
Hgobs = allDat['Hg']

Hgobs_corr = np.copy(Hgobs)
lnHg = Hg_all(SpCpars, Eventpars, LEN, Event, SpC, ID)
Hgmod = ((np.exp(lnHg)-1.0)/1000.0)


# apply the inverse Mills correction for the nondetects
DLinds = np.nonzero(DL)[0]
for i in DLinds:
    x,Eres = Inverse_Mills_Ratio_Correction(Hgobs[i],sigmaIM,Hgmod[i])
    Hgobs_corr[i] = x

ofp = open('AllHg_reference.dat','w')
ofp.write('%s,%s,%s,%s,%s\n' %('ID','DL','Hgmod','Hgobs','Hgobscorr'))
for i,CID in enumerate(ID):
    ofp.write('%d,%d,%f,%f,%f\n' %(CID,DL[i],Hgmod[i],Hgobs[i],Hgobs_corr[i]))
ofp.close()