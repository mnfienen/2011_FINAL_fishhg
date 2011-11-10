'''
This code duplicates the sequence outlined in Prep.bat
'''


from sorting_gqsort import gqsort_Hgdat, gqsort_Hgdatext
from forward_NDMMF import calc_Hg
#debug
from LOO import drop_ID
#
import numpy as np
import os
import time

# ##
# This function calculates the Hg concentration for a left out observation
# Parameters are estimated with the observation (and any broken connections) left out
# then, the Hg concentration is calculated with those new, estimated parameters
# ##
Masterfile = 'NatfishFinalAllobs_20110617_MNF.csv'
Connectionsfile = 'comparable_calcs_20110630.dat'
ID_to_drop = 1

datfile = 'Hgdata.dat'
srtfile = 'Hgdata.srt'
extfile = 'Hgdata.datext'
extsrtfile  = 'Hgdata.datextsrt'
ndxfile = 'Hgdata.ndx'


# read in and parse the entire validation data set
MasterData = np.loadtxt(Masterfile,skiprows=1,delimiter=',')
headers = open(Masterfile,'r').readline().strip().split(',')
MKind = headers.index('ID')
ID = MasterData[:,MKind].astype(int)
SpC = MasterData[:,headers.index('SpC')].astype(int)
event = MasterData[:,headers.index('Event')].astype(int)
length = MasterData[:,headers.index('length')].astype(float)
Hg_obs = MasterData[:,headers.index('Hg')].astype(float)

# adjust the weights by a log transformation
Wt = MasterData[:,headers.index('Wt')].astype(float)
Wt[Wt<2]=1
Wt = np.log(Wt) + 1.0
MasterData[:,headers.index('Wt')]=Wt
# debug printing of weights
#np.savetxt('debug_weights',MasterData)

# first find the event, species, and length of the ID that is to be dropped - this is required
# later for the forward calculation
index_to_drop = np.nonzero(ID==ID_to_drop)[0]
csp = SpC[index_to_drop]
cev = event[index_to_drop]
clen = length[index_to_drop]
obsHg = Hg_obs[index_to_drop]
# find all IDs that need to be dropped
dropIDs = drop_ID(ID_to_drop,Connectionsfile)

MKlist = np.setdiff1d(ID,dropIDs)

# make a set out of the ID field
# SORT MasterData by ID
MasterData = MasterData[MasterData[:,MKind].argsort()]
MasterData[:,MKind] = MasterData[:,MKind].astype(int)



# excellent way to dereference the Master Data - requires that the MasterData matrix
# be sorted by MKind (done above outside the loop)
# for more details see: 
#http://stackoverflow.com/questions/5505380/most-efficient-way-to-pull-specified-rows-from-a-2-d-array
CurrMasterData = MasterData[np.searchsorted(MasterData[:,MKind],MKlist),:]
    

# now, write out the datfile in proper formats
# from Donato's code:
# SPC, Event, length, Result, DL, WT
ofp = open(datfile,'w')
for line in CurrMasterData:
    ofp.write('%3d %7d %13.8f %13.8f %2d %13.8f\n' 
              %(line[6],line[7],line[1],line[4],line[3],line[2]))
ofp.close()

# trim away any orphaned SpC or Events
allSPC = np.unique(CurrMasterData[:,6])
allEVENT = np.unique(CurrMasterData[:,7])
# now read in the parameter starting values files and trim out irrelevant parameters
spcdat = np.loadtxt('Hgspc.srt.master')
currspc = spcdat[np.searchsorted(spcdat[:,0],allSPC),:] # see above for ref. on this technique
ofp = open('Hgspc.srt','w')
for line in currspc:
    ofp.write('%10d %20f\n' %(line[0],line[1]))
ofp.close()

eventdat = np.loadtxt('Hgevent.srt.master')
currevent = eventdat[np.searchsorted(eventdat[:,0],allEVENT),:] # see above for ref. on this technique
ofp = open('Hgevent.srt','w')
for line in currevent:
    ofp.write('%10d %20f\n' %(line[0],line[1]))
ofp.close()


# memory cleanup
del CurrMasterData

# gqsort on Hgdata.dat, sorting by event, SPC, and DL
gqsort_Hgdat(datfile,srtfile)

# append a sequence number after sorting
# analagous to MLEprep01.c and write out to extfile

indat = np.loadtxt(srtfile)
ofp = open(extfile,'w')
i = 0
for line in indat:
    i += 1
    ofp.write('%3d %7d %13.8f %13.8f %2d  %13.8f %8d\n'
              %(line[0],
                line[1],
                line[2],
                line[3],
                line[4],
                line[5],
                i))
ofp.close()

# now sort by SPC
gqsort_Hgdatext(extfile,extsrtfile)

# Finally, strip out the index from extsrtfile and save in the .ndx file
indat = np.loadtxt(extsrtfile)
ofp = open(ndxfile,'w')
for line in indat:
    ofp.write('%8d\n' %(line[-1]))
ofp.close()



os.exit()
# call the external C-code Newton-Raphson parameter estimation code
os.system('./NRparest')    

# finally, read in the results and make the Hg prediction for the left-out value
# SpC parameters
SpCpars = np.loadtxt('BestSPs')
# Event parameters
Eventpars = np.loadtxt('BestEPs')

# pull the parameter values necessary for calculating Hg for the left-out value
spcind = np.nonzero(SpCpars[:,0]==csp)[0]
evind = np.nonzero(Eventpars[:,0]==cev)[0]

# calculate mercury for this index
cHg = calc_Hg(SpCpars[spcind,1],Eventpars[evind,1],clen)
cHg = (np.exp(cHg)-1)/1000.0
# return the left-out modeled Hg concentraion.
#return cHg, obsHg

    

  
    
    
    

