import comps_utils as cu
import numpy as np
import time
import pickle
import zipfile

ubersttime = time.time()
 
# initializations 
infile = 'NatfishFinalAllobs_20110513.CSV' # filename for the main input
inparfile = 'comps1.par'
# read parfile
compsrange = np.loadtxt(inparfile,dtype=int)
# set the range of all indices to knock out
indsrange = [92861]

# read and parse the input file
DL, MASTERKEY, SpC, Event = cu.read_and_parse_input(infile)

# convert to arrays the lists which need to be referenced as arrays
SpC = np.array(SpC,dtype='int')
Event = np.array(Event,dtype='int')
MASTERKEY = np.array(MASTERKEY,dtype='int')
DL = np.array(DL, dtype='int')

total_length = len(MASTERKEY)
all_inds = np.arange(0,total_length)
currresults = []
#------------------------------------------------------------------------------------
# Now, run through the indices to knock out and save down the comparable MASTERKEYS
#------------------------------------------------------------------------------------
for cind in indsrange:
    currinds = np.nonzero(all_inds!=cind)[0]
    currresults.append(
        cu.find_comps(SpC[currinds],Event[currinds],MASTERKEY[currinds],DL[currinds],cind,MASTERKEY[cind])
        )
endtime = time.time()
print '\n' + '*'*25
print 'Final Complete elapsed time = {0} minutes'.format((endtime-ubersttime)/60) 2