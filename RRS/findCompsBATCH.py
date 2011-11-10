from FishHgValidationUtilitiesBatch import read_and_parse_input, find_comps
import numpy as np
import time
import pickle
import zipfile

ubersttime = time.clock()
 
# initializations 
infile = 'NatfishFinalAllobs_20100526.CSV' # filename for the main input
inparfile = 'comps1.par'
# read parfile
compsrange = np.loadtxt(inparfile,dtype=int)
# set the range of all indices to knock out
indsrange = np.arange(compsrange[0],compsrange[1]+1)

# read and parse the input file
DL, MASTERKEY, SpC, Event = read_and_parse_input(infile)

# convert to arrays the lists which need to be referenced as arrays
SpC = np.array(SpC,dtype='int')
Event = np.array(Event,dtype='int')
MASTERKEY = np.array(MASTERKEY,dtype='int')
DL = np.array(DL, dtype='int')

total_length = len(MASTERKEY)
all_inds = np.arange(0,total_length)

#------------------------------------------------------------------------------------
# Now, run through the indices to knock out and save down the comparable MASTERKEYS
#------------------------------------------------------------------------------------
for cind in indsrange:
    currinds = np.nonzero(all_inds!=cind)[0]
    find_comps(SpC[currinds],Event[currinds],MASTERKEY[currinds],DL[currinds],cind)

endtime = time.clock()
print '\n' + '*'*25
print 'Final Complete elapsed time = {0} minutes'.format((endtime-ubersttime)/60) 