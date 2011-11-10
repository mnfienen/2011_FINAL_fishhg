import numpy as np
import pickle
from Random_Sampling_Utils import sample_without_replacement

# set up a list of indices to throw out for repeated random sampling

# output filename for indices to dumpe out
outfile = 'all_out_inds.pkl'

#set number of realizations for repeated random sampling
num_realizations = 1000 
# set percent to leave out at each iteration - between 0 and 1
perc_out = 0.1
# infile from which to determine total number of indices available
infile = 'NatfishFinalAllobs_20110617.CSV'
indat = np.genfromtxt(infile,dtype=None,names=True,delimiter=',')

len_all_inds = len(indat['ID'])
all_inds = np.arange(len_all_inds)

cind_length = int(np.floor(len(all_inds)*perc_out))

all_out_indies = np.zeros((num_realizations,cind_length))

# for speedup, alias the function
_randsample = sample_without_replacement

for i in np.arange(num_realizations):
    print i
    all_out_indies[i,:] = np.array(list(_randsample(len_all_inds,cind_length)))

ofp = open(outfile,'wb')
pickle.dump(all_out_indies,ofp)
ofp.close()
