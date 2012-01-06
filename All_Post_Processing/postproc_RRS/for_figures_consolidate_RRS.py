import numpy as np
import sys



# set number of realizations
num_rel = 1000


# parse all the results files

# first read in the master ID list
ID_master = np.loadtxt('../../RRS/RRS_RESULTS_12_15_2011/all_IDs.dat',skiprows = 1, dtype=int)


# initialize a dictionary to hold the results
allHg = dict([(x,list()) for x in ID_master])

for crel in xrange(num_rel):
    print 'Running realization %d of %d:' %(crel+1,num_rel)
    indat = np.genfromtxt('../../RRS/RRS_RESULTS_12_15_2011/Hg_%d.dat' %(crel),dtype=None,names=True,skiprows=1)
    for crow in indat:
        allHg[crow[0]].append(crow[1])


# now go through and calculate min, max, mean, sd, and n for both transformed and untransformed results

# also calculate the statistics on the residuals using inverse Mills ratio for NDs
ofp = open('RRS_predicted.csv','w')
ofp.write('%s,%s,%s\n' %('ID','Realization','Hg_predicted_RRS'))

for cid in allHg:
    k = 0
    curr_mercury = np.array(allHg[cid])
    if len(curr_mercury) == 0:
        curr_mercury = np.array(list())
    else:
        curr_mercury_log = np.log((curr_mercury*1000.) + 1.0)
    if len(curr_mercury) == 0:
        pass
    else:
        for chg in curr_mercury:
            k+=1
            ofp.write('%d,%d,%f\n' %(cid,k,chg))
ofp.close()