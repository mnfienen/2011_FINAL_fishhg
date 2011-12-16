import numpy as np
import sys



# set number of realizations
num_rel = 1000
allHg_file = 'allHg.pkl'
# determine whether to read and parse the Hg data, or pull it from a pickled instance
parse_all_flag = 1

# parse all the results files

# first read in the master ID list
ID_master = np.loadtxt('../RRS_RESULTS_12_15_2011/all_IDs.dat',skiprows = 1, dtype=int)


# initialize a dictionary to hold the results
allHg = dict([(x,list()) for x in ID_master])

for crel in xrange(num_rel):
    print 'Running realization %d of %d:' %(crel+1,num_rel)
    indat = np.genfromtxt('../RRS_RESULTS_12_15_2011/Hg_%d.dat' %(crel),dtype=None,names=True,skiprows=1)
    for crow in indat:
        allHg[crow[0]].append(crow[1])


# now go through and calculate min, max, mean, sd, and n for both transformed and untransformed results
ofpHg = open('Hgstats_ppm.dat','w')
ofpHg.write('%12s%12s%12s%12s%12s%12s\n' %('ID','N','min','mean','max','SD'))

ofpHg_log = open('Hgstats_log.dat','w')
ofpHg_log.write('%12s%12s%12s%12s%12s%12s\n' %('ID','N','min','mean','max','SD'))

for cid in allHg:
    curr_mercury = np.array(allHg[cid])
    if len(curr_mercury) == 0:
        curr_mercury = np.array(list())
    else:
        curr_mercury_log = np.log((curr_mercury*1000.) + 1.0)
    if len(curr_mercury) == 0:
        ofpHg.write('%12d%12d%12d%12d%12d%12d\n' %(cid,0,-999,-999,-999,-999))
        ofpHg_log.write('%12d%12d%12d%12d%12d%12d\n' %(cid,0,-999,-999,-999,-999))
    else:
        nhg   = len(curr_mercury)
        minhg = np.min(curr_mercury)
        maxhg = np.max(curr_mercury)
        meanhg = np.mean(curr_mercury)
        sdhg = np.std(curr_mercury)
        minhg_log = np.min(curr_mercury_log)
        maxhg_log = np.max(curr_mercury_log)
        meanhg_log = np.mean(curr_mercury_log)
        sdhg_log = np.std(curr_mercury_log)
        ofpHg.write('%12d%12d%12.5f%12.5f%12.5f%12.5f\n' %(cid,nhg,minhg,meanhg,maxhg,sdhg))
        ofpHg_log.write('%12d%12d%12.5f%12.5f%12.5f%12.5f\n' %(cid,nhg,minhg_log,meanhg_log,maxhg_log,sdhg))

ofpHg.close()
ofpHg_log.close()