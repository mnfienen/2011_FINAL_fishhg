import numpy as np
import sys



# set number of realizations
num_rel = 1000
allHg_file = 'allHg.pkl'
# determine whether to read and parse the Hg data, or pull it from a pickled instance
parse_all_flag = 1

# parse all the results files

# first read in the master ID list
ID_master = np.loadtxt('../../RRS/RRS_RESULTS_12_15_2011/all_IDs.dat',skiprows = 1, dtype=int)


# read in the reference values fo Hg, corrected for NDs with Inverse Mills Ratio
indat_ref = np.genfromtxt('AllHg_reference.dat',names=True,dtype=None)
ID_ref = indat_ref['ID']
Hg_obs_ref = indat_ref['Hgobscorr']

# initialize a dictionary to hold the results
allHg = dict([(x,list()) for x in ID_master])

for crel in xrange(num_rel):
    print 'Running realization %d of %d:' %(crel+1,num_rel)
    indat = np.genfromtxt('../../RRS/RRS_RESULTS_12_15_2011/Hg_%d.dat' %(crel),dtype=None,names=True,skiprows=1)
    for crow in indat:
        allHg[crow[0]].append(crow[1])


# now go through and calculate min, max, mean, sd, and n for both transformed and untransformed results
'''
ofpHg = open('Hgstats_ppm.dat','w')
ofpHg.write('%12s%12s%12s%12s%12s%12s\n' %('ID','N','min','mean','max','SD'))

ofpHg_log = open('Hgstats_log.dat','w')
ofpHg_log.write('%12s%12s%12s%12s%12s%12s\n' %('ID','N','min','mean','max','SD'))
'''
# also calculate the statistics on the residuals using inverse Mills ratio for NDs
ofpHg_log_res = open('Hgstats_log_res.dat','w')
ofpHg_log_res.write('%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s\n' %('ID','N','min_abs_res','min_res','mean_res_abs','mean_res','max_abs_res','max_res','SD_res','SD_res_abs'))

for cid in allHg:
    curr_mercury = np.array(allHg[cid])
    if len(curr_mercury) == 0:
        curr_mercury = np.array(list())
    else:
        curr_mercury_log = np.log((curr_mercury*1000.) + 1.0)
    if len(curr_mercury) == 0:
        #ofpHg.write('%12d%12d%12d%12d%12d%12d\n' %(cid,0,-999,-999,-999,-999))
        #ofpHg_log.write('%12d%12d%12d%12d%12d%12d\n' %(cid,0,-999,-999,-999,-999))
        #ofpHg_log_res.write('%12d%12d%12d%12d%12d%12d%12d%12d\n' %(cid,0,-999,-999,-999,-999,-999,-999))
        pass
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
        #ofpHg.write('%12d%12d%12.5f%12.5f%12.5f%12.5f\n' %(cid,nhg,minhg,meanhg,maxhg,sdhg))
        #ofpHg_log.write('%12d%12d%12.5f%12.5f%12.5f%12.5f\n' %(cid,nhg,minhg_log,meanhg_log,maxhg_log,sdhg))
        # now handles the residuals
        hg_obs_curr = Hg_obs_ref[np.nonzero(ID_ref==cid)[0][0]]
        # log transform
        hg_obs_curr = np.log((hg_obs_curr * 1000.0) + 1.0)
        hg_res_curr = (hg_obs_curr-curr_mercury_log)
        hg_res_curr_abs = np.abs(hg_res_curr)
        min_hg_res_abs = np.min(hg_res_curr_abs)
        max_hg_res_abs = np.max(hg_res_curr_abs)
        min_hg_res = np.min(hg_res_curr)
        max_hg_res = np.max(hg_res_curr)
        mean_hg_res = np.mean(hg_res_curr)
        mean_hg_res_abs = np.mean(hg_res_curr_abs)
        sdhg_res = np.std(hg_res_curr)
        sdhg_res_abs = np.std(hg_res_curr_abs)

        
        ofpHg_log_res.write('%16d%16d%16.5f%16.5f%16.5f%16.5f%16.5f%16.5f%16.5f%16.5f\n' 
                            %(cid,
                              nhg,
                              min_hg_res_abs,
                              min_hg_res,
                              mean_hg_res_abs,
                              mean_hg_res,
                              max_hg_res_abs,
                              max_hg_res,
                              sdhg_res,
                              sdhg_res_abs))
        
ofpHg.close()
ofpHg_log.close()
ofpHg_log_res.close()
