import numpy as np
import time
import os
import pickle
import zipfile


def find_comps(SpC,Event,MASTERKEY,DL,outIndex):
    # initializations
    sttime = time.clock()
    max_iter = 15  # maximum allowable iterations
    c_iter = 0     # iteration counter starting at zero
    converged = False # flag to indicate when convergence has been achieved
    compcount = [np.inf] # a "comparable_count" variable compcount that keeps track of how many samples are comparable

    outfile = 'comps' + str(outIndex) + '.pkl'
    # determine which is the most abundant species-cut combination (SpC)
    eunich_SpC = np.unique(SpC)
    spec_counts = dict(zip(eunich_SpC,np.zeros_like(eunich_SpC)))
    for curr_SpC in eunich_SpC:
        cinds = np.nonzero(SpC == curr_SpC)
        spec_counts[curr_SpC] = len(cinds[0])
    # find the index of the most abundant SpC code
    most_abund_SpC =  max(spec_counts, key=spec_counts.get)
    
    # initialize COMPARABLE (Wente's 'TF1') which will hold all comparable samples as a vector
    COMPARABLE = np.zeros(len(MASTERKEY)).astype('int')
    COMPARABLE2  = np.zeros_like(COMPARABLE)
    # comp_Spc --> the SpC code with which all is compared.
    # initialize to most_abundant SpC code to get the party started
    comp_SpC = [most_abund_SpC]
    master_SpC_comp = []
    master_event_comp = []
    # assign comp_inds as indices of all comparable samples --> we start by assuming all fish of the most abundant
    # species-cut combination are comparable (with themselves)
    
    
    while (c_iter <= max_iter) and (converged == False):
        # ######################################## #
        # # # #  Evaluate all the SpC values # # # #
        # ######################################## #
        tmp_comp2 = np.zeros_like(COMPARABLE2)
        for curr_SpC in comp_SpC:
            if curr_SpC not in master_SpC_comp:
                master_SpC_comp.append(curr_SpC)
                comp_inds = np.nonzero(SpC == curr_SpC)[0]
                COMPARABLE[comp_inds]=1
                subsamples = Event[comp_inds]
                eunich_event = np.unique(subsamples)
                
                for curr_event in eunich_event:
                    if curr_event not in master_event_comp:
                        master_event_comp.append(curr_event)
                        # determine the indices with both SpC and event overlap AND non-censored
                        comp2_inds = np.nonzero((Event==curr_event) & (DL<1))[0] 
                        
                        COMPARABLE2[comp2_inds] = 1
                        tmp_comp2[comp2_inds] = 1
                    else:
                        continue
                comp_SpC = np.unique(SpC[np.nonzero(tmp_comp2==1)[0]])
            else:
                continue
            # use set intersection to find all non-zero indices
            COMPARABLE = COMPARABLE | COMPARABLE2
            print 'total comparable -->  {0} --> iteration {1}'.format(sum(COMPARABLE),c_iter)
    
    
        # do the convergence check and advance the iteration number
        compcount.append(sum(COMPARABLE))
        if compcount[-1]-compcount[-2] == 0:
            converged = True
        c_iter += 1
        print '\n'*2 
        print 'Iteration: ' + str(c_iter)
        print compcount
        endtime = time.clock()
        print '\n' + '*'*25
        print 'Tooooootal elapsed time = {0} minutes'.format((endtime-sttime)/60) 
        
    ofp = open(outfile,'wb')
    pickle.dump(MASTERKEY[np.nonzero(COMPARABLE)[0]],ofp)
    ofp.close()
    a = zipfile.ZipFile(outfile + '.zip','w',compression=zipfile.ZIP_DEFLATED)
    a.write(outfile)
    a.close()
    os.remove(outfile)
    endtime = time.clock()
    print '\n' + '*'*25
    print 'Tooooootal elapsed time = {0} minutes'.format((endtime-sttime)/60) 
    

def read_and_parse_input(infile):
    '''
    this code reads in and parses the input file for Fish Hg model validation.
    Note that all lists are left as string lists.

    '''
    
    time_start = time.clock()

    indat = open(infile,'r').readlines()
    headers = indat.pop(0).strip().split(',')
    
 
    DL = []
    MASTERKEY = []
    SpC = []
    Event = []

    
    
    for line in indat:
        tmp=line.strip().split(',')
        DL.append(tmp[3])
        MASTERKEY.append(tmp[6])
        SpC.append(tmp[8])
        Event.append(tmp[9])
    
    del indat
    time_end = time.clock()
    print ' elapsed time for file reading/parsing is: {0:6f} seconds'.format(time_end-time_start)
    return  DL,  MASTERKEY,  SpC, Event
