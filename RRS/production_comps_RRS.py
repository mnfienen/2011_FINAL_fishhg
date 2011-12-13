#!/usr/bin/python
import numpy as np
import sys
import pickle

def read_and_parse_input(infile):
    '''
    this code reads in and parses the input file for Fish Hg model validation.
    Note that all lists are left as string lists.

    '''

    indat = open(infile,'r').readlines()
    headers = indat.pop(0).strip().split(',')
 
    DL = []
    ID = []
    SpC = []
    Event = []

    # note - updated for new Access-derived data file format with ID instead of MASTERKEY
    for line in indat:
        tmp=line.strip().split(',')
        DL.append(tmp[4])
        ID.append(tmp[0])
        SpC.append(tmp[8])
        Event.append(tmp[9])
    
    del indat
    return  DL,  ID,  SpC, Event


infile = 'NatfishFinalAllobs_20111116_MNF.csv' # filename for the main input
allID = np.loadtxt('IDs.csv',skiprows=1,dtype=int)

# read in all possible left-out indices
ifp = open('all_out_inds.pkl','rb')
all_out_IDS = pickle.load(ifp)
ifp.close()

dropind_row = 0#int(sys.argv[1])
dropinds = all_out_IDS[dropind_row,:]
knock_out_ID = allID[dropinds]

ofpINDS = open('left_out_IDS.dat','w')
ofpINDS.write('left out ID codes\n')
for ci in knock_out_ID:
    ofpINDS.write('%d\n' %(ci))
ofpINDS.close()

# initializations
max_iter = 15  # maximum allowable iterations
c_iter = 0     # iteration counter starting at zero
converged = False # flag to indicate when convergence has been achieved
compcount = [np.inf] # a "comparable_count" variable compcount that keeps track of how many samples are comparable

# read and parse the input file
DL, ID, SpC, Event = read_and_parse_input(infile)

# convert to arrays the lists which need to be referenced as arrays
SpC = np.array(SpC,dtype='int')
Event = np.array(Event,dtype='int')
ID = np.array(ID,dtype='int')
DL = np.array(DL, dtype='int')


#currinds = np.nonzero(all_inds!=knock_out_ID)[0]
all_inds = np.arange(len(ID))
currinds = np.setdiff1d(all_inds,dropinds)

# drop the indices requested
SpC = SpC[currinds]
Event = Event[currinds]
ID = ID[currinds]
DL = DL[currinds]

# determine which is the most abundant species-cut combination (SpC)
eunich_SpC = np.unique(SpC)
spec_counts = dict(zip(eunich_SpC,np.zeros_like(eunich_SpC)))
for curr_SpC in eunich_SpC:
    cinds = np.nonzero(SpC == curr_SpC)
    spec_counts[curr_SpC] = len(cinds[0])
# find the index of the most abundant SpC code
most_abund_SpC =  max(spec_counts, key=spec_counts.get)

# initialize COMPARABLE (Wente's 'TF1') which will hold all comparable samples as a vector
COMPARABLE = np.zeros(len(ID)).astype('int')
COMPARABLE2  = np.zeros_like(COMPARABLE)
# comp_Spc --> the SpC code with which all is compared.
# initialize to most_abundant SpC code to get the party started
comp_SpC = [most_abund_SpC]
master_SpC_comp = []
master_event_comp = []
# assign comp_inds as indices of all comparable samples --> we start by assuming all fish of the most abundant
# species-cut combination are comparable (with themselves)

while (c_iter <= max_iter) and (converged == False):
    #print 'iteration number --> ' + str(c_iter)
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
                    # must have at least two samples meeting the comparable criterion
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
    
# del statements below are to keep memory in check
del SpC, COMPARABLE2, comp2_inds,tmp_comp2
# now check on remaining events w.r.t detection limit
compinds = np.nonzero(COMPARABLE)[0]
cEvents = Event[compinds]
# set a variable with all events
eunich_event = np.unique(Event)
del Event
cDL = DL[compinds]
del DL
del comp_inds
# find only events that have at least one value with DL==0
DL_0_Events = cEvents[cDL==0]
# set a new variable with only the DL==0 events
eunich_cEvents = np.unique(DL_0_Events)
del DL_0_Events

# find the events that havse been removed 
kill_events = np.setxor1d(eunich_cEvents,eunich_event)
for cev in kill_events:
    kill_events_inds = np.nonzero(Event==cev)[0]
    COMPARABLE[kill_events_inds] = 0 

# notcomp are the values with COMPARABLE==0
notcomp = ID[np.nonzero(COMPARABLE==0)[0]]

outfile = 'comps_%d.dat' %(dropind_row)
ofp = open(outfile,'w')
ofp.write('%d: row knocked out from all_out_inds.pkl\n' %(dropind_row))
for cv in notcomp:
    ofp.write('%d\n' %(cv))
ofp.write('\n')
ofp.close()
