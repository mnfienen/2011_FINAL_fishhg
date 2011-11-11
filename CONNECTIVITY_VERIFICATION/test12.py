import numpy as np
import comps_utils as cu
import os, time, socket


# initializations 
infile = 'NatfishFinalAllobs_20111110_MNF.CSV' # filename for the main input
#infile = 'NatfishFinalAllobs_20110617_MNF.CSV'
cind = 2653


# initializations
sttime = time.time()
max_iter = 15  # maximum allowable iterations
c_iter = 0     # iteration counter starting at zero
converged = False # flag to indicate when convergence has been achieved
compcount = [np.inf] # a "comparable_count" variable compcount that keeps track of how many samples are comparable

# read and parse the input file

DL, ID, SpC, Event = cu.read_and_parse_input(infile)

# convert to arrays the lists which need to be referenced as arrays
SpC = np.array(SpC,dtype='int')
Event = np.array(Event,dtype='int')
ID = np.array(ID,dtype='int')
DL = np.array(DL, dtype='int')

total_length = len(ID)
all_inds = np.copy(ID)
lenallinds = len(all_inds)
currinds = np.nonzero(all_inds!=cind)[0]
SpC = SpC[currinds]
Event = Event[currinds]
ID = ID[currinds]
cc = ID[ID==currID]
DL = DL[currinds]
lenDL = len(DL)
outVAL = cind
os.exit()
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
    print 'iteration number --> ' + str(c_iter)
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
                   # else:
                    #    COMPARABLE[comp2_inds] = 0
                    #MNF DEBUG/
                    if (curr_SpC in [168,240,332]):
                        print str(curr_event)
                        print str(curr_SpC)
                        jjj = 1
                        currIDS = ID[comp2_inds]
                        currDLS = DL[comp2_inds]
                        #/MNF DEBUG


                     #MNF DEBUG/
                    print 'currSpC = ' + str(curr_SpC)
                    print 'curr_event = ' + str(curr_event)
                    
                    cc = COMPARABLE[cind-5:cind+1]
                    print cc
                    cc2 = COMPARABLE2[cind-5:cind+1]
                    print cc2
                    idid = ID[cind-5:cind+1]
                    print idid
                    #/MNF DEBUG
                else:
                    continue
                  
            comp_SpC = np.unique(SpC[np.nonzero(tmp_comp2==1)[0]])
        else:
            continue
        # use set intersection to find all non-zero indices
        aaaaa = sum(COMPARABLE)
        bbbbb = sum(COMPARABLE2)
        COMPARABLE = COMPARABLE | COMPARABLE2
        ccccccc= sum(COMPARABLE)
        print 'total comparable -->  {0} --> iteration {1}'.format(sum(COMPARABLE),c_iter)
     

    # do the convergence check and advance the iteration number
    compcount.append(sum(COMPARABLE))
    if compcount[-1]-compcount[-2] == 0:
        converged = True
    c_iter += 1
    
COMPARABLE = COMPARABLE & COMPARABLE2
cc = COMPARABLE[1614:1617]
cc2 = COMPARABLE2[1614:1617]
idid = ID[1614:1617]  
notcomptmp = ID[np.nonzero(COMPARABLE==0)[0]]
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
# find only events that have at leasts one value with DL==0
DL_0_Events = cEvents[cDL==0]
# set a new variable with only 
eunich_cEvents = np.unique(DL_0_Events)
del DL_0_Events
# find the events that havse been removed 
kill_events = set(eunich_cEvents) ^ set(eunich_event)
for cev in kill_events:
    kill_events_inds = np.nonzero(cEvents==cev)[0]
    COMPARABLE[kill_events_inds] = 0 

# return only the IDs for which COMPARABLE==0
notcomp = ID[np.nonzero(COMPARABLE==0)[0]]
eltime = time.time() - sttime
try:
    ip = socket.gethostbyname(socket.gethostname())
except:
    ip = 'cannot read IP'