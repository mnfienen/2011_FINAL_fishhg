import numpy as np
from forward_NDMMF import calc_Hg
import os
from scipy.stats import linregress as lm

ols_infile = 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation_REGRESSIONS.dat'
alldat_infile= 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation.csv'
main_output_file = 'LOO_Baby_Model.dat'

class all_data:
    def __init__(self,SpC,Event,length,Hg,DL,Wt,ID):
        self.SpC = SpC
        self.Event = Event
        self.length = length
        self.Hg = Hg
        self.DL = DL
        self.Wt = Wt
        self.ID = ID
        
# get the staring parameters that were generated using scipy.linregress
ols_in = np.genfromtxt(ols_infile,names=True,dtype=None)
allSpC_EVENT_indies = []
for i,cspcevent in enumerate(ols_in['SpC_EVENT']):
    if ols_in['r_squared'][i] > 0.0:
        allSpC_EVENT_indies.append(cspcevent)



alldat = np.genfromtxt(alldat_infile,delimiter = ',',dtype = None,names = True)

    
# adjust the weights by a log transformation
Wt = alldat['Wt'].astype(float)
Wt[Wt<2]=1
Wt = np.log(Wt) + 1.0
alldat['Wt']=Wt

# make a dictionary of all_data class objects containing the info
# necessary to make predictions

# remember that we are only using the ones that have r^2 != 0.0
# the 0.0 r^2 comes from either all lengths or all Hg being the same
# for a specific SpC_EVENT

allobs = []
for i,cspcevent in enumerate(allSpC_EVENT_indies):
    print i
    currinds = np.nonzero(alldat['SpC_EVENT'] == cspcevent)[0]
    allobs.append(all_data(alldat['SpC'][currinds[0]],
                           alldat['Event'][currinds[0]],
                           alldat['length'][currinds],
                           alldat['Hg'][currinds],
                           alldat['DL'][currinds],
                           alldat['Wt'][currinds],
                           alldat['ID'][currinds]))
    
    
allobs = dict(zip(allSpC_EVENT_indies,allobs))
# ########## start LOO ########## #
# now perform leave one out analaysis


# open the output file, write the header, then close it back down
ofp = open(main_output_file,'w')
ofp.write('%8s%8s%8s%17s%17s%17s%17s%10s%10s%17s%17s%17s%17s%17s%6s%6s' %('ID',
                     'SpC',
                     'Event',
                     'Hg_obs',
                     'Hg_LOO',
                     'sigma',
                     'max_log_like',
                     'num_iters',
                     'best_iter',
                     'lmSpC',
                     'lmEvent',
                     'stderrlm',
                     'NDMMF_SpC',
                     'NDMMF_Event',
                     'N',
                     'N_ND') + '\n')
ofp.close()

k = 0
ll = len(alldat)
for i,cID in enumerate(alldat['ID']):
    # here cID is the left-out ID code
    k+=1
    print 'running ID --> %d ' %(cID) + '%8d of %8d' %(k,ll)
    cspcevent = alldat['SpC_EVENT'][i]
    
    # write the NDMMF input files ,first parsing the data
    # ## calibration data
    dat_ofp = open('Hgdata.srt','w')
    numdat = 0
    
    try:
        cID in allobs[cspcevent].ID
        cDL = 0 # keep track of how many NDs there are
        len_LOO = 0
        ccallengths = []
        ccalHgs = []
        for j,IDS in enumerate(allobs[cspcevent].ID):
            if IDS != cID:
                numdat+=1
                dat_ofp.write('%d %d %f %f %d %d %d' %(allobs[cspcevent].SpC,
                                                       allobs[cspcevent].Event,
                                                       allobs[cspcevent].length[j],
                                                       allobs[cspcevent].Hg[j],
                                                       allobs[cspcevent].DL[j],
                                                       allobs[cspcevent].Wt[j],
                                                       allobs[cspcevent].ID[j]) + '\n')
                len_LOO += 1
                cDL += allobs[cspcevent].DL[j]
                ccallengths.append(allobs[cspcevent].length[j])
                ccalHgs.append(allobs[cspcevent].Hg[j])
            else:
                cHg_obs = allobs[cspcevent].Hg[j]
                clen = allobs[cspcevent].length[j]
        dat_ofp.close()
        # perform a linear regression to obtain initial parameters
        x = np.array(ccallengths)
        x = np.log(x+1.0)
        y = np.array(ccalHgs)
        y = np.log((y*1000.0)+1)
        # check to be sure of uniqueness in both length and Hg
        if ((len(x) - len(np.unique((x)))) >= 1):
            if ((len(y) - len(np.unique((y)))) >= 1):
            
                cspclm, ceventlm, r_value, p_value, cstderrlm = lm(x,y)
                sigma_calc = np.std(y)
                # ## sigma value, using STD of Hg
                sig_ofp = open('Hgsigma.dat','w')
                sig_ofp.write('%f\n' %(sigma_calc))
                sig_ofp.close()
                # ## SpC parameter value
                spc_ofp = open('Hgspc.srt','w')
                spc_ofp.write('%d %f\n' %(allobs[cspcevent].SpC,cspclm))
                spc_ofp.close()
                # ## Event parameter value
                event_ofp = open('Hgevents.srt','w')
                event_ofp.write('%d %f\n' %(allobs[cspcevent].Event,ceventlm))
                event_ofp.close()
                # ## Write the index file
                ndx_ofp = open('Hgdata.ndx','w')
                for i in np.arange(numdat):
                    ndx_ofp.write('%d\n' % (i+1))
                ndx_ofp.close()
            
                if ((len_LOO - cDL) > 1):
                    # call the external C-code Newton-Raphson parameter estimation code
                    os.system('./NRparest > nul')    
                    
                    # finally, read in the results and make the Hg prediction for the left-out value
                    # SpC parameters
                    SpCpars = np.loadtxt('BestSPs')
                    # Event parameters
                    Eventpars = np.loadtxt('BestEPs')
                    
                    # calculate mercury for this index
                    cHg = calc_Hg(SpCpars[1],Eventpars[1],clen)
                    cHg = (np.exp(cHg)-1)/1000.0
                    
                    # read in sigma from the summary data file
                    set1 = open('summaryRESULTS.dat','r').readlines()
                    tmp = set1[1].strip().split()
                    # tmp is [max_sig  max_loglike  total_iters best_iteration]
                    ofp = open(main_output_file,'a')
                    ofp.write('%8d%8d%8d%17.8e%17.8e%17.8e%17.8e%10d%10d%17.8e%17.8e%17.8e%17.8e%17.8e%6d%6d' %(cID,
                                                               allobs[cspcevent].SpC,
                                                               allobs[cspcevent].Event,
                                                               cHg_obs,
                                                               cHg,
                                                               float(tmp[0]),
                                                               float(tmp[1]),
                                                               float(tmp[2]),
                                                               float(tmp[3]),
                                                               cspclm,
                                                               ceventlm,
                                                               sigma_calc,
                                                               SpCpars[1],
                                                               Eventpars[1],
                                                               len(y),
                                                               cDL) + '\n')
                    ofp.close()
    except KeyError:
        continue