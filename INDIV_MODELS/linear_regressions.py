import numpy as np
from scipy.stats import linregress as lm

infile = 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation.csv'
outfile = 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation_REGRESSIONS.dat'
# read in the trimmed data set
indat = np.genfromtxt(infile,delimiter = ',', dtype=None,names=True)
Hg = indat['Hg']
lens = indat['length']
DL = ['DL']
SpC_Event = indat['SpC_EVENT']

# set NDs as hald detection limit for linear regressions
NDs = np.nonzero(DL==1)[0]
for cind in NDs:
    Hg[cind] = 0.5*Hg[cind]

# perform the log transformations
Hg = np.log((Hg*1000.0)+1)
lens = np.log(lens+1.0)

# now perform all the linear regressions, writing the results to a file
ofp = open(outfile,'w')
ofp.write('%20s'*6 %('SpC_EVENT','SpC_par','Event_par','sigma','r_squared','N') + '\n')
allSpC_Events = np.unique(SpC_Event)
k=0
allK = len(allSpC_Events)
for cspcev in allSpC_Events:
    k+=1
    print "rockin' " + cspcev + '-->%d of %d' %(k,allK)
    cind = np.nonzero(SpC_Event == cspcev)[0]
    y = Hg[cind]
    x = lens[cind]
    slope, intercept, r_value, p_value, std_err = lm(x,y)
    sigma_calc = np.std(y)
    ofp.write('%19s %19f %19f %19f %19f %19d\n' %(cspcev,slope,intercept,sigma_calc,r_value**2,len(x)))
ofp.close()