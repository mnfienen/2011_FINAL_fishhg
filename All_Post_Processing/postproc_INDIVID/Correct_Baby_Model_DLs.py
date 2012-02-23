import numpy as np
from inv_mills_corrections import Inverse_Mills_Ratio_Correction

# load all the data
infile = 'LOO_Baby_Model_w_DL.dat'

indat = np.genfromtxt(infile,names=True,dtype=None)
allsigma = indat['sigma']
allHgmod = indat['Hg_LOO']
allHgobs = indat['Hg_obs']
allDL = indat['DL']

inlines = open(infile,'r').readlines()
header = inlines.pop(0).strip() + ' Hg_obs_corrected Hg_mod_baby\n'

ofp = open(infile + 'corr.dat','w')
ofp.write(header)

numlines = len(indat)

for i,cline in enumerate(inlines):
    print 'processing %d of %d' %(i,numlines)
    if allDL[i] == 1:
        HgCorr,junkus = Inverse_Mills_Ratio_Correction(allHgobs[i],allsigma[i],allHgmod[i])
    else:
        HgCorr = allHgobs[i]
    ofp.write(cline.strip() + ' %f %f \n' %(HgCorr, allHgmod[i]))
ofp.close()