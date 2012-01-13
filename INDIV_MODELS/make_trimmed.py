import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

infile_data = 'NatfishFinalAllobs_20111116_BABY_MODEL.csv'
infile_sum = 'NatfishFinalAllobs_20111116_BABY_SUMMARY.csv'
outfile_all = infile_data[:-4] + '_trimmed.csv'
outfile_val = infile_data[:-4] + '_trimmed_validation.csv'

indat = np.genfromtxt(infile_sum,dtype=None,names=True)
SpCE_kounts = dict(zip(indat['SpC_EVENT'],indat['Count']))
del indat
indat = open(infile_data,'r').readlines()
headers = indat.pop(0)
ofp_trim = open(outfile_all,'w')
ofp_val = open(outfile_val,'w')
ofp_trim.write(headers)
ofp_val.write(headers)
nlines = len(indat)
k = 0
for line in indat:
    k+=1
    print 'processing %d of %d' %(k,nlines)
    tmp = line.strip().split(',')
    ckount = SpCE_kounts[tmp[-1]]
    if ckount > 2:
        ofp_val.write(line)
    if ckount > 1:
        ofp_trim.write(line)

    
ofp_trim.close()
ofp_val.close()