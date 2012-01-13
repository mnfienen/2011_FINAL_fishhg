import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

infile = 'NatfishFinalAllobs_20111116_BABY_SUMMARY.csv'

indat = np.genfromtxt(infile,dtype=None,names=True)
plt.figure()
n,bins,patches =  plt.hist(indat['Count'],bins=np.linspace(0,np.max(indat['Count']+1),np.max(indat['Count'])+2))
plt.xlim((0,50))
print '%d' %(np.max(indat['Count']))
bins = bins.astype(int)

hist_results = dict(zip(bins,n))

ofp = open('binned_summary.dat','w')
ofp.write('%20s%20s\n' %('No. in Spc_Event','Frequency'))
for i in hist_results:
    ofp.write('%20d%20d\n' %(i,hist_results[i]))
ofp.close()
plt.savefig('histogram.pdf')