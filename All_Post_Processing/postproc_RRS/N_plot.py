import numpy as np
import matplotlib.pyplot as plt

infile = 'Hgstats_log.dat'
indat = np.genfromtxt(infile,dtype=None,names=True)

plt.figure
plt.bar(indat['ID'],indat['N'])
plt.savefig('N_plot.pdf')