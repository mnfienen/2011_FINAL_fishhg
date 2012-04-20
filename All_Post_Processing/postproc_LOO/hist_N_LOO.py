import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42

infile = 'LOO_output_consolidated_corr.csv'

indat = np.genfromtxt(infile,delimiter = ',',names = True,dtype = None)

SpC = indat['SpC']
eunich = np.unique(SpC)

N_spc = []
for i in eunich:
    N_spc.append(len(SpC[SpC==i]))
N_spc = np.array(N_spc)
N_spc_all = N_spc.copy()
numbins = 400
len_exclude = len(N_spc[N_spc>numbins])
N_spc[N_spc>=numbins]=0
h = plt.hist(N_spc,bins = numbins+1)
plt.xlim(1,numbins)
plt.xlabel('Number of Observations in Species-cut')
plt.ylabel('Frequency of Species-cut')
plt.title('Observations in Species-cut')
tx = numbins*0.4
ty = np.max(h[0])*0.95
plt.text(tx,ty,'%d Species-Cuts have N > %d' %(len_exclude,numbins))

plt.savefig('Fig_S_2.pdf')


infile2='NatfishFinalAllobs_20111116.CSV'
indat = np.genfromtxt(infile2,delimiter = ',',names = True,dtype = None)

Event = indat['Event']
SpC_all = indat['SpC']
SpCinds = np.nonzero(np.in1d(SpC_all,SpC))[0]

Event = Event[SpCinds]
eunich = np.unique(Event)
plt.figure()
Nevent = []
for i in eunich:
    Nevent.append(len(Event[Event==i]))
Nevent = np.array(Nevent)
Nevent_all = Nevent.copy()
numbins = 60
len_exclude = len(Nevent[Nevent>numbins])
Nevent[Nevent>=numbins]=0

h = plt.hist(Nevent,bins = numbins+1)
plt.xlim(1,numbins)
plt.xlabel('Number of Observations in Event')
plt.ylabel('Frequency of Event')
plt.title('Observations in Event')
tx = numbins*0.4
ty = np.max(h[0])*0.95
plt.text(tx,ty,'%d Events have N > %d' %(len_exclude,numbins))
plt.savefig('Fig_S_3.pdf')


#plt.show()