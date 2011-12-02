
import numpy as np

infile = 'output_LOO.dat'

indat = open(infile,'r').readlines()
indat = np.array(indat)
'''
modmeas = indat[4::5]
droppedIDs = indat[::5]
sigs = indat[2::5]
'''
allIDS = []
allSIG = []
allMEAS = []
allMOD = []
missing_modeled = []

ofp = open('LOO_output_consolidated.csv','w')
ofp.write('ID,measured,LOO_modeled,sigma\n')

for i,line in enumerate(indat):
    tmpline = line.strip().split()
    if tmpline[0] == 'dropped_ID':
        tmp_MM = indat[i+1].strip().split()
        try:            
            allIDS.append(int(tmp_MM[0]))
            allMEAS.append(float(tmp_MM[2]))
            allMOD.append(float(tmp_MM[1]))
            allSIG.append(float(indat[i-1].strip().split()[0]))
            ofp.write('%d,%f,%f,%f\n' %(allIDS[-1],allMEAS[-1],allMOD[-1],allSIG[-1]))
        except:
            missing_modeled.append(int(indat[i-3].strip().split()[1]))

ofp.close()            