import numpy as np

infile = 'NatfishFinalAllobs_20111116_MNF.csv'
outfile = 'NatfishFinalAllobs_20111116_BABY_MODEL.csv'
indat = open(infile,'r').readlines()

header = indat.pop(0)
ofp = open(outfile,'w')
header = header.strip().split(',')
for cv in header:
    ofp.write('%s,' %(cv))
ofp.write('SpC_EVENT\n')
SpC_EVENT = []
for line in indat:
    outline = line.strip()
    ofp.write(outline)
    tmp = outline.split(',')
    SpC_EVENT.append(tmp[6] + '_' + tmp[7])
    ofp.write(SpC_EVENT[-1] + '\n')
ofp.close()

# now count up the number of occurrences in each SpC_EVENT
SpC_EVENT = np.array(SpC_EVENT)
eunich_SpCE = np.unique(SpC_EVENT)
outfile_summary = 'NatfishFinalAllobs_20111116_BABY_SUMMARY'
ofp = open(outfile_summary,'w')
ofp.write('%20s%20s\n' %('SpC_EVENT','Count'))
for i in eunich_SpCE:
    cinds = np.nonzero(SpC_EVENT == i)[0]
    ofp.write('%20s%20d\n' %(i,len(cinds)))
ofp.close()
