import numpy as np
import comps_utils as cu

# initializations 
infile = 'NatfishFinalAllobs_20111110_MNF.CSV' # filename for the main input
#infile = 'NatfishFinalAllobs_20110617_MNF.CSV'
knockout = 1616



[outVAL,notcomp,eltime, ip, c_iter]=cu.find_comps_all_local(infile,knockout)

print 'knocked out --> ' + str(knockout)
for i in notcomp:
    print i
    
print '-'*5
print outVAL