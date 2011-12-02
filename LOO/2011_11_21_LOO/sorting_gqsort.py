'''
Sorting program to sort HGdata.dat into HGdata.srt

This emulates the C code written by David Donato calls gqsort

sorting using np.lexsort found on Stack Overflow:
http://stackoverflow.com/questions/2706605/sorting-a-2d-numpy-array-by-multiple-axes
'''

import numpy as np
def gqsort_Hgdat(infile,outfile):
    
    tmp = np.loadtxt(infile)
    
    # sort by columns in the order 1,0,4 (note apparent reversal of order)
    ind = np.lexsort((tmp[:,4],
                      tmp[:,0],
                      tmp[:,1]))
    
    tmp = tmp[ind]
    
    ofp = open(outfile,'w')
    for line in tmp:
        ofp.write('%3d %7d %13.8f %13.8f %2d %13.8f\n'
                  %(line[0],
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                    line[5]))
     
    ofp.close()
def gqsort_Hgdatext(infile,outfile):
    
    tmp = np.loadtxt(infile)
    
    # sort by columns in the order 1,0,4 (note apparent reversal of order)
    ind = np.lexsort((tmp[:,4],
                      tmp[:,1],
                      tmp[:,0]))
    
    tmp = tmp[ind]
    
    ofp = open(outfile,'w')
    for line in tmp:
        ofp.write('%3d %7d %13.8f %13.8f %2d  %13.8f %8d\n'
                  %(line[0],
                    line[1],
                    line[2],
                    line[3],
                    line[4],
                    line[5],
                    line[6]))
     
    ofp.close()