# driver for reading, parsing, and plotting RRS data

import numpy as np
import os

from summary_plotting_RRS import RRS_post, killer

# user-specified variables
infile = 'RRS_residuals_all.csv'
clean_slate = True

# logging information
if os.path.exists(os.path.join(os.getcwd(),'figures')) == False:
    os.mkdir(os.path.join(os.getcwd(),'figures'))
logfile = open('figures/LOO_postprocALLONLY.log','w')
logfile.write('Postprocessing of RRS analysis\n\n\n')

# initialize a class for output
RRS = RRS_post(infile,logfile)


# load the data
RRS.load_data()

# parse the data
RRS.parse_input()

# make the barcharts
#RRS.barcharts()

# make the boxplots for SpC codes
'''
RRS.boxplotSpC(abs_flag=False)
RRS.boxplotSpC(abs_flag=True)
'''

# make the plots by data source
RRS.boxplot_DataSource()
