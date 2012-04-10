# driver for reading, parsing, and plotting LOO data

import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages
from LOO_postproc import LOO_post, killer

# user-specified variables
infile = 'LOO_for_plotting_4_2_2012.csv'
clean_slate = True
maincutoff = 30
# logging information
if os.path.exists(os.path.join(os.getcwd(),'figures')) == False:
    os.mkdir(os.path.join(os.getcwd(),'figures'))
logfile = open('figures/LOO_postprocSpC.log','w')
logfile.write('Postprocessing of LOO analysis for species-cut combos\n\n\n')

# initialize a class for output
LOO = LOO_post(infile,logFlag = True)

# load the data
LOO.load_data()

# parse the data
LOO.parse_input()


# make individual species_cut_combo-specific misfit plots
curr_plots = 'figures/species_cut'
if clean_slate:
    killer(curr_plots)
LOO.plot_each_series_misfit(logfile, cutoff = maincutoff,series = 'SpC',subf = curr_plots)

logfile.close()
