# driver for reading, parsing, and plotting LOO data

import numpy as np
import os

from LOO_postproc import LOO_post, killer

# user-specified variables
infile = 'LOO_output_consolidated_corr.csv'
clean_slate = True

# logging information
if os.path.exists(os.path.join(os.getcwd(),'figures')) == False:
    os.mkdir(os.path.join(os.getcwd(),'figures'))
logfile = open('figures/LOO_postprocALLONLY.log','w')
logfile.write('Postprocessing of LOO analysis\n\n\n')

# initialize a class for output
LOO = LOO_post(infile,logFlag = True)

# load the data
LOO.load_data()

# parse the data
LOO.parse_input()

# make misfit plots - one plot with all values
curr_plots = 'figures/consolidated'
if clean_slate:
    killer(curr_plots)
LOO.plot_consolidated_misfit(logfile,subf = curr_plots)

maincutoff = 20
# make individual species-specific misfit plots
curr_plots = 'figures/species'
if clean_slate:
    killer(curr_plots)
LOO.plot_each_series_misfit(logfile, cutoff = maincutoff,series = 'CommonName',subf = curr_plots)

# make individual species-specific misfit plots
curr_plots = 'figures/species_cut'
if clean_slate:
    killer(curr_plots)
LOO.plot_each_series_misfit(logfile, cutoff = maincutoff,series = 'SpC',subf = curr_plots)

# make individual species-specific misfit plots
curr_plots = 'figures/data_source'
if clean_slate:
    killer(curr_plots)
LOO.plot_each_series_misfit(logfile, cutoff = maincutoff,series = 'DataSource',subf = curr_plots)

logfile.close()

# make boxplots for the top 5 species, grouped by abundance
LOO.boxplots(toprank = 5)
