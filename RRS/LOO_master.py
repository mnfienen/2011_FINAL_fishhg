#! /usr/local/bin/python
import numpy as np
import sys
from LOO import LOO_Hg

# N.B. this is to run with no observations dropped

####################################################################################
# Main Function
####################################################################################

Masterfile = 'NatfishFinalAllobs_20110617_MNF.csv'
Connectionsfile = 'comparable_calcs_20110630.dat'
ID_to_drop = -99999


cHg,obsHg = LOO_Hg(ID_to_drop,Masterfile,Connectionsfile)

