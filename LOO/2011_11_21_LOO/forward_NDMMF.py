#! /usr/local/bin/python

'''
Forward version of the NDDMF using event and SpC parameters

Mike Fienen
mnfienen@usgs.gov

'''
import numpy as np



########################################################################
# FUNCTION TO CALCULATE MERCURY CONCENTRATION:                         #
#  the formula is: ln(C_ijk+1) = a_k + ln(length_ijk+ 1) + b_j + e_ijk #
#      where C_ijk is the ith sample of the jth event for kth SpC      #
#      a_k are the SpCpars and b_j are the Eventpars                   #
########################################################################
def calc_Hg(SpCpar,Eventpar,length):
    lnC = (SpCpar * np.log(length + 1)) + Eventpar
    return lnC

########################################################################
# FUNCTION TO CALCULATE ALL MERCURY CONCENTRATIONS:                    #
#       uses calc_Hg for all samples                                   #
########################################################################
def Hg_all(SpCpars, Eventpars, length, event, SpC, ID):
    lnC = list()         
    for i,ind in enumerate(ID):
        csp = SpC[i]
        cev = event[i]
        spcind = np.nonzero(SpCpars[:,0]==csp)[0]
        evind = np.nonzero(Eventpars[:,0]==cev)[0]
        lnCtmp = calc_Hg(SpCpars[spcind,1],Eventpars[evind,1],length[i])
        lnC.append(lnCtmp[0])
    return np.array(lnC)

