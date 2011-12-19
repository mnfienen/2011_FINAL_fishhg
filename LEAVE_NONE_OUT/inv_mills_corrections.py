'''
Code to correct the residuals for values that include censored values
using the inverse Mills ratio

a m!ke@usgs joint
mnfienen@usgs.gov
'''
import numpy as np
from scipy.special import erf

def Inverse_Mills_Ratio_Correction(alphaIM,sigmaIM,y):
    
    # calcualte the Inverse Mills Ratio as a function of detection limit (alphaIM)
    # and standard error of the regression (sigmaIM) 
    # y is the modeled value
    # the output takes advantage of the relation that x = E(res) + y
    #sigmaIM = (np.exp(sigmaIM)-1.0)

    alphaIM = np.log((alphaIM*1000.0)+1.0)
    y = np.log((y*1000.0)+1.0)
    
    # note that, here, we are really after the expected value of the residual!!!!!!!
    alphaIM = alphaIM - y
    # first calcualte the PDF
    PDF = 1.0/(sigmaIM*np.sqrt(2*np.pi))* np.exp(-(1.0/2.0)*((alphaIM)**2/(sigmaIM**2)))
    
    # next calculate the CDF
    CDF = (1.0/2.0)*(1+(erf(alphaIM/(alphaIM*np.sqrt(2)))))
    
    # calculate the expected value of the residual, defined as sigmaIM*(-PDF/CDF)
    E_res = sigmaIM * (-PDF/CDF)
    x = y + E_res
    x = (np.exp(x)-1.0)/1000.0
    return x,E_res

