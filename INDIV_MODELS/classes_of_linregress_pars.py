import numpy as np

ols_infile = 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation_REGRESSIONS.dat'
alldat_infile= 'NatfishFinalAllobs_20111116_BABY_MODEL_trimmed_validation.csv'

class linreg_pars:
    def __init__(self,

# get the staring parameters that were generated using scipy.linregress
ols_in = np.genfromtxt(ols_infile,delimiter = ',',names=True,dtype=None)



