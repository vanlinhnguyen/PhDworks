"""
Created on Aug 03 2016

@author: Linh Van Nguyen (linh.van.nguyen@hotmail.com)
"""
import numpy as np
#from netCDF4 import Dataset
import regressionUtils as ru

# sspacing = 8, tspacing = 6
n_alphas = 100
alphas = np.append(0,np.linspace(1e-10, 100, n_alphas))
RR_alpha_opt = ru.RR_cv_estimate_alpha(6,8,alphas)
print('\n Optimal lambda:', RR_alpha_opt)