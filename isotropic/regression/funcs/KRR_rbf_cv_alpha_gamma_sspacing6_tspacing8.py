# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:32:11 2016

@author: nguyen
"""

import numpy as np
from netCDF4 import Dataset

# Constants
Nh = 96
Nt = 37
sspacing = 6
tspacing = 8

HTLS_sknots = np.arange(0,Nh,sspacing)
HTHS_sknots = np.arange(0,Nh,1)
LTHS_tknots = np.arange(0,Nh,tspacing)
Nl = len(HTLS_sknots)
Ns = len(LTHS_tknots)

Dh = Nh*Nh
Dl = Nl*Nl

N = Nt*Ns

#Load all training data
Xh_tr = np.zeros((N, Dh))
Xl_tr = np.zeros((N, Dl))
ncfile1 = Dataset('/data/ISOTROPIC/data/data_downsampled4.nc','r')
for t in range(Nt):
    count = 0
    for i in LTHS_tknots:
        xh = np.array(ncfile1.variables['velocity_x'][t,0:Nh,0:Nh,i])
        xl = xh[0:-1:sspacing,0:-1:sspacing]
        Xh_tr[t*Ns + count,:] = np.reshape(xh,(1, Dh))
        Xl_tr[t*Ns + count,:] = np.reshape(xl,(1, Dl))
        count = count + 1
ncfile1.close()


# normalized: centered, variance 1
mea_l = np.zeros(Dl)
sig_l = np.zeros(Dl)
for k in range(Dl):
    mea_l[k] = np.mean(Xl_tr[:,k])
    sig_l[k] = np.std(Xl_tr[:,k])
    Xl_tr[:,k] = (Xl_tr[:,k]-mea_l[k])/sig_l[k] 
    
mea_h = np.zeros(Dh)
sig_h = np.zeros(Dh)
for k in range(Dh):
    mea_h[k] = np.mean(Xh_tr[:,k])
    sig_h[k] = np.std(Xh_tr[:,k])
    Xh_tr[:,k] = (Xh_tr[:,k]-mea_h[k])/sig_h[k]     
    
# Use GridSearchCV to get a rough idea of parameters
from sklearn.grid_search import GridSearchCV
from sklearn.kernel_ridge import KernelRidge

alphas = np.logspace(-5, -1, 10)
gammas = np.logspace(-7, 3, 10)
kr = GridSearchCV(KernelRidge(kernel='rbf'), cv=6,
                  param_grid={"alpha": alphas, "gamma": gammas})
kr.fit(Xl_tr, Xh_tr)

print('\n Best estimator:', kr.best_estimator_)

# save to .mat file
import scipy.io as io
io.savemat('/data/ISOTROPIC/regression/KRR_rbf_cv_alpha_gamma_sspacing6_tspacing8.mat', 
           dict(alphas=alphas, gammas=gammas, KRR_alpha_opt=kr.best_estimator_.alpha,
                KRR_gamma_opt=kr.best_estimator_.gamma))
                