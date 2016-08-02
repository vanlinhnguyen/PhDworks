"""
Created on Aug 02 2016

@author: Linh Van Nguyen (linh.van.nguyen@hotmail.com)
"""

import numpy as np
from netCDF4 import Dataset

def RR_cv_estimate_alpha(sspacing, tspacing, alphas):
    # Constants
    Nh = 96
    Nt = 37
    
    # Position of measurements in space-time
    HTLS_sknots = np.arange(0,Nh,sspacing)
    LTHS_tknots = np.arange(0,Nh,tspacing)
    Nl = len(HTLS_sknots)
    Ns = len(LTHS_tknots)
    
    # Dimension of HTLS and LTHS
    P = Nh*Nh
    Q = Nl*Nl
    M = Nt*Ns
    
    #Load all training data
    Xh_tr = np.zeros((M, P))
    Xl_tr = np.zeros((M, Q))
    ncfile1 = Dataset('/data/ISOTROPIC/data/data_downsampled4.nc','r')
    for t in range(Nt):
        count = 0
        for i in LTHS_tknots:
            xh = np.array(ncfile1.variables['velocity_x'][t,0:Nh,0:Nh,i])
            xl = xh[0:-1:sspacing,0:-1:sspacing] # xh[np.meshgrid(HTLS_sknots,HTLS_sknots)]
            Xh_tr[t*Ns + count,:] = np.reshape(xh,(1, P))
            Xl_tr[t*Ns + count,:] = np.reshape(xl,(1, Q))
            count = count + 1
    ncfile1.close()
    
    
    # normalized: centered, variance 1
    mea_l = np.zeros(Q)
    sig_l = np.zeros(Q)
    for k in range(Q):
        mea_l[k] = np.mean(Xl_tr[:,k])
        sig_l[k] = np.std(Xl_tr[:,k])
        Xl_tr[:,k] = (Xl_tr[:,k]-mea_l[k])/sig_l[k] 
        
    mea_h = np.zeros(P)
    sig_h = np.zeros(P)
    for k in range(P):
        mea_h[k] = np.mean(Xh_tr[:,k])
        sig_h[k] = np.std(Xh_tr[:,k])
        Xh_tr[:,k] = (Xh_tr[:,k]-mea_h[k])/sig_h[k]     
    
    # RidgeCV
    from sklearn.linear_model import RidgeCV    
    ridge = RidgeCV(alphas = alphas, cv = 10, fit_intercept=False, normalize=False)
    ridge.fit(Xl_tr, Xh_tr)
    
    RR_lambda_opt = ridge.alpha_
    
    print('\n Optimal lambda:', RR_lambda_opt)
    
    # save to .mat file
    import scipy.io as io
    filename = "".join(['/data/PhDworks/isotropic/regerssion/RR_cv_alpha_sspacing',
                        str(sspacing),'_tspacing',str(tspacing),'.mat'])
    io.savemat(filename, dict(alphas=alphas, RR_lambda_opt=RR_lambda_opt))
    
    # return
    return RR_lambda_opt