# -*- coding: utf-8 -*-
"""
Created on Thu May  5 11:00:38 2016

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
        xl = xh[0:-1:sspacing,0:-1:sspacing] # xh[np.meshgrid(HTLS_sknots,HTLS_sknots)]
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

############## Ridge Regression ###############################################
#from sklearn.linear_model import Ridge
import scipy.io as sio
mf = sio.loadmat('/data/ISOTROPIC/regression/RR_cv_alpha_sspacing6_tspacing8.mat', 
                 squeeze_me=True, struct_as_record=False)
RR_alpha_opt = mf['RR_lambda_opt']
print'Optimal lambda:', RR_alpha_opt

from sklearn.linear_model import Ridge
ridge = Ridge(alpha=RR_alpha_opt, fit_intercept=False, normalize=False)
ridge.fit(Xl_tr, Xh_tr)
print np.shape(ridge.coef_)

# save to file
#import pickle
#with open('/data/ISOTROPIC/regression/RRmodel_sspacing3_tspacing4.pickle', 'w') as f:
#    pickle.dump([ridge], f)

# Getting back the objects:
#with open('/data/ISOTROPIC/regression/RRmodel_sspacing3_tspacing4.pickle') as f:
#    ridge = pickle.load(f)


############## Prediction and save to file ####################################
import os 
try:
   os.remove('/data/ISOTROPIC/data/RR_sspacing6_tspacing8.nc')
except OSError:
   pass
ncfile2 = Dataset('/data/ISOTROPIC/data/RR_sspacing6_tspacing8.nc', 'w') 

ncfile1 = Dataset('/data/ISOTROPIC/data/data_downsampled4.nc','r')

# create the dimensions
ncfile2.createDimension('Nt',Nt)
ncfile2.createDimension('Nz',Nh)
ncfile2.createDimension('Ny',Nh)
ncfile2.createDimension('Nx',Nh)
# create the var and its attribute
var = ncfile2.createVariable('Urec', 'd',('Nt','Nz','Ny','Nx'))

for t in range(Nt):
    print('3D snapshot:',t)
    for i in range(Nh):
        xl = np.array(ncfile1.variables['velocity_x'][t,0:Nh:sspacing,0:Nh:sspacing,i]) # load only LR
        xl = np.divide(np.reshape(xl,(1, Nl*Nl)) - mea_l, sig_l) #pre-normalize
        xrec = np.multiply(ridge.predict(xl), sig_h) + mea_h # re-normalize the prediction
        var[t,:,:,i] = np.reshape(xrec, (Nh,Nh)) # put to netcdf file

# Close file
ncfile1.close()
ncfile2.close()