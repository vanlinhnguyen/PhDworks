"""
Created on Aug 02 2016

@author: Linh Van Nguyen (linh.van.nguyen@hotmail.com)
"""

import numpy as np
from netCDF4 import Dataset

def data_preprocess(sspacing, tspacing):
    """
    Load coupled input-output of LR and HR from file and normalize to zero-mean
    and one- standard deviation

    Parameters
    ----------
    sspacing : 2D subsampling ratio in space (in one direction)

    tspacing : 1D subsampling ratio in time
    
    """
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
    return (Xl_tr, mea_l, sig_l, Xh_tr,mea_h,sig_h)    
    
    
####################### RIDGE REGRESSION ######################################
def RR_cv_estimate_alpha(sspacing, tspacing, alphas):
    """
    Estimate the optimal regularization parameter using grid search from a list
    and via k-fold cross validation

    Parameters
    ----------
    sspacing : 2D subsampling ratio in space (in one direction)

    tspacing : 1D subsampling ratio in time

    alphas : list of regularization parameters to do grid search
    
    """
    #Load all training data
    (Xl_tr, mea_l, sig_l, Xh_tr,mea_h,sig_h) =  data_preprocess(sspacing, tspacing)  
    
    # RidgeCV
    from sklearn.linear_model import RidgeCV    
    ridge = RidgeCV(alphas = alphas, cv = 10, fit_intercept=False, normalize=False)
    ridge.fit(Xl_tr, Xh_tr)
    
    RR_alpha_opt = ridge.alpha_
    
    print('\n Optimal lambda:', RR_alpha_opt)
    
    # save to .mat file
    import scipy.io as io
    filename = "".join(['/data/PhDworks/isotropic/regerssion/RR_cv_alpha_sspacing',
                        str(sspacing),'_tspacing',str(tspacing),'.mat'])
    io.savemat(filename, dict(alphas=alphas, RR_alpha_opt=RR_alpha_opt))
    
    # return
    return RR_alpha_opt
    
    
    
def RR_allfields(sspacing, tspacing, RR_alpha_opt):   
    """
    Reconstruct all fields using RR and save to netcdf file

    Parameters
    ----------
    sspacing : 2D subsampling ratio in space (in one direction)

    tspacing : 1D subsampling ratio in time

    RR_alpha_opt : optimal regularization parameter given from RR_cv_estimate_alpha(sspacing, tspacing, alphas)
    
    """
    # Constants
    Nh = 96
    Nt = 37
       
    #Load all training data
    (Xl_tr, mea_l, sig_l, Xh_tr,mea_h,sig_h) =  data_preprocess(sspacing, tspacing)  
    
    # Ridge Regression
    from sklearn.linear_model import Ridge
    ridge = Ridge(alpha=RR_alpha_opt, fit_intercept=False, normalize=False)
    ridge.fit(Xl_tr, Xh_tr)
    print np.shape(ridge.coef_)
       
    # Prediction and save to file
    filename = "".join(['/data/PhDworks/isotropic/regerssion/RR_sspacing',
                        str(sspacing),'_tspacing',str(tspacing),'.nc'])
    import os 
    try:
       os.remove(filename)
    except OSError:
       pass
    ncfile2 = Dataset(filename, 'w') 
    
    ncfile1 = Dataset('/data/PhDworks/isotropic/refdata_downsampled4.nc','r')
    
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
            xl = np.divide(np.reshape(xl,(1, xl.size)) - mea_l, sig_l) #pre-normalize
            xrec = np.multiply(ridge.predict(xl), sig_h) + mea_h # re-normalize the prediction
            var[t,:,:,i] = np.reshape(xrec, (Nh,Nh)) # put to netcdf file
    
    # Close file
    ncfile1.close()
    ncfile2.close()
    

####################### OTHER FUNCTIONS #######################################
def plot_learning_curve(estimator, plt, X, y, ylim=None, cv=None, n_jobs=1, 
                        train_sizes=np.linspace(.1, 1.0, 5)):
    """
    Generate a simple plot of the test and traning learning curve.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        An object of that type which is cloned for each validation.

    plt : current matplotlib plot

    X : array-like, shape (n_samples, n_features)
        Training vector, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape (n_samples) or (n_samples, n_features), optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    ylim : tuple, shape (ymin, ymax), optional
        Defines minimum and maximum yvalues plotted.

    cv : integer, cross-validation generator, optional
        If an integer is passed, it is the number of folds (defaults to 3).
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    n_jobs : integer, optional
        Number of jobs to run in parallel (default 1).
    """
    if ylim is not None:
        plt.ylim(*ylim)
    plt.xlabel("Number of training examples")
    plt.ylabel("Score")
    
    from sklearn.learning_curve import learning_curve
    train_sizes, train_scores, test_scores = learning_curve(estimator, X, y, cv=cv, 
                                                            n_jobs=n_jobs, train_sizes=train_sizes)
    
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)

    plt.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    plt.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    
    plt.grid()
    plt.legend(loc="best")
    return plt
    
    

def interp2 (x, y, z, xnew, ynew, kind='cubic'):
    from scipy import interpolate
    f = interpolate.interp2d(x, y, z, kind=kind)
    return f(xnew, ynew)

def NRMSE (xref, xrec):
    err = np.sqrt(np.sum(np.square(xref.ravel()-xrec.ravel())))/np.sqrt(np.sum(np.square(xref.ravel())))
    return err