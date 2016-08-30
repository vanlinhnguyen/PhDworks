"""
Created on Aug 02 2016

@author: Linh Van Nguyen (linh.van.nguyen@hotmail.com)
"""
import numpy as np
#from netCDF4 import Dataset
import regressionUtils as ru

print('\n =================================================================\n') 
print('\n ======================= Ridge Regression ========================\n') 
print('\n =================================================================\n') 
n_alphas = 100
alphas = np.append(0,np.linspace(1e-10, 100, n_alphas))

# Case 5 
print('\n ========================= Start case 5 ============================') 
sspacing = 3
tspacing = 4
RR_alpha_opt = ru.RR_cv_estimate_alpha(sspacing,tspacing,alphas)

text = "".join(['\n Optimal lambda for spacing ',
                        str(sspacing),' (space) and ',str(tspacing),' (time) is ', "%.1f" % RR_alpha_opt])
print(text) 
ru.RR_allfields(sspacing, tspacing, RR_alpha_opt)

# Case 6
print('\n ========================= Start case 6 ============================') 
sspacing = 4
tspacing = 6
RR_alpha_opt = ru.RR_cv_estimate_alpha(4,6,alphas)
text = "".join(['\n Optimal lambda for spacing ',
                        str(sspacing),' (space) and ',str(tspacing),' (time) is ', "%.1f" % RR_alpha_opt])
print(text)                        
ru.RR_allfields(sspacing, tspacing, RR_alpha_opt)

# Case 7
print('\n ========================= Start case 7 ============================') 
sspacing = 6
tspacing = 8
RR_alpha_opt = ru.RR_cv_estimate_alpha(6,8,alphas)
text = "".join(['\n Optimal lambda for spacing ',
                        str(sspacing),' (space) and ',str(tspacing),' (time) is ', "%.1f" % RR_alpha_opt])
print(text)                        
ru.RR_allfields(sspacing, tspacing, RR_alpha_opt)


print('\n =================================================================\n') 
print('\n ==================== Kernel Ridge Regression ====================\n') 
print('\n =================================================================\n') 