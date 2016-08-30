# PhD works - isotropic turbulence dataset

This repos synthesizes all codes and data to reproduce our results in the thesis: 

*Reconstruction of finely resolved velocity fields in turbulent flows from low resolution measurements*

## Main parts

The repos is comprised by the following main part

* Regression: linear and nonlinear regression; require [scikit-learn](http://scikit-learn.org/stable/)
* Dictionary Learning: statistics of different learning approaches, couple dictionary learning by three post-processing techniques; require [SPAMS]{http://spams-devel.gforge.inria.fr/}
* NLM: Matlab and C code with openMP to speedup the computation
* Bayesian fusion: code to run all cases, with step to estimate the weights
* comparison_plot: to compare all methods via NRMSE and plots
* datastats: statistics of the data

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/linhvannguyen/PhDworks/LICENSE.md) file for details
