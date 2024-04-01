# Readme for codes accompanying "Technical Note: Altitude scaling of Cl-36 production from Fe"

# Introduction 
This package contains codes and accompanying data that were used to calibrate production rates and to make the figures in the manuscript “Technical Note: Altitude scaling of Cl-36 production from Fe”. These codes have been tested in Octave version 7.3 and require the “io” and “optim” packages.
 
# File Index:
## Folder – Calibration
This folder contains the codes used to calibrate production rates. 

_Scaling factors.m_ – Code used to calculate time-averaged, reaction-specific scaling factors. 

### subfolder – attenuation
_Argento_attenuation.m_ – Fits an exponential function to subsurface production polynomials from Argento et al. (2015) for Cl-36 from K and Be-10 in quartz in granite and calculates the attenuation length ratio.

_Fe_attenuation.m_ – Fits an exponential function to the Be-10 in quartz and Cl-36 from Fe production profiles in the atmosphere calculated using LSDscaling_mod and calculates the attenuation length ratio.  

### subfolder – constant exposure
_constant_exp_FevBe.m_ – Calibrates the production rate of Cl-36 from Fe against Be-10 in quartz using a constant-exposure geomorphic model that considers differential decay. 

_constant_exp_FevK.m_ – Calibrates the production rate of Cl-36 from Fe against Cl-36 in feldspar using a constant-exposure geomorphic model.

_constant_exp_K.m_ – Calibrates the production rate of Cl-36 from neutron capture against Be-10 in quartz using a constant-exposure geomorphic model that considers differential decay.

_Owens_FevBe.m_ – Runs constant_exp_FevBe for the Owens Valley samples and propagates error using a finite difference approximation to the partial derivatives. 

_Owens_FevK.m_ – Runs constant_exp_FevK for the Owens Valley samples and propagates error using a finite difference approximation to the partial derivatives. 

_Owens_K.m_ – Runs constant_exp_K for the Owens Valley samples and propagates error using a finite difference approximation to the partial derivatives. 

### subfolder – steady erosion
#### sub-subfolder muons
_muon_parameters_ – function for setting-up muon parameters using a 3-exponential term model

_comppars36_muons_ – calculates muon production parameters for Cl-36 (code modified from Marrero et al. (2016) to use a 3-exponential term muon model)

_comppars1026_muons_ – calculates muon production parameters for Be-10 (code modified from Marrero et al. (2016) to use a 3-exponential term muon model)

#### sub-subfolder – neutrons
_Mt_Evans.m_ – Runs steady_E_cal for the Mt. Evans sample and propagates error using a finite-difference approximation to the partial derivatives. 

_steady_E_cal.m_ – Calibrates the production rates of Cl-36 from Fe against Cl-36 in feldspar and Be-10 in quartz and of Cl-36 from K against Be-10 in quartz for the Mt. Evans sample using a steady-state erosion geomorphic model.

## Folder – Data
This folder contains spreadsheets with input data used for the production rate calibrations and in making Figures 4 and 5.  

_Mt Evans.xlsx_ – compositional data for determining low-energy neutron production of 36Cl in the sample from Mt. Evans

_Input data.xlsx_ – measured nuclide concentrations and other input parameters used for the production rate calibrations and making Figures 4 and 5 

## Folder – Figures
This folder contains the codes used to generate the figures in the main text:

_Banana.m_ – Makes Figure 5, a two-isotope Cl-36/Be-10 plot for the Mt. Evans sample showing that this sample falls on or near the steady-state erosion line. 

_Final_Plot_ – Makes Figure 4, a plot of calibrated scaling factor ratios against the modeled scaling factor ratios with altitude. 

_Fluxes_Plot_ – Makes Figure 1, a plot of excitation functions and normalized cosmic ray fluxes.

_Nuclide_Altitude_Plot_ – Makes Figure 2, a plot of the deviation of production ratios from the sea level ratios with increasing altitude. 

## Folder – Production Modeling
This folder contains codes that are used to model production rates and scaling factors. 

_air_pressure.m_ – Calculates the air pressure at an input location using the ERA40 model.   

_avg_RC.m_ - Makes a linearly spaced vector of cutoff rigidity (Rc) and solar modulation potential (SPhi) between the present and an input time for a specified location using the models described in Lifton et al. (2014).

_LSDscaling_mod.m_ – Calculates reaction-specific scaling factors.   

_Neutrons_mod.m_ - Calculates theoretical spallation production rates from neutrons.

_Protons_mod.m_ - Calculates theoretical spallation production rates from protons. 

### subfolder - Shared files
This folder contains files from Lifton et al. (2014) and modeled reference sea level and high latitude production rates (ref.mat & Reference.mat). 

### subfolder - Muon scaling
This subfolder contains codes for modeling muon production.

_Interpolate.m_ - interpolation function from Marrero et al. (2016)

_Muonfluxsato.m_ - Calculates production rates from muons using code from Marrero et al. (2016).

_Muon_model.m_ - Calculates site-specific coefficients and attenuation lengths for a single-exponential term model of muon depth-attenuation that can be used in the low-energy neutron flux model. 

_Muon_model_muons.m_ - Calculates site-specific coefficients and attenuation lengths for a three-exponential term model of muon depth-attenuation (two terms for slow negative muons and one for fast muons), that can be used to model muon production of Be-10 in quartz and Cl-36 from K and Ca in feldspar with greater accuracy than the single-term model. 

_MuonsX.m_ - muon model from Lifton et al., (2014) as included in Marrero et al. (2016) 

# References:
Argento, D.C., Stone, J.O., Reedy, R.C., O’Brien, K., 2015. Physics-based modeling of cosmogenic nuclides part II – Key aspects of in-situ cosmogenic nuclide production. Quaternary Geochronology, The CRONUS-EARTH Volume: Part I 26, 44–55. https://doi.org/10.1016/j.quageo.2014.09.005

Lifton, N., Sato, T., Dunai, T.J., 2014. Scaling in situ cosmogenic nuclide production rates using analytical approximations to atmospheric cosmic-ray fluxes. Earth and Planetary Science Letters 386, 149–160. https://doi.org/10.1016/j.epsl.2013.10.052

Marrero, S.M., Phillips, F.M., Borchers, B., Lifton, N., Aumer, R., Balco, G., 2016. Cosmogenic nuclide systematics and the CRONUScalc program. Quaternary Geochronology 31, 160–187.



