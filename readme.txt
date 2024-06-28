IMPORTANT: CURRENTLY ONLY AVAILABLE FOR REVIEW BY COMPUTERS AND GEOSCIENCES
           SOFTWARE SHOULD NOT BE USED PRIOR TO PUBLICATION.
-------------------------------------------------------------------------------

Description:
-----------------------
These Matlab scripts are used to assimilate data from GRACE and GRACE-FO satellites. 
There are four data assimilation (DA) schemes available: EnKF 1D, EnKF 3D, EnKS 1D, and EnKS 3D.
EnKF and EnKS differ in how the monthly GRACE update is distributed throughout the month,
whereas the 1D and 3D are associated with the inclusion of spatial correlation errors in the update calculation. 

DA Configuration:
-----------------------
Users can customize the DA settings, which include:
   nens: number of ensembles
   para_err: error size of parameters (% of the nominal value)
   forc_err: error size of forcing data (%P,T_err,%ETpot)
   restart_err: error size of restart
   forc_err_range: magnitude of the perturbation
   pforcing_clen: correlation length for forcing data error (P,T,ETpot)
   pparameters_clen: correlation length for parameter error
   restart_clen: correlation length for restart error
   da_type: data assimilation approach (EnKF, EnKS)
   da_dim: inclusion of spatial correlation error (1D, 3D)
   grace_error_scale: scale factor of GRACE error
   grace_addmean: add model mean field to GRACE data

Functions:
-----------------------
The provided function will perform DA based on the settings. Specifically:
   HBV_lumpedx: modified HBV model
   HBV_lumped_update: HBV model with model state update
   func_perturb_parameters: generates ensemble members of a model parameter
   func_perturb_forcing: generates ensemble members of a meteorological field
   func_perturb_restart: generates ensemble members of a model initial state
   step01_load_input: reads and obtain dimensions from the input file
   step02_run_ol: performs the HBV model’s open-loop simulation
   step03_1_perturbation: applies random noises to model parameters, initial states, and forcing data
   step03_2_grace_da: assimilates GRACE/GRACE-FO data into the HBV model

Test Case:
-----------------------
The test case is included. Simply open main.m to start the DA in Matlab, then press the run or F5 key.

Setting Up Your Own DA:
-----------------------
To perform DA in other locations, configure the data array name (and size) as follows:
   P (ntime, ngrid): precipitation time series (daily)
   T (ntime, ngrid): temperature time series (daily)
   ETpot (ntime, ngrid): potential evapotranspiration time series (daily)
   forcing_date (ntime, 3): date of forcing data (daily)
   forcing_coord (ngrid, 2): coordinate of forcing data (degree)
   lat (nlat): latitude (degree) of the study domain
   lon (nlon): longitude (degree) of the study domain
   parameters (npara, ngrid): model parameters
   parameters_name	(npara, 1): name of model parameters
   grace_data (nmonth, ngrid): GRACE data (monthly)
   grace_error (nmonth, ngrid): GRACE data error (monthly)
   grace_date (nmonth, 2): Date of GRACE data (monthly)

HBV Model
-----------------------
This DA software package includes the original HBV Matlab code and ensemble parameters with catchment shapefile provided by Beck et al. (2016) (see hbv_model_beck directory).
- Beck, H.E., Dijk, A.I.J.M. van, Roo, A. de, Miralles, D.G., McVicar, T.R., Schellekens, J., Bruijnzeel, L.A., 2016. Global-scale regionalization of hydrologic model parameters. Water Resources Research 52, 3599–3622. https://doi.org/10.1002/2015WR018247

Important Note:
-----------------------
Despite the ease of use of the software, this open-source code is provided "as is," with no performance guarantees or warranties. The authors and contributors make no representations or warranties regarding the software's functionality, reliability, or performance, which may vary depending on factors such as configurations and input data. The authors and contributors are not responsible for any damages caused by using this code. Before using this code in your research, it is strongly recommended that you thoroughly validate it.

Suggested Reference:
-----------------------
Tangdamrongsub, N., Beck, H., Dong, J. (2023) A step-by-step MATLAB implementation of GRACE and GRACE-FO satellite data assimilation for enhancing terrestrial hydrology analysis. Curerntly under reviewed by Computers and Geosciences.

Author: 
-----------------------
Natthachet Tangdamrongsub @ Asian Institute of Technology (AIT), Thailand