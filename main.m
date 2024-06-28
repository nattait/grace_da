%{
    Description:  
	    GRACE data assimilation process with example settings
    Inputs:
	    data\input.mat - An example forcing data and model parameters
    Required functions
        step01_load_input     - Reads and obtain dimension observation from the input file
        step02_run_ol         - Performs the HBV model's openloop simulation
        step03_1_perturbation - Applies random noises to model parameters, initial states and forcing data
        step03_2_grace_da     - Assimilates GRACE/GRACE-FO data into the HBV model
    Output:
        DA - Data assimilation results
    Examples: 
        In Terminal:   Type "matlab main.m"
        In Matlab GUI: Click Run or F5
    
    -----------------------------
    Description:
        Refer to function/description.m
    -----------------------------
%}

%------------- BEGIN CODE --------------

addpath('function')

filename = 'data_test\input.mat';

%% Step 1: Prepare model inputs

input = step01_load_input(filename); % Load input

%% Step 2: Perform an open-loop run (OL)

OL = step02_run_ol(input); % rename from the above

%% Step 3.1: Apply perturbation

% Configure how model and input are perturbed
pert_configure.nens              = 100;
pert_configure.para_err          = 10;            % error size of parameters (percent of nominal value)
pert_configure.forc_err          = [5 2 5];       % error size of forcing data (%P,T_err,%ETpot)
pert_configure.forc_err_range    = 'minmax';
pert_configure.pforcing_clen     = [1.0 1.0 1.0]; % correlation length for forcing correlation error (P,T,ETpot)
pert_configure.pparameters_clen  = 1.0;           % correlation length for parameters correlation error
pert_configure.restart_clen      = 1.0;
pert_configure.restart_err       = 10;

PERTURB = step03_1_perturbation(OL,input,pert_configure);

%% Step 3.2: GRACE DA

% Prepare data for DA
data.P               = PERTURB.P_ens;
data.T               = PERTURB.T_ens;
data.ETpot           = PERTURB.ETpot_ens;
data.Qsim_init       = PERTURB.Qsim_init;
data.SM_init         = PERTURB.SM_init;
data.SUZ_init        = PERTURB.SUZ_init;
data.SLZ_init        = PERTURB.SLZ_init;
data.SNOWPACK_init   = PERTURB.SNOWPACK_init;
data.ETact_init      = PERTURB.ETact_init;
data.MELTWATER_init  = PERTURB.MELTWATER_init;
data.parameters      = PERTURB.para_ens;  

data.OL_Qsim         = OL.Qsim;
data.OL_SM           = OL.SM;
data.OL_SUZ          = OL.SUZ;
data.OL_SLZ          = OL.SLZ;
data.OL_SNOWPACK     = OL.SNOWPACK;
data.OL_ETact        = OL.ETact;
data.OL_MELTWATER    = OL.MELTWATER;

data.parameters_name = input.parameters_name;
data.coordinate      = input.forcing_coord;
data.forcing_date    = input.forcing_date;
data.grace_date      = input.grace_date;
data.grace_data      = input.grace_data;
data.grace_error     = input.grace_error;

% Configure DA parameters
configure.nens              = 100;
configure.da_type           = 'EnKS';
configure.da_dim            = '3D';
configure.grace_error_scale = 1; % only for test, grace_error_scale = 1 in real world case
configure.grace_addmean     = 1;

% Run GRACE DA
DA = step03_2_grace_da(data,configure);

%------------- END OF CODE --------------