function PERTURB = step03_1_perturbation(OL,input,configure)
    %{
    perturbation applies random noises to model parameters, initial states and forcing data
    Syntax:  
		PERTURB = step03_1_perturbation(OL,input,configure)
    Inputs:
		OL        - A variable contained model outputs from the openloop (e.g., from step02_run_ol)
        input     - A variable contained all essential input data (e.g., from step01_load_input)
        configure - Pertubation setting
                  - Required parameters
					- nens             - Number of ensemble members (e.g., 100)
                    - para_err         - Error size of parameters (percent of nominal value) (e.g., 10)
                    - forc_err         - Error size of forcing data (P (%), T (not %), ETpot (%)) (e.g., [5 2 5])
                    - forc_err_range   - Magnitude of the perturbation
                                         options: 'median','minmax','nominal'
                    - pforcing_clen    - Correlation length for forcing correlation error in degree (P,T,ETpot) (e.g., [1.0 1.0 1.0])
                    - pparameters_clen - Correlation length for parameters correlation error (e.g., 1.0)
                    - restart_clen     - Correlation length for initial state's correlation error in degree (e.g., 1.0)
                    - restart_err      - Error size of initial states (percent of nominal value) (e.g., 10)
    Output:
		PERTURB - A variable contained ensemble fields (i.e., forcing, parameters, initial states)
    Example: 
		PERTURB = step03_1_perturbation(OL,input,configure);

    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
	%}
	
	%------------- BEGIN CODE --------------
    % Perturb model parameters
    PERTURB.para_ens = repmat(input.parameters,1,1,configure.nens);
    PERTURB.para_ens(1,:,:)  = func_perturb_parameters(input.parameters(1,:),input,configure); % BETA
    PERTURB.para_ens(3,:,:)  = func_perturb_parameters(input.parameters(3,:),input,configure); % FC
    PERTURB.para_ens(4,:,:)  = func_perturb_parameters(input.parameters(4,:),input,configure); % K0
    PERTURB.para_ens(5,:,:)  = func_perturb_parameters(input.parameters(5,:),input,configure); % K1
    PERTURB.para_ens(6,:,:)  = func_perturb_parameters(input.parameters(6,:),input,configure); % K2
    PERTURB.para_ens(7,:,:)  = func_perturb_parameters(input.parameters(7,:),input,configure); % LP
    PERTURB.para_ens(9,:,:)  = func_perturb_parameters(input.parameters(9,:),input,configure); % PERC

    % Perturb forcing data
    PERTURB.P_ens     = func_perturb_forcing(input,configure,'P');
    PERTURB.T_ens     = func_perturb_forcing(input,configure,'T');
    PERTURB.ETpot_ens = func_perturb_forcing(input,configure,'ETpot');

    PERTURB.P_ens(PERTURB.P_ens < 0) = 0;
    PERTURB.ETpot_ens(PERTURB.ETpot_ens < 0) = 0;

    % Ensemble runs to obtain initial states
    sim_start = find(input.forcing_date(:,1)==input.grace_date(1,1) & ...
                          input.forcing_date(:,2)==input.grace_date(1,2));

    PERTURB.Qsim_init      = func_perturb_restart(OL.Qsim(sim_start(1),:),input,configure);
    PERTURB.SM_init        = func_perturb_restart(OL.SM(sim_start(1),:),input,configure);
    PERTURB.SUZ_init       = func_perturb_restart(OL.SUZ(sim_start(1),:),input,configure);
    PERTURB.SLZ_init       = func_perturb_restart(OL.SLZ(sim_start(1),:),input,configure);
    PERTURB.SNOWPACK_init  = func_perturb_restart(OL.SNOWPACK(sim_start(1),:),input,configure);
    PERTURB.ETact_init     = func_perturb_restart(OL.ETact(sim_start(1),:),input,configure);
    PERTURB.MELTWATER_init = func_perturb_restart(OL.MELTWATER(sim_start(1),:),input,configure);
    %------------- END OF CODE --------------
end