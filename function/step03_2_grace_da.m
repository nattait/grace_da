function DA = step03_2_grace_da(data,configure)
    %{
    grace_da assimilates GRACE/GRACE-FO data into the HBV model
    Syntax:  
		DA = step03_2_grace_da(data,configure)
    Inputs:
		data - Essential inputs
			 - Required parameters
				  - coordinate:     Coordinates (degree) of the study domain (longitude, latitude)
				  - forcing_date:   Calendar date associated with the forcing data
				  - OL_SM:          Soil moisture estimates from the openloop simulation
				  - OL_SLZ:         Lower zone storage estimates from the openloop simulation
				  - OL_SUZ:         Upper zone storage estimates from the openloop simulation
				  - P:              Precipitation field
				  - T:              Temperature field
				  - ETpot:          Potential evapotranspiration field
				  - SM_init:        Initial state of SM
				  - SUZ_init:       Initial state of SUZ
				  - SLZ_init:       Initial state of SLZ
				  - SNOWPACK_init:  Initial state of SNOWPACK
				  - MELTWATER_init: Initial state of MELTWATER
				  - parameters:     Ensemble parameters
				  - grace_data:     GRACE/GRACE-FO data
                  - grace_error:    Error of GRACE/GRACE-FO data
				  - grace_date:     Calendar date associated with GRACE/GRACE-FO data
		configure - DA Configurations/settings
                  - Required parameters
                      - nens:       Number of ensemble members to generate (e.g., 100)
                      - da_dim:     DA scheme in spatial domain
                                    options: '1D','3D'
					  - da_type:    DA scheme in temporal domain
                                    options: 'EnKF','EnKS'
                      - grace_error_scale: Scale factor of GRACE error used in 1D scheme (e.g., 1.2)
                      - grace_addmean: Option to add TWS Long-term mean from model to GRACE
                                       1 = Add TWS long-term mean computed from model to GRACE
                                       0 = Remove long-term mean value from model TWS components
    Output:
       DA - GRACE DA results contaning model variables with dimension (ntime,ngrid,2). 
            The third dimension store mean (1) and standard deviation (2) values.
    Examples: 
       DA = step03_2_grace_da(data,configure);   
    
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
    %}

    %------------- BEGIN CODE --------------
    nstate_upd = 3; % Update SM, SUZ, and SLZ
    ngrid = size(data.coordinate,1);
    ntime_DA = size(data.grace_date,1);
    nens = configure.nens; % Number of ensemble members

    %%%% Declare arrays to save DA results
    DA.Qsim      = zeros(ntime_DA,ngrid,2);
    DA.TWS       = zeros(ntime_DA,ngrid,2);
    DA.SM        = zeros(ntime_DA,ngrid,2);
    DA.SUZ       = zeros(ntime_DA,ngrid,2);
    DA.SLZ       = zeros(ntime_DA,ngrid,2);
    DA.SNOWPACK  = zeros(ntime_DA,ngrid,2);
    DA.ETact     = zeros(ntime_DA,ngrid,2);           
    DA.MELTWATER = zeros(ntime_DA,ngrid,2);
    DA.date = zeros(ntime_DA,3);
       
    %%%% Compute long-term mean of storage from OL with respect to the assimilation period
    SM_mean  = mean(data.OL_SM,1);
    SLZ_mean = mean(data.OL_SLZ,1);
    SUZ_mean = mean(data.OL_SUZ,1);

    %%%% Assimilae monthly TWS data into model
    tstop = 0;
    for run_id=1:size(data.grace_date,1) 
        %%%% Search for simulation period consistent with GRACE data period
        simp      = find(data.forcing_date(:,1)==data.grace_date(run_id,1) & ...
                    data.forcing_date(:,2)==data.grace_date(run_id,2)); 
        P_sim     = data.P(simp(1):simp(end),:,:);
        T_sim     = data.T(simp(1):simp(end),:,:);
        ETpot_sim = data.ETpot(simp(1):simp(end),:,:);
        ntime_sim = size(P_sim,1);
        DA.date(run_id,:) = data.forcing_date(simp(end),:);

        %%%%% Declare arrays for saving monthly assimilation results
        Qsim_sim      = zeros(ntime_sim,ngrid,nens);
        SM_sim        = zeros(ntime_sim,ngrid,nens);
        SUZ_sim       = zeros(ntime_sim,ngrid,nens);
        SLZ_sim       = zeros(ntime_sim,ngrid,nens);
        SNOWPACK_sim  = zeros(ntime_sim,ngrid,nens);
        ETact_sim     = zeros(ntime_sim,ngrid,nens);
        MELTWATER_sim = zeros(ntime_sim,ngrid,nens);

        %%%%% Forecast step (model simulation without DA)
        for ens = 1:nens
            %%%%% Define initial states
            if run_id == 1 % Use model states from OL for the first month
                initial_states.SM        = reshape(data.SM_init(:,ens),1,ngrid);
                initial_states.SUZ       = reshape(data.SUZ_init(:,ens),1,ngrid);
                initial_states.SLZ       = reshape(data.SLZ_init(:,ens),1,ngrid);
                initial_states.SNOWPACK  = reshape(data.SNOWPACK_init(:,ens),1,ngrid);
                initial_states.MELTWATER = reshape(data.MELTWATER_init(:,ens),1,ngrid);
            else %%%%% Use model states from the last step of previous simulation period
                initial_states.SM        = reshape(SM_sim_state(:,ens),1,ngrid);
                initial_states.SUZ       = reshape(SUZ_sim_state(:,ens),1,ngrid);
                initial_states.SLZ       = reshape(SLZ_sim_state(:,ens),1,ngrid);
                initial_states.SNOWPACK  = reshape(SNOWPACK_sim_state(:,ens),1,ngrid);
                initial_states.MELTWATER = reshape(MELTWATER_sim_state(:,ens),1,ngrid);
            end

            %%%%% Setting daily increment to zeros in forecast step
            initial_states.SM_Dinc        = zeros(1,ngrid);
            initial_states.SUZ_Dinc       = zeros(1,ngrid);
            initial_states.SLZ_Dinc       = zeros(1,ngrid);
            initial_states.SNOWPACK_Dinc  = zeros(1,ngrid);
            initial_states.MELTWATER_Dinc = zeros(1,ngrid);

            [Qsim_sim(:,:,ens), SM_sim(:,:,ens), SUZ_sim(:,:,ens), SLZ_sim(:,:,ens), ...
                SNOWPACK_sim(:,:,ens), ETact_sim(:,:,ens), MELTWATER_sim(:,:,ens)] = ...
                HBV_lumped_update(P_sim(:,:,ens),ETpot_sim(:,:,ens),T_sim(:,:,ens), ...
                data.parameters(:,:,ens), initial_states);
        end

        %%%%% Setting DA record time
        tstart = tstop + 1;
        tstop  = (tstart - 1) + length(simp);

        %%%%% Check whether GRACE observation is available
        obs  = data.grace_data(run_id,:);

        %%%%% Perfroming GRACE DA when the observation is available (analysis step)
        if ~all(isnan(obs))
            %%%% Use grace_error_scale only for test, set grace_error_scale = 1 in real case
            obs_err = data.grace_error(run_id,:) * configure.grace_error_scale;

            if strcmp(configure.da_dim,'3D')
                %%%%% Perturb observation
                Y      = obs_err(1)*randn(1,nens); % Generate observation error
                if configure.grace_addmean == 1
                    Dmat   = obs(1) + mean(SM_mean+SLZ_mean+SUZ_mean) + Y;
                else
                    Dmat   = obs(1) + Y;
                end

                %%%% Collect model states for entire region
                %%%% Note: only SM,SUZ,SLZ are considered here since snow is not presented in the study area
                Upd    = zeros(nstate_upd,ngrid,nens);
                Aamat  = zeros(nstate_upd,ngrid,nens);
                Amat = zeros(nstate_upd*ngrid,nens);
                Hmat = (1/ngrid) * ones(1,nstate_upd*ngrid);
                ji = 0;
                for gi = 1:ngrid
                    if configure.grace_addmean == 1
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SM_sim(:,gi,:)));
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SUZ_sim(:,gi,:)));
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SLZ_sim(:,gi,:)));
                    else
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SM_sim(:,gi,:)))  - SM_mean(gi);  % monthly SM anomaly
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SUZ_sim(:,gi,:))) - SLZ_mean(gi); % monthly SUZ anomaly
                        ji = ji+1; Amat(ji,:) = squeeze(mean(SLZ_sim(:,gi,:))) - SUZ_mean(gi); % monthly SLZ anomaly
                    end
                end
                Dp = Dmat - (Hmat * Amat); % (1,nens)
                Ap = Amat - repmat(mean(Amat,2),1,size(Amat,2)); % (nstate_upd,nens)
                Pe = Ap * Ap';
                Re = Y * Y';
                K  = (Pe * Hmat')/(Hmat * (Pe * Hmat') + Re);
                iUpd = K * Dp;

                ji=0;
                for gi = 1:ngrid
                    if configure.grace_addmean == 1
                        ji = ji+1; Upd(1,gi,:) = iUpd(ji,:); Aamat(1,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:); % monthly update SM
                        ji = ji+1; Upd(2,gi,:) = iUpd(ji,:); Aamat(2,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:); % monthly update SUZ
                        ji = ji+1; Upd(3,gi,:) = iUpd(ji,:); Aamat(3,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:); % monthly update SLZ
                    else
                        ji = ji+1; Upd(1,gi,:) = iUpd(ji,:); Aamat(1,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:) + SM_mean(gi); % monthly update SM
                        ji = ji+1; Upd(2,gi,:) = iUpd(ji,:); Aamat(2,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:) + SLZ_mean(gi); % monthly update SUZ
                        ji = ji+1; Upd(3,gi,:) = iUpd(ji,:); Aamat(3,gi,:) = squeeze(Amat(ji,:)) + iUpd(ji,:) + SUZ_mean(gi); % monthly update SLZ
                    end
                end

            elseif strcmp(configure.da_dim,'1D')
                %%%% Perturbing observation
                Y      = obs_err'.*randn(ngrid,nens);
                if configure.grace_addmean == 1
                    Dmat   = repmat((obs+SM_mean+SLZ_mean+SUZ_mean)',1,nens) + Y;
                else
                    Dmat   = repmat(obs',1,nens) + Y;
                end

                %%%%% Perform DA 1D (grid-by-grid basis)
                Hmat  = ones(1,nstate_upd); % Measurement operator
                Upd   = zeros(nstate_upd,ngrid,nens);
                Aamat = zeros(nstate_upd,ngrid,nens);
                for gi = 1:ngrid
                    %%%%% Collecting model states
                    Amat = zeros(nstate_upd,nens);
                    if configure.grace_addmean == 1
                        Amat(1,:) = squeeze(mean(SM_sim(:,gi,:)));
                        Amat(2,:) = squeeze(mean(SUZ_sim(:,gi,:)));
                        Amat(3,:) = squeeze(mean(SLZ_sim(:,gi,:)));
                    else
                        Amat(1,:) = squeeze(mean(SM_sim(:,gi,:)))  - SM_mean(gi);  % monthly SM anomaly
                        Amat(2,:) = squeeze(mean(SUZ_sim(:,gi,:))) - SLZ_mean(gi); % monthly SUZ anomaly
                        Amat(3,:) = squeeze(mean(SLZ_sim(:,gi,:))) - SUZ_mean(gi); % monthly SLZ anomaly
                    end

                    Dp    = Dmat(gi,:) - (Hmat * Amat); % (1,nens)
                    Ap    = Amat - repmat(mean(Amat,2),1,size(Amat,2)); % (nstate_upd,nens)
                    Pe    = Ap * Ap'; % (nstate_upd,nstate_upd)
                    Re    = Y(gi,:) * Y(gi,:)'; % (1,1)
                    Kmat  = (Pe * Hmat')/(Hmat * (Pe * Hmat') + Re); % kalman gain

                    %%%%% Collecting update
                    Upd(:,gi,:)   = Kmat * Dp;
                    if configure.grace_addmean == 1
                        Aamat(1,gi,:) = Amat(1,:) + squeeze(Upd(1,gi,:))';
                        Aamat(2,gi,:) = Amat(2,:) + squeeze(Upd(2,gi,:))';
                        Aamat(3,gi,:) = Amat(3,:) + squeeze(Upd(3,gi,:))';
                    else
                        Aamat(1,gi,:) = Amat(1,:) + squeeze(Upd(1,gi,:))' + SM_mean(gi);
                        Aamat(2,gi,:) = Amat(2,:) + squeeze(Upd(2,gi,:))' + SLZ_mean(gi);
                        Aamat(3,gi,:) = Amat(3,:) + squeeze(Upd(3,gi,:))' + SUZ_mean(gi);
                    end
                end
            end

            if strcmp(configure.da_type,'EnKS')
                %%%%% Compute daily increment (nstate_upd, ngrid, nens)
                Upd_Dinc = Upd/ntime_sim;

                %%%%% Repeat Run (with adding daily increment)
                for ens = 1:nens
                    %%%%% Reinitialize initial states
                    if run_id == 1 % Use model states from OL for the first month
                        initial_states.SM        = reshape(data.SM_init(:,ens),1,ngrid);
                        initial_states.SUZ       = reshape(data.SUZ_init(:,ens),1,ngrid);
                        initial_states.SLZ       = reshape(data.SLZ_init(:,ens),1,ngrid);
                        initial_states.SNOWPACK  = reshape(data.SNOWPACK_init(:,ens),1,ngrid);
                        initial_states.MELTWATER = reshape(data.MELTWATER_init(:,ens),1,ngrid);
                    else %%%%% Use model states from the last step of previous simulation period
                        initial_states.SM        = reshape(SM_sim_state(:,ens),1,ngrid);
                        initial_states.SUZ       = reshape(SUZ_sim_state(:,ens),1,ngrid);
                        initial_states.SLZ       = reshape(SLZ_sim_state(:,ens),1,ngrid);
                        initial_states.SNOWPACK  = reshape(SNOWPACK_sim_state(:,ens),1,ngrid);
                        initial_states.MELTWATER = reshape(MELTWATER_sim_state(:,ens),1,ngrid);
                    end

                    %%%%% Setting daily increment
                    initial_states.SM_Dinc        = reshape(Upd_Dinc(1,:,ens),1,ngrid);
                    initial_states.SUZ_Dinc       = reshape(Upd_Dinc(2,:,ens),1,ngrid);
                    initial_states.SLZ_Dinc       = reshape(Upd_Dinc(3,:,ens),1,ngrid);
                    initial_states.SNOWPACK_Dinc  = zeros(1,ngrid);
                    initial_states.MELTWATER_Dinc = zeros(1,ngrid);

                    %%%%% Run model with adding daily increment
                    [Qsim_sim(:,:,ens), SM_sim(:,:,ens), SUZ_sim(:,:,ens), SLZ_sim(:,:,ens), ...
                        SNOWPACK_sim(:,:,ens), ETact_sim(:,:,ens), MELTWATER_sim(:,:,ens)] = ...
                        HBV_lumped_update(P_sim(:,:,ens),ETpot_sim(:,:,ens),T_sim(:,:,ens), ...
                        data.parameters(:,:,ens), initial_states);
                end

                %%%%% Update initial states for next simulation period,
                %%%%% using the last day of the repeated run
                SM_sim_state        = squeeze(SM_sim(end,:,:)); % (ngrid,nens)
                SUZ_sim_state       = squeeze(SUZ_sim(end,:,:));
                SLZ_sim_state       = squeeze(SLZ_sim(end,:,:));
                SNOWPACK_sim_state  = squeeze(SNOWPACK_sim(end,:,:));
                MELTWATER_sim_state = squeeze(MELTWATER_sim(end,:,:));

            elseif strcmp(configure.da_type,'EnKF')
                %%%%% Apply the monthly update to the last day of the month
                %%%%% No repeat run is required
                SM_sim_state        = squeeze(Aamat(1,:,:)); % (ngrid,nens)
                SUZ_sim_state       = squeeze(Aamat(2,:,:));
                SLZ_sim_state       = squeeze(Aamat(3,:,:));
                SNOWPACK_sim_state  = squeeze(SNOWPACK_sim(end,:,:));
                MELTWATER_sim_state = squeeze(MELTWATER_sim(end,:,:));
            end
        else
            %%%%% When observations are not available
            %%%%% Update initial states for next simulation period
            %%%%% using the last day of the forecast run
            SM_sim_state        = squeeze(SM_sim(end,:,:)); % (ngrid,nens)
            SUZ_sim_state       = squeeze(SUZ_sim(end,:,:));
            SLZ_sim_state       = squeeze(SLZ_sim(end,:,:));
            SNOWPACK_sim_state  = squeeze(SNOWPACK_sim(end,:,:));
            MELTWATER_sim_state = squeeze(MELTWATER_sim(end,:,:));
        end

        %%%%% Save DA results
        TWS_sim = SM_sim + SUZ_sim + SLZ_sim;
        DA.Qsim(tstart:tstop,:,1)      = mean(Qsim_sim,3);      DA.Qsim(tstart:tstop,:,2)      = std(Qsim_sim,0,3);
        DA.TWS(tstart:tstop,:,1)       = mean(TWS_sim,3);       DA.TWS(tstart:tstop,:,2)       = std(TWS_sim,0,3);
        DA.SM(tstart:tstop,:,1)        = mean(SM_sim,3);        DA.SM(tstart:tstop,:,2)        = std(SM_sim,0,3);
        DA.SUZ(tstart:tstop,:,1)       = mean(SUZ_sim,3);       DA.SUZ(tstart:tstop,:,2)       = std(SUZ_sim,0,3);
        DA.SLZ(tstart:tstop,:,1)       = mean(SLZ_sim,3);       DA.SLZ(tstart:tstop,:,2)       = std(SLZ_sim,0,3);
        DA.SNOWPACK(tstart:tstop,:,1)  = mean(SNOWPACK_sim,3);  DA.SNOWPACK(tstart:tstop,:,2)  = std(SNOWPACK_sim,0,3);
        DA.ETact(tstart:tstop,:,1)     = mean(ETact_sim,3);     DA.ETact(tstart:tstop,:,2)     = std(ETact_sim,0,3);
        DA.MELTWATER(tstart:tstop,:,1) = mean(MELTWATER_sim,3); DA.MELTWATER(tstart:tstop,:,2) = std(MELTWATER_sim,0,3);
    end
	%------------- END OF CODE --------------
end
