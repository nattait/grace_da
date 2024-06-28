function [Qsim, SM, SUZ, SLZ, SNOWPACK, ETact, MELTWATER] = HBV_lumped_update(P,ETpot,T,parameters,initial_states)
    %{ 
     Qsim = HBV_lumped(P,ETpot,T,parameters)
    [Qsim, SM, SUZ, SLZ, SNOWPACK, ETact] = HBV_lumped(P,ETpot,T,parameters)
    
    Runs the HBV hydrological model (Seibert, 2005). NaN values have to be
    removed from input arrays. The first 10 years of the record is used as
    initialization period, before the entire record is run (i.e., if the
    record is 12 years long, the first 10 years of the record are run to
    initialize the model, after which the full 12-year record is run.
    
    Input:
      P = array with daily values of precipitation (mm/d)
      ETpot = array with daily values of potential evaporation (mm/d)
      T = array with daily values of air temperature (deg C)
      parameters = array with parameter values having the following structure:
                   [BETA; CET; FC; K0; K1; K2; LP; MAXBAS; PERC; UZL; PCORR;...
                    TT; CFMAX; SFCF; CFR; CWH]
      initial_states = initial states (struct) ! Natt (added on 03 March 2021)
    
    Output:
      Qsim = array with daily values of simulated discharge (mm)
      SM (optional) = soil storage (mm)
      SUZ (optional) = upper zone storage (mm)
      SNOWPACK (optional) = snow pack depth (mm)
      ETact (optional) = actual evaporation (mm)
    
    ---------------------------
    Notes:
      Modifying from:
      Matlab implementation of HBV by Hylke Beck
      (hylke.beck@gmail.com).
    
    Last modified:
       Sep 2015: First release (Hylke)
       May 2023: Added variable updates (Natt)
    %}
    
    if nargin~=5
        error('specify correct number of input arguments')
    end

    parBETA = parameters(1,:);
    parCET = parameters(2,:);
    parFC = parameters(3,:);
    parK0 = parameters(4,:);
    parK1 = parameters(5,:);
    parK2 = parameters(6,:);
    parLP = parameters(7,:);
    parMAXBAS = parameters(8,:);
    parPERC = parameters(9,:);
    parUZL = parameters(10,:);
    parPCORR = parameters(11,:);
    parTT = parameters(12,:);
    parCFMAX = parameters(13,:);
    parSFCF = parameters(14,:);
    parCFR = parameters(15,:);
    parCWH = parameters(16,:);

    SNOWPACK = nan(size(P));  SNOWPACK=[initial_states.SNOWPACK;SNOWPACK];
    MELTWATER = nan(size(P)); MELTWATER=[initial_states.MELTWATER;MELTWATER];
    SM = nan(size(P));        SM=[initial_states.SM;SM];
    SUZ = nan(size(P));       SUZ=[initial_states.SUZ;SUZ];
    SLZ = nan(size(P));       SLZ=[initial_states.SLZ;SLZ];
    ETact = nan(size(P));     ETact=[nan(1,size(ETact,2));ETact];
    Qsim = nan(size(P));      Qsim=[nan(1,size(Qsim,2));Qsim];

    P     = [zeros(1,size(P,2));P];
    ETpot = [zeros(1,size(ETpot,2));ETpot];
    T     = [zeros(1,size(T,2));T];

    % Daily increment from analysis
    SM_Dinc        = initial_states.SM_Dinc;
    SUZ_Dinc       = initial_states.SUZ_Dinc;
    SLZ_Dinc       = initial_states.SLZ_Dinc;
    SNOWPACK_Dinc  = initial_states.SNOWPACK_Dinc;
    MELTWATER_Dinc = initial_states.MELTWATER_Dinc;
        
    for ii = length(P(1,:)):-1:1
        P(:,ii) = parPCORR(ii).*P(:,ii);
        SNOW(:,ii) = P(:,ii);
        RAIN(:,ii) = P(:,ii);
        RAIN(T(:,ii)<parTT(ii),ii) = 0;
        SNOW(T(:,ii)>=parTT(ii),ii) = 0;
        SNOW(:,ii) = SNOW(:,ii).*parSFCF(ii);
    end

    for t = 2:size(P,1)
        % Snow
        SNOWPACK(t,:) = max(SNOWPACK(t-1,:)+SNOWPACK_Dinc,0)+SNOW(t,:); % Natt: adding daily increment (if available)
        melt = parCFMAX .* (T(t,:)-parTT);
        melt(melt<0) = 0;
        melt = min(melt,SNOWPACK(t,:));
        MELTWATER(t,:) = (max(MELTWATER(t-1,:)+MELTWATER_Dinc,0))+melt; % Natt: adding daily increment (if available)
        SNOWPACK(t,:) = SNOWPACK(t,:)-melt;
        refreezing = parCFR .* parCFMAX .* (parTT-T(t,:));
        refreezing(refreezing<0) = 0;
        refreezing = min(refreezing,MELTWATER(t,:));
        SNOWPACK(t,:) = SNOWPACK(t,:)+refreezing;
        MELTWATER(t,:) = MELTWATER(t,:)-refreezing;
        tosoil = MELTWATER(t,:) - (parCWH .* SNOWPACK(t,:));
        tosoil(tosoil<0) = 0;
        MELTWATER(t,:) = MELTWATER(t,:)-tosoil;

        % Soil and evaporation
        SM(t-1,:) = max(SM(t-1,:) + SM_Dinc,0); % Natt: adding daily increment (if available)
        soil_wetness = real((SM(t-1,:)./parFC).^parBETA);
        soil_wetness(soil_wetness<0) = 0;
        soil_wetness(soil_wetness>1) = 1;
        recharge = (RAIN(t,:)+tosoil) .* soil_wetness;
        SM(t,:) = SM(t-1,:)+RAIN(t,:)+tosoil-recharge;
        excess = SM(t,:)-parFC;
        excess(excess<0) = 0;
        SM(t,:) = SM(t,:)-excess;
        evapfactor = SM(t,:) ./ (parLP .* parFC);
        evapfactor(evapfactor<0) = 0;
        evapfactor(evapfactor>1) = 1;
        ETact(t,:) = ETpot(t,:).*evapfactor;
        ETact(t,:) = min(SM(t,:), ETact(t,:));
        SM(t,:) = SM(t,:)-ETact(t,:);

        % Groundwater boxes
        SUZ(t,:) = (max(SUZ(t-1,:)+SUZ_Dinc,0))+recharge+excess; % Natt: adding daily increment (if available)
        PERC = min(SUZ(t,:), parPERC);
        SUZ(t,:) = SUZ(t,:)-PERC;
        Q0 = parK0 .* max(SUZ(t,:)-parUZL, 0.0);
        SUZ(t,:) = SUZ(t,:)-Q0;
        Q1 = parK1.*SUZ(t,:);
        SUZ(t,:) = SUZ(t,:)-Q1;
        SLZ(t,:) = (max(SLZ(t-1,:)+SLZ_Dinc,0))+PERC; % Natt: adding daily increment (if available)
        Q2 = parK2.*SLZ(t,:);
        SLZ(t,:) = SLZ(t,:)-Q2;
        Qsim(t,:) = Q0+Q1+Q2;
    end

    % Add delay to simulated runoff
    parMAXBAS = round(parMAXBAS*100)/100;
    Qsim_smoothed = zeros(size(Qsim));
    for ii = length(P(1,:)):-1:1
        window = parMAXBAS(ii)*100;
        w = zeros(round(window),1);
        for x = 1:round(window)
            w(x) = window/2 - abs(window/2-x+0.5);
        end
        w = [w; zeros(200,1)];
        w_small = zeros(ceil(parMAXBAS(ii)),1);
        for x = 1:ceil(parMAXBAS(ii))
            w_small(x) = sum(w(x*100-99:x*100));
        end
        w_small = w_small./sum(w_small);
        for bb = 1:length(w_small)
            Qsim_smoothed(bb:size(Qsim_smoothed,1),ii) = Qsim_smoothed(bb:size(Qsim_smoothed,1),ii)...
                +Qsim(1:size(Qsim_smoothed,1)-bb+1,ii)*w_small(bb);
        end
    end
    Qsim = Qsim_smoothed;

    % Remove initialization period
    Qsim(1,:) = [];
    SM(1,:) = [];
    SUZ(1,:) = [];
    SLZ(1,:) = [];
    SNOWPACK(1,:) = [];
    ETact(1,:) = [];
    MELTWATER(1,:) = [];
    
end