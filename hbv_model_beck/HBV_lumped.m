function [Qsim, SM, SUZ, SLZ, SNOWPACK, ETact] = HBV_lumped(P,ETpot,T,parameters)
% Qsim = HBV_lumped(P,ETpot,T,parameters)
% [Qsim, SM, SUZ, SLZ, SNOWPACK, ETact] = HBV_lumped(P,ETpot,T,parameters)
%
% Runs the HBV hydrological model (Seibert, 2005). NaN values have to be
% removed from input arrays. The first 10 years of the record is used as
% initialization period, before the entire record is run (i.e., if the
% record is 12 years long, the first 10 years of the record are run to
% initialize the model, after which the full 12-year record is run.
%
% Input:
%   P = array with daily values of precipitation (mm/d)
%   ETpot = array with daily values of potential evaporation (mm/d)
%   T = array with daily values of air temperature (deg C)
%   parameters = array with parameter values having the following structure:
%                [BETA; CET; FC; K0; K1; K2; LP; MAXBAS; PERC; UZL; PCORR;...
%                 TT; CFMAX; SFCF; CFR; CWH]
%
% Output:
%   Qsim = array with daily values of simulated discharge (mm)
%   SM (optional) = soil storage (mm)
%   SUZ (optional) = upper zone storage (mm)
%   SNOWPACK (optional) = snow pack depth (mm)
%   ETact (optional) = actual evaporation (mm)
%
%   Matlab implementation of HBV by Hylke Beck
%   (hylke.beck@gmail.com).
%
%   Last modified 13 September 2015

if nargin~=4
    error('specify correct number of input arguments')
end

% Add initialization period of 10 years to the start of the record
if length(P(:,1))<(10*365.25)
    repeats = ceil((10*365.25)/length(P(:,1)));
    init_period = int32(repeats*length(P(:,1)));
    P = [repmat(P,repeats,1); P];
    ETpot = [repmat(ETpot,repeats,1); ETpot];
    T = [repmat(T,repeats,1); T];
else
    init_period = int32(10*365.25);
    P = [P(1:init_period,:); P];
    ETpot = [ETpot(1:init_period,:); ETpot];
    T = [T(1:init_period,:); T];
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

SNOWPACK = nan(size(P)); SNOWPACK(1,:) = 0.0001;
MELTWATER = nan(size(P)); MELTWATER(1,:) = 0.0001;
SM = nan(size(P)); SM(1,:) = 0.0001;
SUZ = nan(size(P)); SUZ(1,:) = 0.0001;
SLZ = nan(size(P)); SLZ(1,:) = 0.0001;
ETact = nan(size(P)); ETact(1,:) = 0.0001;
Qsim = nan(size(P)); Qsim(1,:) = 0.0001;

for ii = length(P(1,:)):-1:1
    P(:,ii) = parPCORR(ii).*P(:,ii);
    SNOW(:,ii) = P(:,ii);
    RAIN(:,ii) = P(:,ii);
    RAIN(T(:,ii)<parTT(ii),ii) = 0;
    SNOW(T(:,ii)>=parTT(ii),ii) = 0;
    SNOW(:,ii) = SNOW(:,ii).*parSFCF(ii);
end

for t = 2:length(P)
    
    % Snow
    SNOWPACK(t,:) = SNOWPACK(t-1,:)+SNOW(t,:);
    melt = parCFMAX .* (T(t,:)-parTT);
    melt(melt<0) = 0;
    melt = min(melt,SNOWPACK(t,:));
    MELTWATER(t,:) = MELTWATER(t-1,:)+melt;
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
    soil_wetness = (SM(t-1,:)./parFC).^parBETA;
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
    SUZ(t,:) = SUZ(t-1,:)+recharge+excess;
    PERC = min(SUZ(t,:), parPERC);
    SUZ(t,:) = SUZ(t,:)-PERC;
    Q0 = parK0 .* max(SUZ(t,:)-parUZL, 0.0);
    SUZ(t,:) = SUZ(t,:)-Q0;
    Q1 = parK1.*SUZ(t,:);
    SUZ(t,:) = SUZ(t,:)-Q1;
    SLZ(t,:) = SLZ(t-1,:)+PERC;
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
        Qsim_smoothed(bb:length(Qsim_smoothed),ii) = Qsim_smoothed(bb:length(Qsim_smoothed),ii)...
            +Qsim(1:length(Qsim_smoothed)-bb+1,ii)*w_small(bb);
    end
end
Qsim = Qsim_smoothed;

% Remove initialization period
Qsim(1:init_period,:) = [];
SM(1:init_period,:) = [];
SUZ(1:init_period,:) = [];
SLZ(1:init_period,:) = [];
SNOWPACK(1:init_period,:) = [];
ETact(1:init_period,:) = [];