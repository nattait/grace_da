function forcing_ens = func_perturb_forcing(input,configure,forcing_type)
    %{
    func_perturb_forcing generates ensemble members of a meteorological field
    Syntax:  
		forcing_ens = func_perturb_forcing(input,configure,forcing_type)
    Inputs:
		input - Meteorological forcing fields
			  - Required parameters
				  - P: Precipitation
				  - T: Temperature
				  - ETpot: Potential evapotranspiration
				  - forcing_coord: Coordinates (degree) of forcing data (longitude, latitude)
		configure - Configurations/settings described the way forcing data are perturbed
                  - Required parameters
                      - nens: Number of ensemble members to generate (e.g., 100)
                      - forc_err: Error size P (%), T (^oC), and ETpot (%) (e.g., [10 2 10])
                      - pforcing_clen: Correlation length of P, T, and ETpot in degree (e.g., [1.0 1.0 1.0])
                                       pforcing_clen = 0 when no correlation error is required.
                      - forc_err_range: magnitude of the perturbation
                        options: 'median','minmax','nominal'
		forcing_type  - Type of input forcing data
                        options: 'P','T','ETpot'
    Output:
        forcing_ens - Ensemble forcing fields (3-dimension; ntime,ngrid,nens)
    Examples: 
        P_ens     = func_perturb_forcing(input,configure,'P');
        T_ens     = func_perturb_forcing(input,configure,'T');
        ETpot_ens = func_perturb_forcing(input,configure,'ETpot');   
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
    %}

    %------------- BEGIN CODE --------------
        if strcmp(forcing_type,'P'), index=1; forcing = input.P;
        elseif strcmp(forcing_type,'T'), index=2; forcing = input.T;
        elseif strcmp(forcing_type,'ETpot'), index=3; forcing = input.ETpot; 
        end
    
        if strcmp(configure.forc_err_range,'median')
            forcing_range = repmat(median(forcing),size(forcing,1),1); % use median value of each day
        elseif strcmp(configure.forc_err_range,'minmax')
            forcing_range = repmat((max(forcing)-min(forcing)),size(forcing,1),1); % use median value of each day
        elseif strcmp(configure.forc_err_range,'nominal')
            forcing_range = forcing;
        end

        ntime = size(forcing,1);
        ngrid = size(forcing,2);
        if configure.pforcing_clen(index) == 0 % No spatial correlation error
            temp_err = randn(ntime,ngrid,configure.nens);
        else
            Sigma    = func_get_correlation_matrix(input.forcing_coord,configure.pforcing_clen(index)); % including correlation error
            temp_err = zeros(ntime,ngrid,configure.nens);
            for ti=1:ntime % perturb every day
                temp_err(ti,:,:) = mvnrnd(zeros(ngrid,1),Sigma,configure.nens)'; % (ngrid,nens)
            end
        end
        temp_err = temp_err - repmat(mean(temp_err,3),1,1,size(temp_err,3)); % Recenter Gaussian distribution
        forcing_ens = zeros(ntime,ngrid,configure.nens);
        for ens = 1:configure.nens
            if strcmp(forcing_type,'T') % Applying additive noises
                forcing_ens(:,:,ens) = forcing + configure.forc_err(index).*temp_err(:,:,ens);
            elseif strcmp(forcing_type,'P') || strcmp(forcing_type,'ETpot') % Applying multiplicative noises
                mlp = log(1 + ((configure.forc_err(index)/100).*forcing_range).^2);
                forcing_ens(:,:,ens) = forcing .* exp(-0.5*mlp + sqrt(mlp).*temp_err(:,:,ens));
            end
        end
    %------------- END OF CODE --------------
end

function Mcor = func_get_correlation_matrix(coord,cor_length)
	%{
    func_get_correlation_matrix generates a correlation matrix from the given coordinates and correlation length
    Syntax:  
		Mcor = func_get_correlation_matrix(coord,cor_length)
    Inputs:
		coord      - Coordinates of the forcing data (2-dimentional array)
		cor_length - Correlation length of the forcing data (a scalar value)
    Output:
		Mcor - Correlation matrix
    Example: 
		Mcor = func_get_correlation_matrix(coord,0.1);
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
	%}
	
	%------------- BEGIN CODE --------------
    Mcor = zeros(size(coord,1),size(coord,1));
    for p1 = 1:size(coord,1)
        for p2 = 1:p1
            [dis,~] = distance(coord(p1,2),coord(p1,1), ...
                                    coord(p2,2),coord(p2,1));
            Mcor(p1,p2) = exp(-dis/cor_length);
            Mcor(p2,p1) = Mcor(p1,p2);
        end
    end
	%------------- END OF CODE --------------
end