function para_ens = func_perturb_parameters(para_nom,input,configure)
    %{
    func_perturb_parameters generates ensemble members of a model parameter
    Syntax:  
		para_ens = func_perturb_parameters(para_nom,input,configure)
    Inputs:
		para_nom - A model parameter (nominal value)
        input    - Essential input data
			     - Required parameters
				     - forcing_coord: Coordinates (degree) of the study domain (longitude, latitude)
		configure - Configurations/settings described the way a parameter is perturbed
                  - Required parameters
                      - nens: Number of ensemble members to generate (e.g., 100)
                      - para_err: Magnitude of the perturbation in percent (e.g., 10)
                      - pparameters_clen: Correlation length in degree (e.g., 1.0)
                                          pparameters_clen = 0 when no correlation error is required.
    Output:
        para_ens - Ensemble parameters (2-dimension; ngrid,nens)
    Examples: 
        para_ens = func_perturb_parameters(para_nom,input,configure);
    
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
    %}

    %------------- BEGIN CODE --------------
    ngrid = size(para_nom,2);
    if configure.pparameters_clen == 0
        err = randn(ngrid,configure.nens);
    else
        Sigma = func_get_correlation_matrix(input.forcing_coord,configure.pparameters_clen); % including correlation error
        err   = mvnrnd(zeros(ngrid,1),Sigma,configure.nens)'; % (ngrid,nens)
    end
    err = err - repmat(mean(err,2,'omitnan'),1,size(err,2)); % Recenter Gaussian distribution
    para_ens = NaN(ngrid,configure.nens);
    for ens = 1:configure.nens
        para_ens(:,ens) = para_nom + (configure.para_err/100).*para_nom.*err(:,ens)';
        for gi = 1:ngrid
            if para_ens(gi,ens) < 0
                para_ens(gi,ens) = para_nom(gi);
            end
        end
    end
    %------------- END OF CODE --------------
end