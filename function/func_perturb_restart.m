function restart_ens = func_perturb_restart(nom,input,configure)
    %{
    func_perturb_restart generates ensemble members of a model initial state
    Syntax:  
		restart_ens = func_perturb_restart(nom,input,configure)
    Inputs:
		nom   - A model initial state (nominal value)
        input - Essential input data
			  - Required parameters
				  - forcing_coord: Coordinates (degree) of the study domain (longitude, latitude)
		configure - Configurations/settings described the way a parameter is perturbed
                  - Required parameters
                      - nens: Number of ensemble members to generate (e.g., 100)
                      - restart_err: Magnitude of the perturbation in percent (e.g., 10)
                      - restart_clen: Correlation length in degree (e.g., 1.0)
                                      restart_clen = 0 when no correlation error is required.
    Output:
        restart_ens - Ensemble initial states (2-dimension; ngrid,nens)
    Examples: 
        restart_ens = func_perturb_restart(nom,input,configure);
    
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
    %}

    %------------- BEGIN CODE --------------
    ngrid = size(nom,2);
    if configure.restart_clen == 0
        err = randn(ngrid,configure.nens);
    else
        Sigma = func_get_correlation_matrix(input.forcing_coord,configure.restart_clen); % including correlation error
        err   = mvnrnd(zeros(ngrid,1),Sigma,configure.nens)'; % (ngrid,nens)
    end
    err = err - repmat(mean(err,2,'omitnan'),1,size(err,2)); % Recenter Gaussian distribution
    restart_ens = NaN(ngrid,configure.nens);
    for ens = 1:configure.nens
        restart_ens(:,ens) = nom + (configure.restart_err/100).*nom.*err(:,ens)';
        for gi = 1:ngrid
            if restart_ens(gi,ens) < 0
                restart_ens(gi,ens) = nom(gi);
            end
        end
    end
    restart_ens(restart_ens<0) = 0;
    %------------- END OF CODE --------------
end