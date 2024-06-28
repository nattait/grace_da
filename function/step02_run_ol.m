function OL = step02_run_ol(input)
    %{
    run_ol performs the HBV model's openloop simulation
    Syntax:  
		OL = step02_run_ol(input)
    Inputs:
		input - A variable contained all essential input data
    Output:
		OL - A variable stored model outputs
    Example: 
		OL = step02_run_ol(input); % input variable is obtained from step01_load_input

    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
	%}
	
	%------------- BEGIN CODE --------------
    [OL.Qsim, OL.SM, OL.SUZ, OL.SLZ, OL.SNOWPACK, OL.ETact, OL.MELTWATER] = ...
            HBV_lumpedx(input.P,input.ETpot,input.T,input.parameters);
    OL.TWS = OL.SM+OL.SUZ+OL.SLZ;
    %------------- END OF CODE --------------
end