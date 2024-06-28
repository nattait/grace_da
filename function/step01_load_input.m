function input = step01_load_input(filename)
    %{
    load_input reads and obtain dimension from the input file
    Syntax:  
		input = step01_load_input(filename)
    Inputs:
		filename - Name of the input file (e.g., data\input.mat)
    Output:
        input - A variable contained all essential data
    Examples: 
       input = step01_load_input('data\input.mat');   
    
    -----------------------------
    Description:
         Refer to description.m
    -----------------------------
    %}

    %------------- BEGIN CODE --------------
    load(filename,'input');

    %%%%% Identifying simulation period and number of grid cells
    input.ntime = size(input.P,1);
    input.ngrid = size(input.P,2);
    input.nlon  = length(input.lon);
    input.nlat  = length(input.lat);
    %------------- END OF CODE --------------
end