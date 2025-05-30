function [] = RunModel(vary_parameters, ...
    lower_limit_param1, upper_limit_param1, N_param1, ...
    lower_limit_param2, upper_limit_param2, N_param2, ...
    disease_timing, ...
    nu_type, ...
    Dosing, ...
    dose_time, dose_amount)



if strcmp(vary_parameters{1}, 'NA')
    
    % Load Parameter Values
    params = LoadParameters();

    % Specify Parameter Values
    phiA = params(1);
    deltaA = params(2);
    lambdaEL = params(3);
    omegaEL = params(4);
    phiE = params(5);
    deltaE = params(6);
    lambdaR = params(7);
    omegaR = params(8);
    C = params(9);
    phiR = params(10);
    deltaR = params(11);
    kappa = params(15);

    


end