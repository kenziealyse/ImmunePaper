% Clear the workspace
clc
clear
close all

% Which parameters are you varying? (max is 2)

vary_parameters = {'r1'};

% Lower and upper limit and number of sample points of parameter 1 (if not needed, set to 0)

lower_limit_param1 = 0;
upper_limit_param1 = 1;
N_param1 = 20;

% Lower and upper limit and number of sample points of parameter 2 (if not needed, set to 0)

lower_limit_param2 = 0;
upper_limit_param2 = 1;
N_param2= 20;

% Are we looking at early or late disease onset? (if varying r1 or r2, use
% 'NA'

disease_timing = {'Late'};

% Are we looking at high or low nu? (if varying, use 'NA')

nu_type = {'Low'};

% Dosing?

Dosing = {'No'};

% If yes, what amount(s) and what time(s)? (if no dosing, set both to 0)

dose_time = 0;
dose_amount = 0;

% Run the model and save the results
RunModel(vary_parameters, ...
    lower_limit_param1, upper_limit_param1, N_param1, ...
    lower_limit_param2, upper_limit_param2, N_param2, ...
    disease_timing, ...
    nu_type, ...
    Dosing, ...
    dose_time, dose_amount)