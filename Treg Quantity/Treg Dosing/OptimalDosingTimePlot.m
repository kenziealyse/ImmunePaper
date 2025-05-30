%% OptimalDosingTimePlot.m
% This script simulates and evaluates the timing of immunotherapy dosing 
% in a disease model by running ODE simulations with varying dose times.
% It calculates the product of effector peptide (E_P) and beta cell mass (β) at dosing 
% and correlates it with the time to significant beta cell loss (20% mass).
% The results are visualized in a scatter plot and saved as a PDF.

% CLEAR THE WORKSPACE
clc
clear
close all

% User-defined parameters
disease_type = 'late';        % Disease stage: 'early' or 'late'
nu = 1e-5;                   % Model parameter nu
disease_timing = 30;         % Time of disease onset in years
end_time = disease_timing*365; % Simulation end time in days (up to disease onset)

% Dosing Information
Dose_amount = 10^(9);         % Amount of dose added to RL compartment
N = 50;                      % Number of dosing time points to simulate
Dose_Time_vals = linspace(1, end_time, N); % Array of dosing times (days)

% Define r1 and r2 values for early and late disease based on user input
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

if strcmp(disease_type,'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
elseif strcmp(disease_type,'late')
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

% Initial beta cell mass
Init_B = 1*10^6;

% Load fixed parameters from external function
params = LoadParameters();

% Unpack parameters for clarity
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

% Define output figure name
figurename = ['../../Figures/OptimalDosingSimsplot1_', disease_type, '_nuis', num2str(nu), '.pdf'];

% Combine parameters into vector including r1, r2, and nu
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

% Preallocate arrays to store results
T_vals = zeros(length(Dose_Time_vals), 1);  % Time to 20% beta cell mass remaining
EtimesB = zeros(length(Dose_Time_vals), 1); % E_P * beta cell mass at dosing time

% Loop over all dose times
for i = 1:length(Dose_Time_vals)

    Dose_Time = Dose_Time_vals(i);

    %% Pre-dose simulation: from time 0 to Dose_Time
    % Initial conditions for all compartments before dosing
    AL_initcond = 0;
    EL_initcond = 10;
    RL_initcond = 10;
    AP_initcond = 0;
    EP_initcond = 0;
    RP_initcond = 0;
    B_initcond = 1*10^6;

    tspan1 = 0:0.01:Dose_Time;

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                 EP_initcond RP_initcond B_initcond]';

    % Run ODE solver with event detection for 20% beta cell mass
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
    [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

    % Add dose to RL compartment at dosing time
    AL_initcond = Y_predose(end, 1);
    EL_initcond = Y_predose(end, 2);
    RL_initcond = Y_predose(end, 3) + Dose_amount; % Add dose here
    AP_initcond = Y_predose(end, 4);
    EP_initcond = Y_predose(end, 5);
    RP_initcond = Y_predose(end, 6);
    B_initcond = Y_predose(end, 7);

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                 EP_initcond RP_initcond B_initcond]';

    % Post-dose simulation: from Dose_Time to 70 years
    tspan2 = Dose_Time:0.01:70*365;

    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
    [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan2, init_cond, options);

    % Combine pre-dose and post-dose time and solution arrays
    T = [T_predose; T_postdose];
    Y = [Y_predose; Y_postdose];

    % Extract compartments for convenience
    AL = Y(:, 1);
    EL = Y(:, 2);
    RL = Y(:, 3);
    AP = Y(:, 4);
    EP = Y(:, 5);
    RP = Y(:, 6);
    B = Y(:, 7);

    % Store time (in years) to reach 20% beta cell mass
    T_vals(i) = round(T(end)./365);

    % Find index in time vector corresponding to dosing time
    dose_index = find(T == Dose_Time);

    % Calculate E_P * beta cell mass at dose time
    EtimesB(i) = EP(dose_index(1)) * B(dose_index(1));
end

% Plot scatter of E_P * beta cell mass at dosing vs time to beta cell loss
figure(2)
scatter(EtimesB, T_vals, 'filled', 'SizeData', 75, "MarkerEdgeColor","b", ...
    "MarkerFaceColor", [0 0.7 0.7]) 
ylabel('Time to 20% Beta Cell Mass')
xlabel('E_P \times \beta at dose time')
ax = gca;
ax.FontSize = 25; % Set axis font size
set(gcf, 'Position', [100, 300, 900, 400]);
set(gcf, 'Color', 'White')

% Save figure as PDF to Figures folder
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder