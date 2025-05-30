%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: DDosing_withAPCdepletion_diffheatmap.m
% Description:
%   This script simulates the effects of varying both the timing and dose
%   of regulatory T cell therapy and APC depletion. The results are 
%   visualized as a heatmap showing the difference in time
%   to 20% beta cell mass compared to baseline disease progression.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% -------------------------- User Parameters -----------------------------

disease_type = 'late';         % 'early' or 'late' disease
nu = 0.01;                     % Effector activation rate
disease_timing = 64.0012;      % Baseline disease onset in years
end_time = disease_timing*365; % Simulation ends at disease onset

% Dosing setup
N = 300;  % Number of discrete dose timings
Dose_vals = 10.^(0:1:10);    % Dose amounts (log scale)
Dose_Time_vals = linspace(1, end_time, N);  % Dosing time vector

% -------------------------- Disease Parameters --------------------------

% r1 and r2 determine immune regulation thresholds
earlydiseasepoint = [0.2e-5, 0.2e-5];
latediseasepoint  = [0.8e-5, 1.0e-5];

if strcmp(disease_type,'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
elseif strcmp(disease_type,'late')
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

% Initial beta cell count
Init_B = 1e6;

% Load base model parameters
params = LoadParameters();
phiA = params(1); deltaA = params(2); lambdaEL = params(3); omegaEL = params(4);
phiE = params(5); deltaE = params(6); lambdaR = params(7); omegaR = params(8);
C = params(9); phiR = params(10); deltaR = params(11); kappa = params(15);

% Save parameters to vector with updated r1, r2, nu
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

% Set output filename
figurename = ['../../Figures/DosingdiffHeatMapwithDepletion_', ...
              disease_type, '_nuis', num2str(nu), '.pdf'];

% ----------------------- Preallocate Result Matrix ----------------------

T_vals = zeros(length(Dose_vals), length(Dose_Time_vals));

% ----------------------- Simulation Loop -------------------------------

for i = 1:length(Dose_Time_vals)
    Dose_Time = Dose_Time_vals(i);

    for j = 1:length(Dose_vals)
        Dose_amount = Dose_vals(j);

        % Initialize compartments
        AL = 0; EL = 10; RL = 10;
        AP = 0; EP = 0; RP = 0;
        B  = Init_B;

        % Simulate pre-dose phase
        tspan1 = 0:0.01:Dose_Time;
        init_cond = [AL EL RL AP EP RP B]';
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, ...
                  NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) ...
                  NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), ...
                  tspan1, init_cond, options);

        % Add dose to RL at dosing time
        Y_end = Y_predose(end, :);
        Y_end(3) = Y_end(3) + Dose_amount;  % RL += dose
        init_cond = Y_end';

        % Simulate post-dose phase
        tspan2 = Dose_Time:0.01:70*365;
        [T_postdose,Y_postdose] = ode23s(@(t,Y) ...
                  NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), ...
                  tspan2, init_cond, options);

        % Combine time courses
        T = [T_predose; T_postdose];
        Y = [Y_predose; Y_postdose];

        % Record outcome time (to 20% beta cell mass)
        T_vals(j,i) = round(T(end)/365);
    end
end

% ------------------------- Plotting Results -----------------------------

% Create heatmap of dosing outcomes
figure(1)
imagesc(Dose_Time_vals./365, Dose_vals, T_vals - disease_timing);
xlim([min(Dose_Time_vals)./365 max(Dose_Time_vals)./365])
ylim([min(Dose_vals) max(Dose_vals)])
xlabel('Treg Dose and APC Depletion Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
title('Difference in Time to 20% Beta Cell Mass', 'FontSize', 17);
set(gca,'YDir','normal');
set(gca, 'Yscale', 'log')
ax = gca;
ax.FontSize = 21;
yticks([1e5 1e6 1e7 1e8 1e9 1e10])

% Define custom colormap
numColors = 100;
redValues   = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues  = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];
customColorMap = [redValues', greenValues', blueValues'];
colormap(gca, customColorMap);
colorbar;

% Set color axis (clim) limits
zlimit = ceil(max(max(T_vals - disease_timing))/10)*10;
if zlimit == 0
    zlimit = 10;
end
clim([0, zlimit]);

% Set plot display properties
set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')

% Save figure
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', ...
         'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);  % Save figure to file