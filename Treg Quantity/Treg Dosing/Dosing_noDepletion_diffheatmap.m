%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Name: Dosing_noDepletion_diffheatmap.m
%
% Description:
% This script simulates disease progression in a T1D model under varying 
% immune regulatory T-cell (RL) dosing amounts and times. For each combination 
% of dose amount and timing, it computes the resulting delay (or acceleration) 
% in time to 20% beta cell mass loss compared to the untreated case.
%
% Output:
% - A heatmap showing the difference in time to critical beta cell mass loss
%   as a function of dose time and dose amount.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% -----------------------------
% User-defined parameters
% -----------------------------
disease_type = 'early';   % Set disease type: 'early' or 'late'
nu = 0.01;                % Regulatory strength parameter
disease_timing = 5.9966;  % Time to disease in untreated case (years)
end_time = disease_timing * 365; % Simulate dosing up to disease onset

% Dosing ranges
N = 300;
Dose_vals = 10.^(5:1:10);                          % Range of dose sizes
Dose_Time_vals = linspace(1, end_time, N);        % Range of dose times

% -----------------------------
% Set r1 and r2 based on disease type
% -----------------------------
earlydiseasepoint = [0.2e-5, 0.2e-5];
latediseasepoint  = [0.8e-5, 1e-5];

if strcmp(disease_type,'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
else
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

Init_B = 1e6;  % Initial beta cell mass

% Load and unpack fixed model parameters
params = LoadParameters();
phiA = params(1);  deltaA = params(2); lambdaEL = params(3);
omegaEL = params(4); phiE = params(5); deltaE = params(6);
lambdaR = params(7); omegaR = params(8); C = params(9);
phiR = params(10); deltaR = params(11); kappa = params(15);

% Store all parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

% Output file name
figurename = ['../../Figures/DosingdiffHeatMapNoDepletion_', ...
              disease_type, '_nuis', num2str(nu), '.pdf'];

% Initialize results matrix
T_vals = zeros(length(Dose_vals), length(Dose_Time_vals));

% -----------------------------
% Simulate across dose times and amounts
% -----------------------------
for i = 1:length(Dose_Time_vals)
    Dose_Time = Dose_Time_vals(i);
    for j = 1:length(Dose_vals)
        Dose_amount = Dose_vals(j);

        % --- Pre-dose simulation ---
        init_cond = [0; 10; 10; 0; 0; 0; Init_B];  % Initial conditions
        tspan1 = 0:0.01:Dose_Time;
        options = odeset('Events', @(t, Y) ...
            PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_predose, Y_predose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), ...
                                        tspan1, init_cond, options);

        % Apply dose
        AL = Y_predose(end, 1); EL = Y_predose(end, 2); RL = Y_predose(end, 3);
        AP = Y_predose(end, 4); EP = Y_predose(end, 5); RP = Y_predose(end, 6);
        B  = Y_predose(end, 7);
        RL = RL + Dose_amount;

        % --- Post-dose simulation ---
        init_cond = [AL; EL; RL; AP; EP; RP; B];
        tspan2 = Dose_Time:0.01:70*365;
        options = odeset('Events', @(t, Y) ...
            PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_postdose, Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), ...
                                          tspan2, init_cond, options);

        % Combine and record time to disease
        T = [T_predose; T_postdose];
        T_vals(j,i) = round(T(end)./365);  % Convert to years
    end
end

% -----------------------------
% Plot heatmap of delays or accelerations in disease onset
% -----------------------------
figure(1)
imagesc(Dose_Time_vals./365, Dose_vals, T_vals - round(disease_timing));
set(gca,'YDir','normal');  % Standard orientation
xlim([min(Dose_Time_vals)./365 max(Dose_Time_vals)./365])
ylim([1e5 max(Dose_vals)])
xlabel('Dose Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
title('Difference in Time to 20% Beta Cell Mass', 'FontSize', 17);
ax = gca;
ax.FontSize = 21;
yticks([1e5 1e6 1e7 1e8 1e9 1e10])
set(gca, 'Yscale', 'log')

% Custom colormap from red → blue → green
numColors = 100;
redValues   = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues  = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];
customColorMap = [redValues', greenValues', blueValues']; 
colormap(gca, customColorMap);
colorbar;

% Set z-axis color limits based on max deviation
zlimit = ceil(max(max(T_vals - disease_timing))/10)*10;
if zlimit == 0
    zlimit = 10;
end
clim([0, zlimit]);

% Format and save figure
set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure
