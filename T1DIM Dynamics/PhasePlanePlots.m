%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: PhasePlanePlots.m
%
% Description:
% This script simulates and plots phase plane trajectories (EL vs RL) and
% beta cell mass dynamics over time for different initial regulatory T-cell 
% populations (RL_initcond). It uses a predefined ODE model to analyze disease
% progression in terms of beta cell loss, highlighting differences by initial conditions.
%
% The script saves the resulting plots as a PDF file in the Figures folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';        % Choose 'early' or 'late' disease parameter set
nu = 1e-5;                    % Parameter nu value
disease_timing = 34.6806;     % Approximate disease onset in years
end_time = disease_timing*365; % Convert disease onset to days for simulation end
scale = 1e7;                  % Scale factor for plotting or parameters

% Define r1 and r2 Values for early and late disease
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

% Assign r1 and r2 based on disease type selected
if strcmp(disease_type,'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
elseif strcmp(disease_type,'late')
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

% Override r1 and r2 for this run (can comment out if undesired)
r1 = 1;
r2 = 1;

% Initial Beta cell mass
Init_B = 1*10^6;

% Load fixed model parameters
params = LoadParameters();

% Extract individual parameters from vector for clarity
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

% Define the output figure name
figurename = ['../Figures/PhasePlanePlots', '.pdf'];

% Combine parameters including r1, r2, and nu into a single vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Initial values of RL to simulate
RL_vals = [0, 10, 1e3];

% Preallocate vector for final time points (years)
T_vals = zeros(length(RL_vals), 1);

% Time span for ODE simulation: from day 0 to 70 years (in days)
tspan1 = 0:0.01:70*365;

% Initial conditions for compartments: AL, EL, RL, AP, EP, RP, B
AL_initcond = 0;
EL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

% Calculate constants kR and kE for plotting phase plane boundaries
kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));

% Create a new figure for plotting EL vs RL trajectories
figure(1)

% Loop over different initial RL values
for i = 1:length(RL_vals)

        % Set initial RL value for this run
        RL_initcond = RL_vals(i);

        % Initialize state vector with current RL_initcond
        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
       
        % Run ODE simulation using NuRegTcellmodel and given parameters
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond);
        
        % Extract solution compartments for plotting and analysis
        AL = Y(:, 1);
        EL = Y(:, 2);
        RL = Y(:, 3);
        AP = Y(:, 4);
        EP = Y(:, 5);
        RP = Y(:, 6);
        B = Y(:, 7);
        
        % Save final time point in years for this initial condition
        T_vals(i) = round(T(end)./365);

        % Plot EL vs RL phase plane trajectory for this initial RL
        figure(1)
        plot(EL, RL, 'LineWidth', 1.8, 'DisplayName', ['$R_L(0) = ', num2str(RL_initcond), '$'])
        hold on

        % Plot beta cell mass over time in separate figure
        figure(2)
        plot(T./365, B./Init_B, 'LineWidth', 1.8, 'DisplayName', ['$R_L(0) = ', num2str(RL_initcond), '$'])
        hold on

end

% Customize beta cell mass plot (Figure 2)
figure(2)
xlabel('Time, years')
ylabel('\beta/\beta(0)')
ylim([0.9 1])
xlim([0 70])
yticks([0.9 1])
% Add horizontal line marking 20% beta cell mass threshold
yline(0.2*Init_B/Init_B, '--', 'LineWidth', 1.3, 'Color', 0.5*[1 1 1], 'DisplayName',...
    '$20\% \hspace{0.5em} \beta \hspace{0.5em} cells$')

% Set font size for axes
ax = gca;
ax.FontSize = 30;

% Set figure size and background color
set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')

% Configure legend with LaTeX interpreter and increase height
h_leg = legend('Interpreter', 'latex', 'Location', 'Northeast');
HeightScaleFactor = 2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(4) = NewHeight;

% Customize phase plane plot (Figure 1) with boundary lines and markers
figure(1)
plot([kE 0], [0 kR], '--', 'Color', [0.5, 0.5, 0.5],...
    'HandleVisibility', 'on', 'linewidth', 2, 'DisplayName', ...
    '$R_L = -E_L + C(1 - \frac{1}{k_R})$')
plot(kE, 0, 'rp', 'MarkerSize', 35, 'MarkerFaceColor', "#77AC30", ...
    'MarkerEdgeColor', "#77AC30", 'DisplayName', '$C(1 - \frac{1}{k_E})$')
plot(0, kR,  'rp', 'MarkerSize', 35, 'MarkerFaceColor', 	"#7E2F8E", ...
    'MarkerEdgeColor', "#7E2F8E", 'DisplayName', '$C(1 - \frac{1}{k_R})$')

ylabel('R_L')
xlabel('E_L')
legend show

% Configure legend with LaTeX interpreter and increase height
h_leg = legend('Interpreter', 'latex', 'Location', 'Northeast');
HeightScaleFactor = 2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(4) = NewHeight;

hold on

% Set font size for axes
ax = gca;
ax.FontSize = 30;

% Set axis limits for phase plane plot
xlim([0 kE])
ylim([0, inf]);

% Set figure size and background color
set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')

% Adjust figure size for final saving
set(gcf, 'Position', [100, 300, 900, 600]);
set(gcf, 'Color', 'White')

% Save the figure as pdf in the Figures folder
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder