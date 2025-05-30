% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';
nu = 1e-5;
disease_timing = 34.6806;
end_time = disease_timing*365; % Only dose up to time of disease
scale = 1e7;

% Define r1 and r2 Values for  early and late disease and high and low nu
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

if strcmp(disease_type,'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
elseif strcmp(disease_type,'late')
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

r1 = 1;
r2 = 1;

Init_B = 1*10^6;

params = LoadParameters();

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

% Figure Name
figurename = ['DosingNoDepletionFigures/PhasePlanePlots', '.pdf'];

% Save Parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Preallocate Space for T_vals
RL_vals = [0, 10, 1e3];

T_vals = zeros(length(RL_vals), 1);
tspan1 = 0:0.01:70*365;

AL_initcond = 0;
EL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;


kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));

figure(1)

for i = 1:length(RL_vals)

        % Vary RL(0)
        RL_initcond = RL_vals(i);


        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
       
              
        % Run the Model
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond);
        
        % Relabel Compartments to Easily Keep Track
        AL = Y(:, 1);
        EL = Y(:, 2);
        RL = Y(:, 3);
        AP = Y(:, 4);
        EP = Y(:, 5);
        RP = Y(:, 6);
        B = Y(:, 7);
        
        T_vals(i) = round(T(end)./365);

        figure(1)
        plot(EL, RL, 'LineWidth', 1.8, 'DisplayName', ['$R_L(0) = ', num2str(RL_initcond), '$'])
        hold on

        figure(2)
        plot(T./365, B./Init_B, 'LineWidth', 1.8, 'DisplayName', ['$R_L(0) = ', num2str(RL_initcond), '$'])
        hold on

end

figure(2)
xlabel('Time, years')
ylabel('\beta/\beta(0)')
ylim([0.9 1])
xlim([0 70])
yticks([0.9 1])
yline(0.2*Init_B/Init_B, '--', 'LineWidth', 1.3, 'Color', 0.5*[1 1 1], 'DisplayName',...
    '$20\% \hspace{0.5em} \beta \hspace{0.5em} cells$')
ax = gca;
ax.FontSize = 30;% Set Plot Size
set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')
h_leg = legend('Interpreter', 'latex', 'Location', 'Northeast');
HeightScaleFactor = 2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(4) = NewHeight;

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
h_leg = legend('Interpreter', 'latex', 'Location', 'Northeast');
HeightScaleFactor = 2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(4) = NewHeight;
hold on
ax = gca;
ax.FontSize = 30;% Set Plot Size
xlim([0 kE])
ylim([0, inf]);

set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')


set(gcf, 'Position', [100, 300, 900, 600]);
set(gcf, 'Color', 'White')
% set(gca, 'Xscale', 'log')

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

