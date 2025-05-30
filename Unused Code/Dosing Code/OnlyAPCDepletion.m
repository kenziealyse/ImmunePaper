% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';
nu = 0.01;
disease_time = 64.0012;
end_time = disease_time*365; % Only dose up to time of disease

% Dosing Information
N = 50;

Dose_Time_vals = linspace(1, end_time, N);

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
figurename = ['APCDepletionOnlyFigures/TTDvsDepletion_', disease_type, '_nuis', num2str(nu), '.pdf'];

% Save Parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Preallocate Space for T_vals
T_vals = zeros(length(Dose_Time_vals), 1);
tspan = 0:0.01:70*365;

for i = 1:length(Dose_Time_vals)

       Dose_Time = Dose_Time_vals(i);

        %% Pre Dose      
        AL_initcond = 0;
        EL_initcond = 10;
        RL_initcond = 10;
        AP_initcond = 0;
        EP_initcond = 0;
        RP_initcond = 0;
        B_initcond = 1*10^6;

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
               
        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), Init_B));
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), tspan, init_cond, options);

        T_vals(i) = T(end)./365;

end

% Define custom colormap
figure(1)
plot(Dose_Time_vals./365, T_vals, 'Color', 'k', 'LineWidth', 1.3, 'HandleVisibility','on');
hold on
yline(disease_time, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.3)
xlabel('APC Depletion Times, years', 'FontSize', 17);
ylabel('Time to 20% Beta Cell Mass', 'FontSize', 17);
ax = gca;
ax.FontSize = 21;
ylim([0 ceil(max(T_vals/10))*10])
legend('Time to disease \bf{with} APC depletion', 'Time to disease \bf{without} APC depletion', ...
    'FontSize', 17, 'Location', 'best')

% Set Plot Size
set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

