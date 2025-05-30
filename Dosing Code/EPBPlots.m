% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';
nu = 1e-5;
disease_timing = 30;
end_time = disease_timing*365; % Only dose up to time of disease

% Dosing Information
Dose_amount = 10.^(9);
N = 50;
Dose_Time_vals = linspace(1, end_time, N);%[1, 14, 182.5, 547.5]; %

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
figurename = ['DosingNoDepletionFigures/OptimalDosingSimsplot1_', disease_type, '_nuis', num2str(nu), '.pdf'];

% Save Parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Preallocate Space for T_vals
T_vals = zeros(length(Dose_Time_vals), 1);
EtimesB = zeros(length(Dose_Time_vals), 1);

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

        tspan1 = 0:0.01:Dose_Time;

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
       
              
        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

        AL_initcond = Y_predose(end, 1);
        EL_initcond = Y_predose(end, 2);
        RL_initcond = Y_predose(end, 3) + Dose_amount;
        AP_initcond = Y_predose(end, 4);
        EP_initcond = Y_predose(end, 5);
        RP_initcond = Y_predose(end, 6);
        B_initcond = Y_predose(end, 7);

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

        tspan2 = Dose_Time:0.01:70*365;

        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan2, init_cond, options);
                
        T = [T_predose; T_postdose];
        Y = [Y_predose; Y_postdose];

        % Relabel Compartments to Easily Keep Track
        AL = Y(:, 1);
        EL = Y(:, 2);
        RL = Y(:, 3);
        AP = Y(:, 4);
        EP = Y(:, 5);
        RP = Y(:, 6);
        B = Y(:, 7);

        T_vals(i) = round(T(end)./365);
        dose_index = find(T == Dose_Time);

        EtimesB(i) = EP(dose_index(1)).*B(dose_index(1));
end

figure(2)
scatter(EtimesB, T_vals, 'filled', 'SizeData', 75, "MarkerEdgeColor","b", ...
    "MarkerFaceColor", [0 0.7 0.7]) 
ylabel('Time to 20% Beta Cell Mass')
xlabel('E_P \times \beta at dose time')
ax = gca;
ax.FontSize = 25;% Set Plot Size
set(gcf, 'Position', [100, 300, 900, 00]);

% % RP and B
% figure(1)
% yyaxis left
% plot(Dose_Time_vals./365, EtimesB, '-o') 
% xlabel('Dose time, years')
% ylabel('E_P \times \beta at dose time')
% ax = gca;
% ax.FontSize = 25;% Set Plot Size
% 
% yyaxis right
% plot(Dose_Time_vals./365, T_vals, '-o') 
% ylabel('Time to 20% Beta Cell Mass')
% ax = gca;
% ax.FontSize = 25;% Set Plot Size


% set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')

% figure(2)
% % Plot before dosing 'too late'
% scatter(EtimesB, T_vals, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 2)
% % Plot after dosing 'too late'
% scatter(EtimesB, T_vals, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 2)
% ylabel('Time to 20% Beta Cell Mass')
% xlabel('E_P \times \beta')
% ax = gca;
% ax.FontSize = 21;% Set Plot Size

% set(gcf, 'Position', [100, 300, 1400, 800]);
set(gcf, 'Color', 'White')
% set(gca, 'Xscale', 'log')

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

