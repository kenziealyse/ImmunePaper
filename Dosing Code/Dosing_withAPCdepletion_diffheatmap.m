% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';
nu = 0.01;
disease_timing = 64.0012;
end_time = disease_timing*365; % Only dose up to time of disease

% Dosing Information
N = 300;
Dose_vals = [10.^(0:1:10)];
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
figurename = ['Combo/DosingdiffHeatMapwithDepletion_', disease_type, '_nuis', num2str(nu), '.pdf'];

% Save Parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Preallocate Space for T_vals
T_vals = zeros(length(Dose_vals), length(Dose_Time_vals));

for i = 1:length(Dose_Time_vals)

    Dose_Time = Dose_Time_vals(i);

    for j = 1:length(Dose_vals)

        Dose_amount = Dose_vals(j);

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
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), tspan1, init_cond, options);

        AL_initcond = Y_predose(end, 1);
        EL_initcond = Y_predose(end, 2);
        RL_initcond = Y_predose(end, 3)+ Dose_amount;
        AP_initcond = Y_predose(end, 4);
        EP_initcond = Y_predose(end, 5);
        RP_initcond = Y_predose(end, 6);
        B_initcond = Y_predose(end, 7);

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

        tspan2 = Dose_Time:0.01:70*365;

        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), Init_B));
        [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time), tspan2, init_cond, options);
                
        T = [T_predose; T_postdose];
        Y = [Y_predose; Y_postdose];

        T_vals(j,i) = round(T(end)./365);


    end

end

% Define custom colormap
figure(1)
imagesc(Dose_Time_vals./365, Dose_vals, T_vals - disease_timing);
xlim([min(Dose_Time_vals)./365 max(Dose_Time_vals)./365])
ylim([min(Dose_vals) max(Dose_vals)])
xlabel('Treg Dose and APC Depletion Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Difference in Time to 20% Beta Cell Mass', 'FontSize', 17);
ax = gca;
ax.FontSize = 21;
yticks([1e5 1e6 1e7 1e8 1e9 1e10])
set(gca, 'Yscale', 'log')

% Define the number of colors
numColors = 100;

% Define transitions for red, green, and blue
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)]; % Red fades to zero
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)]; % Green grows
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)]; % Blue in middle

% Combine into a custom colormap
customColorMap = [redValues', greenValues', blueValues']; 

% Apply the colormap
colormap(gca, customColorMap);

% Show color bar with custom tick labels
colorbar;

% Set z axis limits
zlimit = ceil(max(max(T_vals - disease_timing))/10)*10;
if zlimit == 0
    zlimit = 10;
end
clim([0, zlimit]);  % Set color axis limits based on data range

% Set Plot Size
set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

