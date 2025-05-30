% CLEAR THE WORKSPACE
clc
clear
close all

% Dosing Information
N = 2;
Dose_vals = [0, 10, 1e2, 1e3, 1e4];
Dose_Time_vals = 0:3*365:70*365

% Define r1 and r2 Values for  early and late disease and high and low nu
r1_vals = [0.2*10^(-5), 0.8*10^(-5)];
r2_vals = [0.2*10^(-5), 1*10^(-5)];

nu_vals = [0.01, 1*10^(-5)];

% Set r1, r2, and nu
r1 = r1_vals(2);
r2 = r2_vals(2);
nu = nu_vals(2);

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
figurename = ['Figures/timetodisease_RLinitconds_samer1andr2_nuis', num2str(nu), '.pdf'];

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
        AP_initcond = 0;
        EP_initcond = 0;
        RP_initcond = 0;
        B_initcond = 1*10^6;

        if Dose_Time == 0
            tspan1 = 0:0.1:70*365;
            RL_initcond = 10 + Dose_amount;
        else
            tspan1 = 0:0.01:Dose_Time;
            RL_initcond = 10;
        end

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
       
              
        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

        if round(T_predose(end)/365,4) == round(Dose_Time/365,4)

            AL_initcond = Y_predose(end, 1);
            EL_initcond = Y_predose(end, 2);
            RL_initcond = Y_predose(end, 3);
            AP_initcond = Y_predose(end, 4);
            EP_initcond = Y_predose(end, 5);
            RP_initcond = Y_predose(end, 6) + Dose_amount;
            B_initcond = Y_predose(end, 7);
    
            init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';
    
            tspan2 = Dose_Time:0.01:70*365;
    
            % Run the Model
            options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
            [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan2, init_cond, options);
                    
            T = [T_predose; T_postdose];
            Y = [Y_predose; Y_postdose];
    
            T_vals(j,i) = round(T(end)./365);
        else

            T_vals(j,i) = round(T_predose(end)./365);
        end


    end

end

% Define custom colormap
figure(1)
imagesc(Dose_Time_vals./365, Dose_vals, T_vals);
xlim([min(Dose_Time_vals)./365 max(Dose_Time_vals)./365])
ylim([min(Dose_vals) max(Dose_vals)])
xlabel('Dose Times, days', 'FontSize', 17);
ylabel('Dose Amounts', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time to 20% Beta Cell Mass', 'FontSize', 17);

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
clim([0, 70]);  % Set color axis limits based on data range

% Set Plot Size
set(gcf, 'Position', [100, 300, 600, 400]);

    
% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder
