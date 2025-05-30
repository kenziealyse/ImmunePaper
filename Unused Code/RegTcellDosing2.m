% CLEAR THE WORKSPACE
clc
clear
close all

% Dosing Information
Dose_vals = 0;%0:1e8:1e10;
Dose_Time_vals = 0.1;%(1:.1:6).*365;

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

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 0.5*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);
nu = 0.01;%1*10^(-5);%

figurename =['Figures\dosetimevsamount', num2str(nu),'.pdf'];

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

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

        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

        tspan1 = 0:0.01:Dose_Time;
              
        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);
        
        % T_predose(end)
        % Dose_Time(end)
        if T_predose(end)/365 < Dose_Time/365

            T_vals(j,i) = round(T_predose(end)./365); 
           
        else 
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

        end

    end

end

% Define custom colormap
figure(1)
imagesc(Dose_Time_vals./365, Dose_vals, T_vals);
xlim([min(Dose_Time_vals)./365 max(Dose_Time_vals)./365])
ylim([min(Dose_vals) max(Dose_vals)])
xlabel('Dose Times', 'FontSize', 17);
ylabel('Dose Amounts', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Difference in Time to 20% Beta Cell Mass');

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
