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
Dose_amount = 10^(10);%[0, 10.^(0:1:10)];
Dose_Time = 10*365;%linspace(1, end_time, N);
frac_vals = 0:.01:1;

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

% Save Parameters in vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

% Preallocate Space for T_vals
T_vals = zeros(length(frac_vals), 1);

for i = 1:length(frac_vals)

    frac = frac_vals(i);

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
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time, frac), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time, frac), tspan1, init_cond, options);

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
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time, frac), Init_B));
        [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, Dose_Time, frac), tspan2, init_cond, options);
                
        T = [T_predose; T_postdose];
        Y = [Y_predose; Y_postdose];

        T_vals(i) = round(T(end)./365);


end

plot(abs(frac_vals - 1)*100, T_vals - disease_timing, 'LineWidth', 1.3, 'Color', 'k')
xlabel('Percent of APC Depleted', 'FontSize', 17)
ylabel('Difference in Time to 20% Beta Cell Mass', 'FontSize', 17)

ax = gca;
ax.FontSize = 21;

set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')


