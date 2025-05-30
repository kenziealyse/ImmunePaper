% CLEAR THE WORKSPACE
clc
clear
close all

% User change these
disease_type = 'late';
nu = 1e-5;
day_vals = [-30, 0, 30];
disease_timing = 34.6806;
end_time = disease_timing*365; % Only dose up to time of disease

for k = 1:length(day_vals)

   days = day_vals(k);
    % Dosing Information
    N = 50;
    Dose_vals = [10.^(0:1:10)];
    Dose_Time_vals = linspace(1*365, end_time, N);
    
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
        APC_DoseTime = Dose_Time_vals(i) + days;
    
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
            options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, APC_DoseTime), Init_B));
            [T_predose,Y_predose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, APC_DoseTime), tspan1, init_cond, options);
    
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
            options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel_withdepletion(t,Y, params, APC_DoseTime), Init_B));
            [T_postdose,Y_postdose] = ode23s(@(t,Y) NuRegTcellmodel_withdepletion(t,Y, params, APC_DoseTime), tspan2, init_cond, options);
                    
            T = [T_predose; T_postdose];
            Y = [Y_predose; Y_postdose];
    
            T_vals(j,i) = round(T(end)./365);
    
    
        end
    
    end
    
    if days == -30
        minus30Days.DoseTime = Dose_Time_vals./365;
        minus30Days.DoseVals = Dose_vals;
        minus30Days.TTD = T_vals - disease_timing;
        save(['ComboTherapyFiles/minus30Days_', num2str(nu), '_', disease_type ,'_.mat'], 'minus30Days');
    elseif days == 30
        plus30Days.DoseTime = Dose_Time_vals./365;
        plus30Days.DoseVals = Dose_vals;
        plus30Days.TTD = T_vals - disease_timing;
        save(['ComboTherapyFiles/plus30Days_', num2str(nu), '_', disease_type ,'_.mat'], 'plus30Days');
    elseif days == 0
        Original.DoseTime = Dose_Time_vals./365;
        Original.DoseVals = Dose_vals;
        Original.TTD = T_vals - disease_timing;
        save(['ComboTherapyFiles/Original_', num2str(nu), '_', disease_type ,'_.mat'], 'Original');
    else
        Test.DoseTime = Dose_Time_vals./365;
        Test.DoseVals = Dose_vals;
        Test.TTD = T_vals - disease_timing;
        save(['ComboTherapyFiles/Test_', num2str(nu), '_', disease_type ,'_.mat'], 'minus30Days');
    end


end