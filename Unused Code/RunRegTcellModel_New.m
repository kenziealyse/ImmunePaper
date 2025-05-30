% CLEAR THE WORKSPACE
clc
clear
close all

% Dosing Information
Dose_vals = 0;%0:1e8:1e10;
Dose_Time_vals = 70*366;%(1:.1:6).*365;

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

r1 = earlydiseasepoint(1);
r2 = earlydiseasepoint(2);
nu = 1*10^(-5);%0.01;%

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

T_vals = zeros(length(Dose_vals), length(Dose_Time_vals));

for i = 1:length(Dose_Time_vals)

    Dose_Time = Dose_Time_vals(i);


    for j = 1:length(Dose_vals)

        Dose_amount = Dose_vals(j);

        %% Pre Dose
        
        AL_initcond = 0;
        RLA_initcond = 1e9;
        RL_initcond = 10;
        EL_initcond = 10;
        AP_initcond = 0;
        EP_initcond = 0;
        RP_initcond = 0;
        B_initcond = 1*10^6;

        init_cond = [AL_initcond RLA_initcond RL_initcond EL_initcond ...
            AP_initcond EP_initcond RP_initcond B_initcond]';

        tspan1 = 0:0.001:Dose_Time;
              
        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, RegTcellDosingModel_New(t,Y, params, Dose_amount, Dose_Time), Init_B));
        [T_predose,Y_predose] = ode23s(@(t,Y) RegTcellDosingModel_New(t,Y, params, Dose_amount, Dose_Time), tspan1, init_cond, options);
        
        % T_predose(end)
        % Dose_Time(end)
        if T_predose(end)/365 < Dose_Time/365

            T_vals(j,i) = round(T_predose(end)./365); 
           
        else 
            AL_initcond = Y_predose(end, 1);
            RLA_initcond = Y_predose(end, 2);
            RL_initcond = Y_predose(end, 3);
            EL_initcond = Y_predose(end, 4);
            AP_initcond = Y_predose(end, 5);
            EP_initcond = Y_predose(end, 6);
            RP_initcond = Y_predose(end, 7);
            B_initcond = Y_predose(end, 8);
    
            init_cond = [AL_initcond RLA_initcond RL_initcond EL_initcond ...
                AP_initcond EP_initcond RP_initcond B_initcond]';
    
            tspan2 = Dose_Time:.11:70*365;
    
            % Run the Model
            options = odeset('Events', @(t, Y) events(t, Y, RegTcellDosingModel_New(t,Y, params, Dose_amount, Dose_Time), Init_B));
            [T_postdose,Y_postdose] = ode23s(@(t,Y) RegTcellDosingModel_New(t,Y, params, Dose_amount, Dose_Time), tspan2, init_cond, options);
                    
            T = [T_predose; T_postdose];
            Y = [Y_predose; Y_postdose];
    
            T_vals(j,i) = round(T(end)./365)

            figure(2)
            subplot(2,4,1)
            plot(T./365, Y(:,1), 'LineWidth', 1.3)
            % set(gca, 'Yscale', 'log')
            set(gca, 'Xscale', 'log')
            % ylim([0 max(Y(:,2))])
            xlim([0 10])
            ylabel('A_L, cells')

            subplot(2,4,2)
            plot(T./365, Y(:,2), 'LineWidth', 1.3)
            ylabel('R_{L,A}, cells')

            subplot(2,4,3)
            plot(T./365, Y(:,3), 'LineWidth', 1.3)
            ylabel('R_{L}, cells')

            subplot(2,4,4)
            plot(T./365, Y(:,4), 'LineWidth', 1.3)
            ylabel('E_{L}, cells')

            subplot(2,4,5)
            plot(T./365, Y(:,5), 'LineWidth', 1.3)
            ylabel('A_P, cells')

            subplot(2,4,6)
            plot(T./365, Y(:,6), 'LineWidth', 1.3)
            ylabel('E_P, cells')

            subplot(2,4,7)
            plot(T./365, Y(:,7), 'LineWidth', 1.3)
            ylabel('R_P, cells')


            subplot(2,4,8)
            plot(T./365, Y(:,8), 'LineWidth', 1.3)
            hold on
            yline(Init_B*0.2, '--')
            ylabel('\beta, cells')

        end


    end

end

    % Define custom colormap
    figure(1)
    imagesc(Dose_Time_vals./365, Dose_vals, T_vals - T_vals(1,1));
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
    caxis([0, 10]);  % Set color axis limits based on data range

function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(8);
%     disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - .2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end