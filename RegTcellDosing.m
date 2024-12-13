% CLEAR THE WORKSPACE
clc
clear
close all

% Dosing Information
Dose_vals = [1e5];%1e4:1000:1e5;
Dose_Time_vals = 0;%(0.1:1:60).*365;


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

T_vals = zeros(length(Dose_vals), length(Dose_Time_vals));

for i = 1:length(Dose_Time_vals)

    Dose_Time = Dose_Time_vals(i);
    tspan = 0:1:365*70;


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
        
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';
        
        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';
              
        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, Dose_NuRegTcellmodel(t,Y, params, Dose_amount, Dose_Time), B_initcond));
        [T,Y] = ode23s(@(t,Y)Dose_NuRegTcellmodel(t,Y, params, Dose_amount, Dose_Time), tspan, init_cond, options);
                
        T_vals(j,i) = T(end)./365;

        figure(1)
        plot(T./365, Y(:,3))
        hold on

    end

end

    % Define custom colormap
    figure(1)
    imagesc(Dose_Time_vals./365, Dose_vals, T_vals);
    xlim([0 max(Dose_Time_vals)./365])
    ylim([0 max(Dose_vals)])
    xlabel('Dose Times', 'FontSize', 17);
    ylabel('Dose Amounts', 'FontSize', 17);
    set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
    title('Time to 20% Beta Cell Mass');
    
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
    B = y(7);
%     disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - .2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end