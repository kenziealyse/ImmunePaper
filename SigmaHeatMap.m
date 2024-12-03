% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:1:70*365;

params = LoadParameters();

phiA = params(1); 
deltaA = params(2);
lambdaEL = params(3);
omegaEL = params(4);
phiE = params(5);
deltaE = params(6);
lambdaR = params(7);
omegaR = params(8);
epsilon = params(9);
phiR = params(10);
deltaR = params(11);
% sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
% r1 = params(16);
r2 = params(17);

AL_initcond = 0;
EL_initcond = 1000;
RL_initcond = 1000;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;


init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

r1_vals = 0:0.01:0.5;
sigma_vals = 0:.01:5;

T_vals = zeros(length(sigma_vals), length(r1_vals));

for i = 1:length(r1_vals)

    r1 = r1_vals(i)

    for j = 1:length(sigma_vals)
        
        sigma = sigma_vals(j);

        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR epsilon ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';

        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, RegTcellModel(t,Y, params), B_initcond));
        [T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond, options);
        
        % Relabel to easily keep track of compartments
        % LN
        AL = Y(:,1);
        EL = Y(:,2);
        RL = Y(:,3);
        
        % Pancreas
        AP = Y(:,4);
        EP = Y(:,5);
        RP = Y(:,6);
        B = Y(:,7); 

        T_vals(j,i) = T(end)./365;

%         figure(2)
%         plot(T./365,B./max(B), 'Color', 'g')
%         hold on
%         yline(.2*B_initcond./max(B), '--', 'Color', 'k')
%         set(gca, 'Xscale', 'log')

    end  

end

% Define custom colormap
figure(1)
imagesc(r1_vals, sigma_vals, T_vals);
% xlim([0 max(r1_vals)])
% ylim([0 max(sigma_vals)])
xlabel('r_1 values');
ylabel('\sigma values');
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time to 20% Beta Cell Mass');

% Define custom colormap
numColors = 100;
greenValues = linspace(0, 1, numColors);  % Green values range from 0 to 1
redValues = linspace(1, 0, numColors);    % Red values range from 1 to 0

% Create custom colormap without blue component
customColorMap = [redValues', greenValues', zeros(numColors, 1)];  % Red and green, blue is zero

% Apply custom colormap to the current figure
colormap(customColorMap);

% Show color bar with custom tick labels
colorbar;
caxis([0, 70]);  % Set color axis limits based on data range


% %% R2
% 
% % % CLEAR THE WORKSPACE
% % clc
% % clear
% % close all
% 
% % Set Model Conditions
% tspan = 0:1:70*365;
% 
% params = LoadParameters();
% 
% phiA = params(1); 
% deltaA = params(2);
% lambdaEL = params(3);
% omegaEL = params(4);
% phiE = params(5);
% deltaE = params(6);
% lambdaR = params(7);
% omegaR = params(8);
% epsilon = params(9);
% phiR = params(10);
% deltaR = params(11);
% % sigma = params(12);
% alpha = params(13);
% deltaP = params(14);
% kappa = params(15);
% r1 = params(16);
% % r2 = params(17);
% 
% AL_initcond = 0;
% EL_initcond = 1000;
% RL_initcond = 1000;
% AP_initcond = 0;
% EP_initcond = 0;
% RP_initcond = 0;
% B_initcond = 1*10^6;
% 
% 
% init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
%          EP_initcond RP_initcond B_initcond]';
% 
% r2_vals = [10e-7 10e-6 10e-5 10e-4 10e-3 10e-2 10e-1];
% sigma_vals = 0:1:1000;
% 
% T_vals = zeros(length(r2_vals), length(sigma_vals));
% 
% for i = 1:length(r2_vals)
% 
%     r2 = r2_vals(i)
% 
%     for j = 1:length(sigma_vals)
%         
%         sigma = sigma_vals(j);
% 
%         params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR epsilon ...
%           phiR deltaR sigma alpha deltaP kappa r1 r2]';
% 
%         % Run the Model
%         options = odeset('Events', @(t, Y) events(t, Y, RegTcellModel(t,Y, params), B_initcond));
%         [T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond, options);
%         
%         % Relabel to easily keep track of compartments
%         % LN
%         AL = Y(:,1);
%         EL = Y(:,2);
%         RL = Y(:,3);
%         
%         % Pancreas
%         AP = Y(:,4);
%         EP = Y(:,5);
%         RP = Y(:,6);
%         B = Y(:,7); 
% 
%         T_vals(i,j) = T(end)./365;
% 
% %         figure(2)
% %         plot(T./365,B./max(B))
% %         hold on
% %         yline(.2*B_initcond./max(B), '--', 'Color', 'k')
% %         set(gca, 'Xscale', 'log')
% 
%     end  
% 
% end
% 
% % Define custom colormap
% figure(1)
% imagesc(r2_vals, sigma_vals, T_vals);
% xlim([0 max(r2_vals)])
% ylim([0 max(sigma_vals)])
% xlabel('r_2 values');
% ylabel('\sigma values');
% set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
% title('Time to 20% Beta Cell Mass');
% 
% % Define custom colormap
% numColors = 100;
% greenValues = linspace(0, 1, numColors);  % Green values range from 0 to 1
% redValues = linspace(1, 0, numColors);    % Red values range from 1 to 0
% 
% % Create custom colormap without blue component
% customColorMap = [redValues', greenValues', zeros(numColors, 1)];  % Red and green, blue is zero
% 
% % Apply custom colormap to the current figure
% colormap(customColorMap);
% 
% % Show color bar with custom tick labels
% colorbar;
% caxis([0, 70]);  % Set color axis limits based on data range
% 
function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
%     disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - .2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end