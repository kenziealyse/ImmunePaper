% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:70*365;

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
r2 = params(17);

AL_initcond = 0;
EL_initcond = 10;
RL_initconds = 10:10000000:1*10^9;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

figurename = 'r1heatmapDynamicslownu.pdf';

r1 =   0.2*10^(-4);%[0.2*10^(-4), 2.2*10^(-4), 5.5*10^(-4)];
nu =   1*10^(-5);%0.01;%

T_vals = zeros(length(RL_initconds), 1);
lambda = zeros(length(RL_initconds), 1);

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';


for i = 1:length(RL_initconds)

    RL_initcond = RL_initconds(i);
          

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';


        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);
        
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

        T_vals(i) = T(end)./365;

    

end

% disp(T_vals)

plot(RL_initconds, T_vals, 'LineWidth', 1.3, 'Color', 'k')
xlabel('RL(0)', 'FontSize', 17)
ylabel('Time to Disease, years', 'FontSize', 17)

function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end