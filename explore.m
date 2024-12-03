%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:1000;

params = LoadParameters();

phiA = params(1);
deltaA = params(2);
lambdaEL = params(3);
omegaEL = params(4);
phiE = params(5);
deltaE = params(6);
lambdaR = params(7);
omegaRL = params(8);
C = params(9);
phiR = params(10);
deltaR = params(11);
% sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
% r1 = params(16);
r2 = 0;

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 100;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

sigma_vals = 0;%0:1:1000;
r1_vals = 1;%[10e-10, 10e-6, 10e-2, 10e-1];
    
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                 EP_initcond RP_initcond B_initcond]';

kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
kR = 1/((omegaRL)/(C*(omegaRL - (phiR+deltaR))));
alphaE = omegaEL - (phiE+deltaE);
alphaR = omegaRL - (phiR+deltaR);
    
lamda = zeros(length(sigma_vals), 1);

colors = ['r', 'g', 'b', 'm'];

for j = 1:length(r1_vals)

    r1 = r1_vals(j);

    for i = 1:length(sigma_vals)
    
        sigma = sigma_vals(i);
    
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
              phiR deltaR sigma alpha deltaP kappa r1 r2]';
        
        % Run the Model
        [T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond);
    
        EL = Y(:,2);
        RL = Y(:,3);
        EP = Y(:,5);
        RP = Y(:,6);
        B = Y(:,7);
        
        if RL(end) <= 10e-5
            RL(end) = 0;
        end
        if EL(end) < 10e-5
            EL(end) = 0;
        end
    
        top = kappa*(phiE/deltaE*(-kE/kR*RL(end) + kE));
        bottom = 1 + r2*(phiR/deltaR*RL(end));
        lamda(i) = top/bottom;
        
%         figure(2) 
%         plot(EL(end)./10e7, RL(end)./10e7, '*')
%         hold on
%         plot([kE./10e7 0], [0 kR./10e7], '--', 'Color', "#EDB120", 'HandleVisibility', 'off', 'linewidth', 2)
%     

figure(2)
plot(T./365,B./max(B))
hold on
set(gca, 'Xscale', 'log')
    end
    
    figure(1)
    subplot(1,2,1)
    plot(sigma_vals, lamda, 'LineWidth', 1.3, 'Color', colors(j))
    hold on
    legendEntries{j} = ['r_1 = ' num2str(r1_vals(j))];
end

title(['EL(0) = ', num2str(EL_initcond), ', RL(0) = ', num2str(RL_initcond), '\newline r2 = ', num2str(r2)])
ylabel('\lambda_{RL}', 'FontSize', 17)
xlabel('\sigma', 'FontSize', 17)
set(gca, 'Xscale', 'log')
legend(legendEntries, 'Location', 'NorthEast', 'FontSize', 18);





%% r2

r2_vals = [10e-10, 10e-6, 10e-2, 10e-1];
r1 = 0;


for j = 1:length(r2_vals)

    r2 = r2_vals(j);

    for i = 1:length(sigma_vals)
    
        sigma = sigma_vals(i);
    
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
              phiR deltaR sigma alpha deltaP kappa r1 r2]';
        
        % Run the Model
        [~,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond);
    
        EL = Y(:,2);
        RL = Y(:,3);
        EP = Y(:,5);
        RP = Y(:,6);
        B = Y(:,7);
        
        if RL(end) <= 10e-5
            RL(end) = 0;
        end
        if EL(end) < 10e-5
            EL(end) = 0;
        end
    
        top = kappa*(phiE/deltaE*(-kE/kR*RL(end) + kE));
        bottom = 1 + r2*(phiR/deltaR*RL(end));
        lamda(i) = top/bottom;
     
    end
    
    figure(1)
    subplot(1,2,2)
    plot(sigma_vals, lamda, 'LineWidth', 1.3, 'Color', colors(j))
    hold on
    legendEntries{j} = ['r_2 = ' num2str(r1_vals(j))];
end

title(['EL(0) = ', num2str(EL_initcond), ', RL(0) = ', num2str(RL_initcond), '\newline r1 = ', num2str(r1)])
ylabel('\lambda_{RL}', 'FontSize', 17)
xlabel('\sigma', 'FontSize', 17)
set(gca, 'Xscale', 'log')
legend(legendEntries, 'Location', 'NorthEast', 'FontSize', 18);

