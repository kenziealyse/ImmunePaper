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
tspan = 0:.1:10000;

params2 = LoadParameters();

phiA = params2(1);
deltaA = params2(2);
lambdaEL = params2(3);
omegaEL = params2(4);
phiE = params2(5);
deltaE = params2(6);
lambdaR = params2(7);
omegaRL = params2(8);
C = params2(9);
phiR = params2(10);
deltaR = params2(11);
sigma = params2(12);
alpha = params2(13);
deltaP = params2(14);
kappa = params2(15);
r1 = params2(16);
r2 = params2(17);

kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
kR = 1/((omegaRL)/(C*(omegaRL - (phiR+deltaR))));
alphaE = omegaEL - (phiE+deltaE);
alphaR = omegaRL - (phiR+deltaR);

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';


params2 = [alphaE, kE, alphaR, kR]';

APIC = 0:200:1500;

colors = ['r', 'k', 'b', 'm', 'c', 'y', 'g', "#7E2F8E", "#D95319"]; % Red, Green, Blue, Black
sim = 1;
entryIndex = 1;

figure(1)

for i = 1:length(APIC)

    AL_initcond = 0;
    EL_initcond = 0;
    RL_initcond = 0;
    AP_initcond = APIC(i);
    EP_initcond = 0;
    RP_initcond = 0;
    B_initcond = 1*10^6;
    
    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                 EP_initcond RP_initcond B_initcond]';
    
    % Run the Model
    options = odeset('Events', @(t, Y) events(t, Y, RegTcellModel(t,Y, params)));
    [~,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond);

    EL = Y(:,2);
    RL = Y(:,3);
    B = Y(:,7);

    
    if RL(end) <= 10e-5
        RL(end) = 0;
    end
    if EL(end) < 10e-5
        EL(end) = 0;
    end

    plot(EL./10e7, RL./10e7, 'Color', colors(i), 'LineStyle', '--', 'HandleVisibility', 'on', 'LineWidth', 1.5)
    hold on
    ylabel('R_L', 'FontSize', 18)
    xlabel('E_L', 'FontSize', 18)  

    legendEntries{i} = ['Full System, AP =' num2str(APIC(i))];

end
    
    
legend(legendEntries, 'Location', 'NorthEast', 'FontSize', 18);
plot([kE./10e7 0], [0 kR./10e7], '--', 'Color', "#EDB120", 'HandleVisibility', 'off', 'linewidth', 2)
% plot(kE./10e7, 0, "pentagram", 'Color', "#77AC30", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#77AC30", 'HandleVisibility', 'off')
% hold on
% plot(0, kR./10e7, "pentagram", 'Color', "#7E2F8E", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#7E2F8E", 'HandleVisibility', 'off')
xlim([0 kE./10e7])
ylim([0 kR./10e7])

figure(2)
plot(tspan, RL./10e7, 'Color', colors(i), 'LineStyle', '--', 'HandleVisibility', 'on', 'LineWidth', 1.5)
