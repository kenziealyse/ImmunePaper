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
tspan = 0:.1:100;

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

EIC = [100, 1000, 10000];
RIC = [0, 5000];


colors = ['r', 'k', 'b', 'm', 'c', 'y', 'g', "#7E2F8E", "#D95319"]; % Red, Green, Blue, Black
sim = 1;
entryIndex = 1;

figure(1)

for i = 1:length(RIC)


    for j = 1:length(EIC)

    AL_initcond = 0;
    EL_initcond = EIC(j);
    RL_initcond = RIC(i);
    AP_initcond = 0;
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


    figure(1)
    plot(EL./10e7, RL./10e7, 'Color', colors(sim), 'LineStyle', '--', 'HandleVisibility', 'on', 'LineWidth', 1.5)
    hold on
    plot(EL(end)./10e7, RL(end)./10e7, '*', 'Color', colors(sim), 'HandleVisibility', 'off', 'MarkerSize', 20, 'linewidth', 6)
    ylabel('R_L', 'FontSize', 16)
    xlabel('E_L', 'FontSize', 16)

    C = RIC(i)/EIC(j);

    Ess = kR/(C + kR/kE);
    Rss = -kR/kE * Ess + kR;

    plot(Ess./10e7, Rss./10e7, 'o', 'MarkerSize', 15, 'LineWidth', 2, 'Color', colors(sim))


    figure(2)
    plot(tspan, B./max(B), 'Color', colors(sim), 'LineStyle', '--')
    hold on

    legendEntries{entryIndex} = ['Full System, E_L = ', num2str(EIC(i)), '  R_L = ', num2str(EIC(j))];
    legendEntries{entryIndex + 1} = ['Not Full System SS, E_L = ', num2str(EIC(i)), '  R_L = ', num2str(EIC(j))];
    sim = sim + 1;
    entryIndex = entryIndex + 2;

    end
    

end

figure(1)
legend(legendEntries, 'Location', 'Best');
plot([kE./10e7 0], [0 kR./10e7], '--', 'Color', "#EDB120", 'HandleVisibility', 'off', 'linewidth', 2)
plot(kE./10e7, 0, "pentagram", 'Color', "#77AC30", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#77AC30", 'HandleVisibility', 'off')
hold on
plot(0, kR./10e7, "pentagram", 'Color', "#7E2F8E", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#7E2F8E", 'HandleVisibility', 'off')
