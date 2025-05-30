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

params2 = [alphaE, kE, alphaR, kR]';


RIC = [1000, 10000, 100000];
EIC = [100000, 100, 10000000];



colors = ['r', 'k', 'b', 'm', 'c', 'y', 'g', "#7E2F8E", "#D95319"]; % Red, Green, Blue, Black
sim = 1;
entryIndex = 1;

figure(1)

for i = 1:length(RIC)


    for j = 1:length(EIC)

    EL_initcond = EIC(j);
    RL_initcond = RIC(i);

    init_cond = [EL_initcond RL_initcond]';
    
    % Run the Model
    [T,Y] = ode23s(@(t,Y) ELRL(t,Y, params2), tspan, init_cond);
    
    % Relabel to easily keep track of compartments
    % LN
    EL2 = Y(:,1);
    RL2 = Y(:,2);
    

    figure(1)
    if RL2(end) <= 10e-5
        RL2(end) = 0;
    end
    if EL2(end) < 10e-5
        EL2(end) = 0;
    end

    plot(EL2./10e7, RL2./10e7, 'Color', colors(sim), 'LineStyle', '-',  'HandleVisibility', 'on', 'LineWidth', 1.5)
    hold on
    
    C = RIC(i)/EIC(j);

    Ess = kR/(C + kR/kE);
    Rss = -kR/kE * Ess + kR;

    plot(Ess./10e7, Rss./10e7, 'o', 'MarkerSize', 15, 'LineWidth', 2, 'Color', colors(sim))

    legendEntries{entryIndex} = ['Simulation E_L = ', num2str(EIC(j)), '  R_L = ', num2str(RIC(i))];
    legendEntries{entryIndex + 1} = ['S.S. Prediction E_L = ', num2str(EIC(j)), '  R_L = ', num2str(RIC(i))];
    sim = sim + 1;
    entryIndex = entryIndex + 2;



    end
end

legend(legendEntries, 'Location', 'Best');
plot([kE./10e7 0], [0 kR./10e7], '--', 'Color', "#EDB120", 'HandleVisibility', 'off', 'linewidth', 2)
hold on
% plot(kE./10e7, 0, "pentagram", 'Color', "#77AC30", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#77AC30", 'HandleVisibility', 'off')
% plot(0, kR./10e7, "pentagram", 'Color', "#7E2F8E", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#7E2F8E", 'HandleVisibility', 'off')
ylim([0 kR./10e7])
xlim([0 kE./10e7])


function DyDt = ELRL(~, Y, params)

% Paramter Values
alphaE = params(1);
kE = params(2);
alphaR = params(3);
kR = params(4);

% Relabel Compartments to Easily Keep Track
EL = Y(1);
RL = Y(2);

% ODE Equations
dELdt = alphaE*EL*(1 - (EL + RL)/kE);
dRLdt = alphaR*RL*(1 - (EL + RL)/kR);

%Output
DyDt = [dELdt; dRLdt];

end
