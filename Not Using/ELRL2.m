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
sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
r1 = params(16);
r2 = params(17);

kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
kR = 1/((omegaRL)/(C*(omegaRL - (phiR+deltaR))));
alphaE = omegaEL - (phiE+deltaE);
alphaR = omegaRL - (phiR+deltaR);

figure(1)
plot(kE./10e7, 0, "pentagram", 'Color', "#77AC30", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#77AC30")
hold on
plot(0, kR./10e7, "pentagram", 'Color', "#7E2F8E", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#7E2F8E")
legend('k_E', 'k_R', 'FontSize', 14)
xlabel('EL', 'FontSize', 18)
ylabel('RL', 'FontSize', 18)
if kR > kE
    title('Not Full Equations, k_R > k_E', 'FontSize', 20)
elseif kR < kE
    title('Not Full Equations, k_R < k_E', 'FontSize', 20)
elseif kR == kE
    title('Not Full Equations, k_R = k_E', 'FontSize', 20)
end


params = [alphaE, kE, alphaR, kR]';

EL_initcondvals = 1:200:400;

for i = 1:length(EL_initcondvals)

    for j = 1:length(EL_initcondvals)

    EL_initcond = EL_initcondvals(i);
    RL_initcond = EL_initcondvals(j);
    
    init_cond = [EL_initcond RL_initcond]';
    
    % Run the Model
    [T,Y] = ode23s(@(t,Y) ELRL(t,Y, params), tspan, init_cond);
    
    % Relabel to easily keep track of compartments
    % LN
    EL = Y(:,1);
    RL = Y(:,2);

    figure(1)
    if RL(end) <= 10e-5
        RL(end) = 0;
    end
    if EL(end) < 10e-5
        EL(end) = 0;
    end

    figure(1)
%     plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 0.8, 'Color', 'r', 'MarkerSize', 21, 'HandleVisibility', 'off')
    plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 1, 'MarkerSize', 21, 'HandleVisibility', 'off')

hold on

%     figure(2)
%     plot(tspan, EL./10e7,'b')
%     hold on
%     plot(tspan, RL./10e7,'r', 'LineStyle', '--')
% 
    end
    

end

% if kR > kE
%     saveas(gcf,'Not Full Equations, k_R > k_E.png')
% elseif kR < kE
%     saveas(gcf,'Not Full Equations, k_R < k_E.png')
% elseif kR == kE
%     saveas(gcf,'Not Full Equations, k_R = k_E.png')
% end



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