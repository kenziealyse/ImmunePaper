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

figure(1)
plot(kE./10e7, 0, "pentagram", 'Color', "#77AC30", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#77AC30")
hold on
plot(0, kR./10e7, "pentagram", 'Color', "#7E2F8E", 'LineWidth', 2.5, 'MarkerSize', 18, 'MarkerFaceColor', "#7E2F8E")
legend('k_E', 'k_R', 'FontSize', 14)
xlabel('E_L', 'FontSize', 18)
ylabel('R_L', 'FontSize', 18)
if kR > kE
    title('Full Equations, k_R > k_E', 'FontSize', 20)
elseif kR < kE
    title('Full Equations, k_R < k_E', 'FontSize', 20)
elseif kR == kE
    title('Full Equations, k_R = k_E', 'FontSize', 20)
end



params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaRL C ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';

EL_initcondvals = 1:100:4000;

for i = 1:length(EL_initcondvals)

    for j = 1:length(EL_initcondvals)

    AL_initcond = 2000;
    EL_initcond = EL_initcondvals(i);
    RL_initcond = EL_initcondvals(j);
    AP_initcond = 0;
    EP_initcond = 0;
    RP_initcond = 0;
    B_initcond = 1*10^6;
    
    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                 EP_initcond RP_initcond B_initcond]';
    
    % Run the Model
    options = odeset('Events', @(t, Y) events(t, Y, RegTcellModel(t,Y, params)));
    [T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond);
    
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

    if RL(end) <= 10e-5
        RL(end) = 0;
    end
    if EL(end) < 10e-5
        EL(end) = 0;
    end

    if EL_initcondvals(i) > EL_initcondvals(j)
        figure(1)
        plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 0.8, 'Color', 'm', 'MarkerSize', 21, 'HandleVisibility', 'off')
        hold on
    elseif EL_initcondvals(i) < EL_initcondvals(j) 
        figure(1)
        plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 0.8, 'Color', 'c', 'MarkerSize', 21, 'HandleVisibility', 'off')
        hold on
    elseif EL_initcondvals(i) == EL_initcondvals(j) 
        figure(1)
        plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 0.8, 'Color', 'g', 'MarkerSize', 21, 'HandleVisibility', 'off')
        hold on
    end

%     figure(1)
%     plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 0.8, 'Color', 'r', 'MarkerSize', 21, 'HandleVisibility', 'off')
%     hold on
% 
%     figure(2)
%     plot(tspan, EL./10e7,'b')
%     hold on
%     plot(tspan, RL./10e7,'r', 'LineStyle', '--')

    end
    

end

% if kR > kE
%     saveas(gcf,'Full Equations, k_R > k_E.png')
% elseif kR < kE
%     saveas(gcf,'Full Equations, k_R < k_E.png')
% elseif kR == kE
%     saveas(gcf,'Full Equations, k_R = k_E.png')
% end
% 


% ylim([0, kR./10e7])
% xlim([0, kE./10e7])
