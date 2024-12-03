% CLEAR THE WORKSPACE
clc
clear
close all

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


EL = 0:10:C;

lambdaEL = (kappa.*phiE*EL)./(deltaE.*(1 + r2.*((kR./kE).*EL + kR)));

RL = 0:10:C;

lambdaRL = (kappa.*phiE.*((kE./kR).*RL + kE))./(deltaE.*(1 + r2.*(phiR./deltaR).*RL));

if kE > kR
    'kE bigger'
elseif kE < kR
    'kR bigger'
end


figure();
subplot(1,2,1)
plot(EL, lambdaEL, 'linewidth', 1.5, 'Color', 'b')
xlabel('E_L', 'FontSize', 18)
ylabel('\lambda_{EL}', 'FontSize', 18)
set(gca, 'Xscale', 'log')
ylim([min(min(lambdaEL, lambdaRL)) max(max(lambdaEL, lambdaRL))])
xlim([0 C])

subplot(1,2,2)
plot(RL, lambdaRL, 'linewidth', 1.5, 'Color', 'r')
xlabel('R_L', 'FontSize', 18)
ylabel('\lambda_{RL}', 'FontSize', 18)
set(gca, 'Xscale', 'log')
ylim([min(min(lambdaEL, lambdaRL)) max(max(lambdaEL, lambdaRL))])
xlim([0 C])

sgtitle(['r_2 = ', num2str(r2)])




