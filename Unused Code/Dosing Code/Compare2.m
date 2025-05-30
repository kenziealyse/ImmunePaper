% CLEAR THE WORKSPACE
clc
clear
close all

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
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);
Dose_Time = 1*365;%5*365;
Dose_Amount = 0;%1e6;
timing = 100*365;
B_initial = 1*10^6;

N = 100;
nu_vals =  linspace(1*10^(-5), 2500*10^(-5), N); 

figurename = 'timetodiseasevsnu.pdf';

tspan1 = 0:0.01:Dose_Time;
tspan2 = Dose_Time:0.01:70*365;


% Preallocate Space
time = zeros(length(nu_vals), 1);
  
for i = 1:length(nu_vals)

    nu = nu_vals(i);

    AL_initcond = 0;
    EL_initcond = 10;
    RL_initcond = 1e6;
    AP_initcond = 0;
    EP_initcond = 0;
    RP_initcond = 0;
    B_initcond = 1*10^6;

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
    EP_initcond RP_initcond B_initcond]';


    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params, timing), B_initial));
    [T1,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params, timing), tspan1, init_cond, options);

    % Relabel Compartments to Easily Keep Track
    AL = Y(:,1);
    EL = Y(:,2);
    RL = Y(:,3);
    AP = Y(:,4);
    EP = Y(:,5);
    RP = Y(:,6);
    B = Y(:,7);

    AL_initcond = AL(end);
    EL_initcond = EL(end);
    RL_initcond = RL(end) + Dose_Amount;
    AP_initcond = AP(end);
    EP_initcond = EP(end);
    RP_initcond = RP(end);
    B_initcond =  B(end);

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
    EP_initcond RP_initcond B_initcond]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params, timing), B_initial));
    [T2,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params, timing), tspan2, init_cond, options);

    T = [T1;T2];

    % Save end time
    time(i) = T(end)./365;

end

figure(1)
plot(nu_vals, time, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-')
hold on

%% Deplete APC
% Preallocate Space
time = zeros(length(nu_vals), 1);
timing = 5*365;
  
for i = 1:length(nu_vals)

    nu = nu_vals(i);

    AL_initcond = 0;
    EL_initcond = 10;
    RL_initcond = 1e6;
    AP_initcond = 0;
    EP_initcond = 0;
    RP_initcond = 0;
    B_initcond = 1*10^6;

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
    EP_initcond RP_initcond B_initcond]';


    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params, timing), B_initial));
    [T1,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params, timing), tspan1, init_cond, options);

    % Relabel Compartments to Easily Keep Track
    AL = Y(:,1);
    EL = Y(:,2);
    RL = Y(:,3);
    AP = Y(:,4);
    EP = Y(:,5);
    RP = Y(:,6);
    B = Y(:,7);

    AL_initcond = AL(end);
    EL_initcond = EL(end);
    RL_initcond = RL(end) + Dose_Amount;
    AP_initcond = AP(end);
    EP_initcond = EP(end);
    RP_initcond = RP(end);
    B_initcond =  B(end);

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
    EP_initcond RP_initcond B_initcond]';
    
    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params, timing), B_initial));
    [T2,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params, timing), tspan2, init_cond, options);

    T = [T1;T2];

    % Save end time
    time(i) = T(end)./365;

end

figure(1)
plot(nu_vals, time, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')
hold on
% xline(0.01)
ax = gca;
ax.FontSize = 21;
ylabel('Time to 20% Beta Cell Mass')
xlabel('\nu', 'FontSize', 30)

legend('R_L(0) = 1 \times 10^6', ['R_L(0) = 1 \times 10^6, APC depletion timing = ', num2str(timing/365)], 'Location', 'best')

set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder