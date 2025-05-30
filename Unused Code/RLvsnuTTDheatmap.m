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

N = 50;
nu_vals =  linspace(0.25*10^(-5), 2500*10^(-5), N);
RL_initconds = [10, 1e2, 1e3, 1e4, 1e5, 1e6];

figurename = 'Figures\timetodiseasevsnuheatmap.pdf';

AL_initcond = 0;
EL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

tspan1 = 0:0.01:70*365;

% Preallocate Space
time = zeros(length(RL_initconds), length(nu_vals));

for i = 1:length(RL_initconds)

    RL_initcond = RL_initconds(i)

    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
                EP_initcond RP_initcond B_initcond]';

    for j = 1:length(nu_vals)

        nu = nu_vals(j);

        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                      phiR deltaR kappa r1 r2 nu]';

        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

        % Save end time
        time(i,j) = T(end)./365;

    end

end

% Define custom colormap
figure(1)
imagesc(nu_vals, RL_initconds, time);
hold on
% set(gca, 'Yscale', 'log')
% ylim([10 1e6])
xlabel('\nu', 'FontSize', 17);
ylabel('R_L(0)', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time to 20% Beta Cell Mass', 'FontSize', 17);
ax = gca;
ax.FontSize = 21;

% Define the number of colors
numColors = 100;

% Define transitions for red, green, and blue
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)]; % Red fades to zero
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)]; % Green grows
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)]; % Blue in middle

% Combine into a custom colormap
customColorMap = [redValues', greenValues', blueValues']; 

% Apply the colormap
colormap(gca, customColorMap);

% Show color bar with custom tick labels
colorbar;
clim([0, 70]);  % Set color axis limits based on data range

set(gcf, 'Position', [100, 300, 800, 500]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder