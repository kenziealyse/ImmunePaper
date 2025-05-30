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
r2 = 0;

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

figurename = 'r1vsNuheatmapDynamicslogscale.pdf';

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

% Set Vectors for r1 and nu values
N = 100; 
r1_vals =  linspace(0, 0.6*10^(-3), N);
nu_vals =  linspace(0.25*10^(-5), 10*10^(-5), N); 
yticks([2*10^(-5) 4*10^(-5) 6*10^(-5) 8*10^(-5) 10*10^(-5)])
% Preallocate Space
T_vals = zeros(length(nu_vals), length(r1_vals));

for i = 1:length(r1_vals)

    r1 = r1_vals(i)

    for j = 1:length(nu_vals)
        
        nu = nu_vals(j);

        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

        % Run the Model
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
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

        T_vals(j,i) = T(end)./365;

    end  
end

% Heat Map Plot
figure(1)
imagesc(r1_vals, nu_vals, T_vals);
hold on
xlim([0 max(r1_vals)])
ylim([0 max(nu_vals)])
yticks([1e-6 1e-5 1e-4 1e-3])
xlabel('r_1', 'FontSize', 17);
set(gca, 'Yscale', 'log')
ylabel('\nu', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time to 20% Beta Cell Mass');
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