% CLEAR THE WORKSPACE
clc
clear
close all

figurename = 'LogScaler1vsnuTTD';

% Load Parameters
load('NuHeatmapr1_vals.mat')
r1_vals = results.r1vals;
nu_vals = results.nuvals;
T_vals = results.TTD;


% Heat Map Plot
figure(1)
imagesc(r1_vals, nu_vals, T_vals);
hold on
xlim([0 max(r1_vals)])
ylim([0 0.5*10^(-5)])
set(gca, 'Yscale', 'log')
xlabel('r_1', 'FontSize', 17);
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