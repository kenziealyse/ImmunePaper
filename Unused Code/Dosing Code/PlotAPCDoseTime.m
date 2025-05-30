% CLEAR THE WORKSPACE
clc
clear
close all

% Load simulations
disease_type = 'early';
nu = 0.01;

load(['ComboTherapyFiles/minus30Days_', num2str(nu), '_', disease_type ,'_.mat']);
load(['ComboTherapyFiles/plus30Days_', num2str(nu), '_', disease_type ,'_.mat']);
load(['ComboTherapyFiles/Original_', num2str(nu), '_', disease_type ,'_.mat']);

%% Minus 30 Days Plot
% Define custom colormap
subplot(1,2,1)
imagesc(Original.DoseTime, Original.DoseVals, Original.TTD - minus30Days.TTD);
xlim([min(Original.DoseTime) max(Original.DoseTime)])
ylim([min(Original.DoseVals) max(Original.DoseVals)])
xlabel('Treg Dose Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title({'30 days prior', 'Difference in Time to 20% Beta Cell Mass', '(Original - Modified)'}, 'FontSize', 12);
ax = gca;
ax.FontSize = 21;
yticks([1e1 1e3 1e5 1e7 1e9])
set(gca, 'Yscale', 'log')
set(gcf, 'Color', 'White')
set(gcf, 'Position', [100 300 800 500])

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

% Set z axis limits
% min(min(Original.TTD - minus30Days.TTD))
% zlimit = ceil(max(max(Original.TTD - minus30Days.TTD))/10)*10/2;
% if zlimit == 0
%     zlimit = 10;
% end
clim([-30, 60]);  % Set color axis limits based on data range

% clim([-5, 5]);  % Set color axis limits based on data range

%% Plus 30 day plot
% Define custom colormap
subplot(1,2,2)
imagesc(Original.DoseTime, Original.DoseVals, Original.TTD - plus30Days.TTD);
xlim([min(Original.DoseTime) max(Original.DoseTime)])
ylim([min(Original.DoseVals) max(Original.DoseVals)])
xlabel('Treg Dose Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title({'30 days after', 'Difference in Time to 20% Beta Cell Mass', '(Original - Modified)'}, 'FontSize', 12);
ax = gca;
ax.FontSize = 21;
yticks([1e1 1e3 1e5 1e7 1e9])
set(gca, 'Yscale', 'log')
set(gcf, 'Color', 'White')

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

% % Set z axis limits
% zlimit = ceil(max(max(Original.TTD - plus30Days.TTD))/10)*10/2;
% if zlimit == 0
%     zlimit = 10;
% end
% clim([-zlimit, zlimit]);  % Set color axis limits based on data range
% % clim([-5, 5]);  % Set color axis limits based on data range
clim([-30, 60]);  % Set color axis limits based on data range


% Set Plot Size
set(gcf, 'Position', [100, 300, 1200, 500]);
set(gcf, 'Color', 'White')
% 
% % Save the figure as pdf
% set(gcf, 'Units', 'Inches');
% pos = get(gcf, 'Position');
% set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
% saveas(gcf, figurename); % Save Figure in Folder


%% Minus 30 Days Plot
% Define custom colormap
figure()
subplot(1,2,1)
imagesc(Original.DoseTime, Original.DoseVals, Original.TTD);
xlim([min(Original.DoseTime) max(Original.DoseTime)])
ylim([min(Original.DoseVals) max(Original.DoseVals)])
xlabel('Treg Dose Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title({'Original'}, 'FontSize', 12);
ax = gca;
ax.FontSize = 21;
yticks([1e1 1e3 1e5 1e7 1e9])
set(gca, 'Yscale', 'log')
set(gcf, 'Color', 'White')
set(gcf, 'Position', [100 300 800 500])

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

% Set z axis limits
% zlimit = ceil(max(max(Original.TTD))/(10/2))*10/2;
% if zlimit == 0
%     zlimit = 10;
% end
% clim([0, zlimit]);  % Set color axis limits based on data range


subplot(1,2,2)
imagesc(minus30Days.DoseTime, minus30Days.DoseVals, minus30Days.TTD);
xlim([min(minus30Days.DoseTime) max(minus30Days.DoseTime)])
ylim([min(minus30Days.DoseVals) max(minus30Days.DoseVals)])
xlabel('Treg Dose Time, years', 'FontSize', 17);
ylabel('Dose Amount', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title({'30 days prior'}, 'FontSize', 12);
ax = gca;
ax.FontSize = 21;
yticks([1e1 1e3 1e5 1e7 1e9])
set(gca, 'Yscale', 'log')
set(gcf, 'Color', 'White')
set(gcf, 'Position', [100 300 800 500])

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

% Set z axis limits
% zlimit = ceil(max(max(Original.TTD))/(10/2))*10/2;
% if zlimit == 0
%     zlimit = 10;
% end
% clim([0, zlimit]);  % Set color axis limits based on data range
