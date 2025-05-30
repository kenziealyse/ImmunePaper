% Clear the workspace
clc
clear
close all

% Load the file
load('r1vsr2heatmap1e-05_vals.mat')

% % Heat Map Plot
% figure(1)
% imagesc(results.r1_vals, results.r2vals, results.TTD);
% hold on
% xlim([0 3*10^(-5)])
% ylim([0 3*10^(-5)])
% xlabel('r_1', 'FontSize', 17);
% ylabel('r_2', 'FontSize', 17);
% set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
% title('Time to 20% Beta Cell Mass');
% ax = gca;
% ax.FontSize = 21;

% Loop through the vector (starting from the second element)
[num_rows, num_columns] = size(results.TTD);

boundary_rows = [];
boundary_cols = [];

for i = 1:num_rows  % Loop through rows
    for j = 1:num_columns
        if results.TTD(i, j) == 70 && results.TTD(i, j-1) ~= 70  % Check if current value is 70 and previous value is not 70
            boundary_rows = [boundary_rows; i];  % Store row index
            boundary_cols = [boundary_cols, j];  % Store column index
        end
    end
end

r1vals_boundary = results.r1_vals(boundary_cols);
r2vals_boundary = results.r1_vals(boundary_rows);

plot(r1vals_boundary,r2vals_boundary, 'Color', 'r', 'LineWidth', 1.5)


close all

%%

% Load the file
load('highRL0r1vsr2heatmap1e-05_vals.mat')

% Heat Map Plot
figure(1)
imagesc(results.r1_vals, results.r2vals, results.TTD);
hold on
xlim([0 3*10^(-5)])
ylim([0 3*10^(-5)])
xlabel('r_1', 'FontSize', 17);
ylabel('r_2', 'FontSize', 17);
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time to 20% Beta Cell Mass');
ax = gca;
ax.FontSize = 21;

plot(r1vals_boundary,r2vals_boundary, 'Color', 'r', 'LineWidth', 1.5)

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

set(gcf, 'Position', [0 0 650 400])