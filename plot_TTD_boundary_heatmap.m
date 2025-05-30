%% README: Boundary Detection and Heatmap Visualization
% This script loads two MATLAB data files containing simulation results
% of TTD (Time to disease) as a function of parameters r1 and r2.
%
% Part 1: Detects boundary points where TTD reaches 70 for the first time along each row.
% Part 2: Plots a heatmap of TTD values and overlays the detected boundary.
%
% Custom colormap is applied for better visualization.
% Author: [Your Name]
% Date: [Today's Date]

%% Clear workspace and load first data file
clc
clear
close all

% Load the first result file (used for boundary detection)
load('Mat Files/r1vsr2heatmap1e-05_vals.mat')

% Get dimensions of the TTD matrix
[num_rows, num_columns] = size(results.TTD);

% Initialize arrays to store boundary row and column indices
boundary_rows = [];
boundary_cols = [];

% Loop through the TTD matrix to detect where value changes to 70
for i = 1:num_rows
    for j = 2:num_columns  % Start from 2 to safely access (j-1)
        if results.TTD(i, j) == 70 && results.TTD(i, j-1) ~= 70
            boundary_rows = [boundary_rows; i];
            boundary_cols = [boundary_cols, j];
        end
    end
end

% Convert boundary indices to corresponding r1 and r2 values
r1vals_boundary = results.r1_vals(boundary_cols);
r2vals_boundary = results.r1_vals(boundary_rows);

% Optional quick plot of the boundary (commented out close all to keep figure)
plot(r1vals_boundary, r2vals_boundary, 'Color', 'r', 'LineWidth', 1.5)

close all

%% Load second file and plot heatmap
load('Mat Files/highRL0r1vsr2heatmap1e-05_vals.mat')

figure(1)
imagesc(results.r1_vals, results.r2vals, results.TTD);
hold on

% Set axes limits and labels
xlim([0 3e-5])
ylim([0 3e-5])
xlabel('r_1', 'FontSize', 17)
ylabel('r_2', 'FontSize', 17)
set(gca, 'YDir', 'normal')  % Standard axis orientation
title('Time to 20% Beta Cell Mass')
ax = gca;
ax.FontSize = 21;

% Overlay the detected boundary on the heatmap
plot(r1vals_boundary, r2vals_boundary, 'Color', 'r', 'LineWidth', 1.5)

%% Create and apply custom colormap
numColors = 100;

% Create color transitions
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];

% Assemble into colormap
customColorMap = [redValues', greenValues', blueValues'];

% Apply to current axes
colormap(gca, customColorMap);

% Display colorbar and set limits
colorbar;
clim([0, 70]);

% Resize figure
set(gcf, 'Position', [0 0 650 400])
