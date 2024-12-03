 % CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:27375;

params = LoadParameters();

alphaER = params(1);
kC = params(2);
deltaER = params(3);
Psi = params(4);
omegaEA = params(5);
phiEA = params(6);
deltaEA = params(7);
phiA = params(8);
deltaA = params(9);
alphaR = params(10);
omegaR = params(11);
phiR = params(12);
deltaR = params(13);
alphaA = params(14);
d = params(15);
r1 = params(16);
r2 = params(17);

kC_vals = 5:5:100;

for i = 1:length(kC_vals)

    params = [alphaER kC_vals(i) deltaER Psi omegaEA phiEA deltaEA phiA deltaA alphaR ...
        omegaR phiR deltaR alphaA d r1 r2]';
    
    B_init_cond = 1.06*10^6;
    Tcell_init_cond = 6*10^8;
    DC_init_cond = 17;

    init_cond = [Tcell_init_cond 0 0 0 0 0 DC_init_cond 0 B_init_cond]';
    
    threshold = 0.2.*B_init_cond;
    
    % Run the Model
    [T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond);
    
    % Relabel to easily keep track of compartments
    % LN
    ER = Y(:,1);
    C = Y(:,2);
    EAL = Y(:,3);
    AL = Y(:,4);
    RL = Y(:,5);
    
    %Pancreas
    EAP = Y(:,6);
    AP = Y(:,7);
    RP = Y(:,8);
    B = Y(:,9);
    
    index = firstBelowThreshold(B, threshold);
    if index == -1
        time(i) = T(end)./365;
    else
        time(i) = T(index);
    end

end
% Define custom colormap
imagesc(alphaER, kC_vals, time);
xlabel('alphaER');
ylabel('kC');
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Time (years) to 20% Beta Cell Mass');

% Define custom colormap
numColors = 100;
greenValues = linspace(0, 1, numColors);  % Green values range from 0 to 1
redValues = linspace(1, 0, numColors);    % Red values range from 1 to 0

% Create custom colormap without blue component
customColorMap = [redValues', greenValues', zeros(numColors, 1)];  % Red and green, blue is zero

% Apply custom colormap to the current figure
colormap(customColorMap);

% Show color bar with custom tick labels
colorbar;
caxis([0, max(time(:))]);  % Set color axis limits based on data range
% colorbarTicks = linspace(min(time(:)), max(time(:)), 5);  % Define color bar tick positions
% colorbarLabels = cellstr(num2str(colorbarTicks(:), '%.2f'));  % Create color bar tick labels
% colorbar('YTick', colorbarTicks, 'YTickLabel', colorbarLabels);
% 
