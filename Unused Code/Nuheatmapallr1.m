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
r2 = params(17);

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

indvplots = 0;
indvplots_compare = 1;
heatmapplots = 0;

figurename = 'r1heatmapDynamicslownu.pdf';

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';


N = 300; % Desired number of points
r1_vals =  5.5*10^(-4);%[0.2*10^(-4), 2.2*10^(-4), 5.5*10^(-4)];%linspace(0, 0.6*10^(-3), N); %
nu_vals =   0.01;%1*10^(-5);%linspace(2.5*10^(-6), 0.05, N); %
colors = {'k'};%{'r', 'b', 'g'};

T_vals = zeros(length(nu_vals), length(r1_vals));
lambda = zeros(length(nu_vals), length(r1_vals));


for i = 1:length(r1_vals)

    r1 = r1_vals(i)

    for j = 1:length(nu_vals)
        
        nu = nu_vals(j);

        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
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

        lambda(j,i) = kappa*max(EP)/(1+r2*min(RP));

        if indvplots == 1
            figure(5)
            subplot(2,4,1)
            plot(T./365,AL, 'k', 'LineWidth', 1.3)   
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e7])
            xlim([0 10])
            ylabel('AL', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,2)
            plot(T./365,EL, 'k', 'LineWidth', 1.3)        
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            xlim([0 10])
            ylim([10e-1 10e7])
            ylabel('EL', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,3)
            plot(T./365,RL, 'k', 'LineWidth', 1.3)        
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e7])
            xlim([0 10])
            ylabel('RL', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,5)
            plot(T./365,AP, 'k', 'LineWidth', 1.3)        
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e7])
            xlim([0 10])
            ylabel('AP', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,6)
            plot(T./365,EP, 'k', 'LineWidth', 1.3)        
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e7])
            xlim([0 10])
            ylabel('EP', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,7)
            plot(T./365,RP, 'k', 'LineWidth', 1.3)        
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e7])
            xlim([0 10])
            ylabel('RP', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
            subplot(2,4,8)
            plot(T./365,B, 'k', 'LineWidth', 1.3) 
            hold on
            yline(.2*B_initcond, '--', 'Color', 'k', 'LineWidth', 1.3);
            set(gca, 'Xscale', 'log')
            set(gca, 'Yscale', 'log')
            xlim([0 100])
            ylim([10^5 10e5])
            % Define custom y-ticks
            yticks([10^2, 10^4, 0.2 * B_initcond, 10^6]);
            % Define custom y-tick labels
            yticklabels({'10^2', '10^4', '20% B_0', '10^6'});
            ylabel('B', 'FontSize', 17)
            xlabel('Time, years', 'FontSize', 17)
        end
        if indvplots_compare == 1
            colors{i}
            figure(5)
            subplot(2,4,1)
            plot(T./365,AL, 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('A_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,2)
            plot(T./365,EL, 'Color', colors{i}, 'LineWidth', 1.3)  
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('E_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            % legend('Early Disease Onset', 'Late Disease Onset', 'No Disease', 'FontSize', 13)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,3)
            plot(T./365,RL, 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('R_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,5)
            plot(T./365,AP, 'Color', colors{i}, 'LineWidth', 1.3)
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('A_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,6)
            plot(T./365,EP, 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('E_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,7)
            plot(T./365,RP, 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('R_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([0 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])
            subplot(2,4,8)
            plot(T./365,B, 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            yline(.2*B_initcond, '--', 'Color', 'k', 'LineWidth', 1.3);
            set(gca, 'Xscale', 'log')
            ylabel('B, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            xlim([0 10e2])
            ylim([10^5 10e5])
            % Define custom y-ticks
            yticks([10^2, 10^4, 0.2 * B_initcond, 10^6]);
            % Define custom y-tick labels
            yticklabels({'10^2', '10^4', '20% B_0', '10^6'});

            figure(6)
            kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
            kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
            plot(EL(1:10),RL(1:10), 'Color', colors{i}, 'LineWidth', 1.3) 
            hold on
            scatter(EL(end),RL(end), 'Color', colors{i}, 'LineWidth', 1.3, 'Marker', '*') 
            plot([kE 0], [0 kR], '--', 'Color', "#EDB120", 'HandleVisibility', 'off', 'linewidth', 2)
            % set(gca, 'Xscale', 'log')
            xlim([0 0.03])
            xline(EL(end)./10e7)
            ylabel('R_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            % set(gca, 'Yscale', 'log')
            % ylim([0 10e8])
            % xlim([0 10e2])
            % yticks([10e-1 10e1 10e3 10e5 10e7])
        end

    end  

    

end


disp(T_vals)

if heatmapplots == 1
    % Define custom colormap
    figure(1)
    subplot(1,2,1)
    imagesc(r1_vals, nu_vals, T_vals);
    hold on
    xlim([0 max(r1_vals)])
    ylim([0 max(nu_vals)])
    xlabel('r_1 values', 'FontSize', 17);
    ylabel('\nu values', 'FontSize', 17);
    set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
    title('Time to 20% Beta Cell Mass');
    
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
    caxis([0, 70]);  % Set color axis limits based on data range

    % Points to plot
    earlydiseasepoint = [0.2*10^(-4), 0.01];
    latediseasepoint = [2.2*10^(-4), 0.01];
    nodiseasepoint = [5.5*10^(-4), 0.01];

    % scatter(earlydiseasepoint(1), earlydiseasepoint(2), 'filled', 'MarkerEdgeColor', 'k',...
    %     'MarkerFaceColor', 'k', 'LineWidth',1.5);
    % scatter(latediseasepoint(1), latediseasepoint(2), 'filled', 'MarkerEdgeColor', 'k', ...
    %     'MarkerFaceColor', 'k', 'LineWidth',1.5)
    % scatter(nodiseasepoint(1), nodiseasepoint(2), 'filled', 'MarkerEdgeColor', 'k', ...
    %     'MarkerFaceColor', 'k', 'LineWidth',1.5)
    % 

    % rectangle('Position', [0, 0, max(r1_vals), 0.5*10^(-3)], 'EdgeColor', 'k', 'LineWidth', 2);

    figure(1)
    subplot(1,2,2)
    imagesc(r1_vals, nu_vals, T_vals);
    hold on
    xlim([0 max(r1_vals)])
    ylim([0 0.5*10^(-3)])
    xlabel('r_1 values', 'FontSize', 17);
    ylabel('\nu values', 'FontSize', 17);
    set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
    title('Time to 20% Beta Cell Mass, Smaller \nu');

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
    caxis([0, 70]);  % Set color axis limits based on data range

    % Points to plot
    nudiseasepoint = [5.5*10^(-4), 1*10^(-5)];

    % scatter(nudiseasepoint(1), nudiseasepoint(2), 'filled', 'MarkerEdgeColor', 'k',...
    %     'MarkerFaceColor', 'k', 'LineWidth',1.5);
    % 
end

set(gcf, 'Position', [100, 300, 1400, 500]);
    % set(gcf, 'Position', [100, 300, 800, 500]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end