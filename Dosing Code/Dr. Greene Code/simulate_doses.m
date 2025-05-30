function [Y, T] = simulate_doses(t0, tf, Y0, dose_times, dose_amounts, params)
% Dosing T_regs in the lymph node

num_time_steps = 100;  % Number of time steps between doses to solve ODE on
init_B = Y0(end);      % Used in event function to stop simulation
options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), init_B));

dose_times_with_init = [t0; dose_times];     % Add t0 to solve before doses
dose_amounts_with_init = [0; dose_amounts];  % No dose at t0 (unless dose_times start with t0)
Y = Y0;
T = t0;

for i = 1:length(dose_times)
    time_pre = dose_times_with_init(i);
    time_post = dose_times_with_init(i+1);
    Y_pre = Y(end, :);
    
    % if first dose is at time t0, then we don't need to solve any ODE initially
    if time_pre == time_post
        Y(end, 3) = Y(end, 3) + dose_amounts_with_init(i+1);
        continue;
    end

    [T_new, Y_new] = ode23s(@(t,Y) NuRegTcellmodel(t, Y, params), linspace(time_pre, time_post, num_time_steps), Y_pre, options);
    Y = [Y; Y_new];
    T = [T; T_new];
    % Break out of loop if event triggered
    if T(end) < time_post
        break;
    end
    % 3rd index in Y is RL.  Add dose before next loop
    Y(end, 3) = Y(end, 3) + dose_amounts_with_init(i+1);
end


% Solve until tf now that done with doses
% Only do if no breaks have occurred (gotten to end of doses with event triggered)
if T(end) == dose_times(end)
    [T_final, Y_final] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), linspace(T(end), tf, num_time_steps), Y(end,:), options);
    Y = [Y; Y_final];
    T = [T; T_final];
end

end