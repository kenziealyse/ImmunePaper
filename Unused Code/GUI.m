% Clear the workspace
clc
clear
close all



    % Create figure window
    f = figure('Position', [100, 100, 400, 300]);

    % Dropdown for selecting number of vectors
    uicontrol('Style', 'text', 'Position', [50, 220, 200, 20], 'String', 'Select Number of Vectors');
    numVectors = uicontrol('Style', 'popupmenu', 'String', {'1', '2'}, ...
                           'Position', [250, 220, 100, 30]);

    % Create text for lower limit input
    lower_limit_label = uicontrol('Style', 'text', 'Position', [50, 160, 150, 30]);
    
    % Create lower limit input
    lower_limit = uicontrol('Style', 'edit', 'Position', [200, 160, 100, 30]);

    % Create text for upper limit input
    upper_limit_label = uicontrol('Style', 'text', 'Position', [50, 120, 150, 30]);

    % Create upper limit input
    upper_limit = uicontrol('Style', 'edit', 'Position', [200, 120, 100, 30]);

    % Create text for number of points input
    num_points_label = uicontrol('Style', 'text', 'Position', [50, 80, 150, 30]);

    % Create number of points input
    num_points = uicontrol('Style', 'edit', 'Position', [200, 80, 100, 30]);

    % Store handles in the figure's userdata
    handles.numVectors = numVectors;
    handles.lower_limit = lower_limit;
    handles.upper_limit = upper_limit;
    handles.num_points = num_points;
    handles.lower_limit_label = lower_limit_label;
    handles.upper_limit_label = upper_limit_label;
    handles.num_points_label = num_points_label;
    guidata(f, handles);

    % Callback for when number of vectors is selected
    set(numVectors, 'Callback', @(src, event) updateUI(f));

    % Create button to generate the vectors
    uicontrol('Style', 'pushbutton', 'String', 'Generate Vectors', ...
              'Position', [150, 20, 100, 40], 'Callback', @(src, event) generateVectors(f));

    function updateUI(f)
        handles = guidata(f);
        numVec = str2double(handles.numVectors.String{handles.numVectors.Value});

        if numVec == 1
            set(handles.lower_limit_label, 'String', 'Lower Limit for r1');
            set(handles.upper_limit_label, 'String', 'Upper Limit for r1');
        else
            set(handles.lower_limit_label, 'String', 'Lower Limit for r1');
            set(handles.upper_limit_label, 'String', 'Upper Limit for r2');
        end
    end

    function generateVectors(f)
    handles = guidata(f);

    % Get values from the UI
    numVec = str2double(handles.numVectors.String{handles.numVectors.Value});
    lower_val_r1 = str2double(handles.lower_limit.String);
    upper_val_r1 = str2double(handles.upper_limit.String);
    points = str2double(handles.num_points.String);

    if numVec == 1
        % Generate 1 vector (r1)
        r1 = linspace(lower_val_r1, upper_val_r1, points);
        disp('Generated r1:');
        disp(r1);
    else
        % Generate 2 vectors: r1 and r2
        lower_val_r2 = str2double(handles.lower_limit.String);  % Use different limits if needed
        upper_val_r2 = str2double(handles.upper_limit.String);  % Use different limits if needed
        r1 = linspace(lower_val_r1, upper_val_r1, points);
        r2 = linspace(lower_val_r2, upper_val_r2, points);  % Generate second vector
        disp('Generated r1 and r2:');
        disp(r1);
        disp(r2);
    end
end

