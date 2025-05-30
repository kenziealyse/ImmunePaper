function index = firstBelowThreshold(vector, threshold)
    % Find the indices where the vector elements are below the threshold
    belowThresholdIndices = find(vector < threshold);
    
    % Check if any indices were found
    if isempty(belowThresholdIndices)
        index = -1;  % No element is below the threshold
    else
        % Get the first index where the element is below the threshold
        index = belowThresholdIndices(1);
    end
end