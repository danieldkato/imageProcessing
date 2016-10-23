function summarizeGrab(activity, trials, conditions, preStim, postStim, outputDirectory)
    
    % Load activity data:
    activity = csvread(activity);
    
    % Load trial data and strip header:
    [n, t, trials] = xlsread(trials);
    trials = trials(7:end, :);
    
    % Load condition settings:
    conditions = importdata(conditions);
    
    % Trim flanking NaN columns from boxcar-averaged dF/F data, discard any
    % trials that occur during frames corresponding to NaN columns, apply
    % offset to trial start frames so that they are relative to start of
    % trimmed data: 
    [activity, trials] = trimExp(activity, trials);
    
    %plotOverview(activity, trials, conditions, outputDirectory);
    plotPerCell(activity, trials, conditions, preStim, postStim, outputDirectory);
    
end