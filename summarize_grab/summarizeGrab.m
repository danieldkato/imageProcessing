% Last updated DDK 2016-10-23

% OVERVIEW: 
% This function creates summary plots for a given acqusition session during
% which multiple trial or stimulus conditions are presented. This
% includes:

% 1) An overview plot that includes full activity traces for every
% identified cell throughout the entire acquisition session. Also includes
% shaded rectangles denoting the time and condition of each stimulus
% period.

% 2) For every identified cell:
%   a) a plot of the mean activity with SEM bars for each condition, and 
%   b) a plot including activity traces or each individual trial, color
%   coded by condition


% REQUIREMENTS:
% 1) The MATLAB function trimExp, available at
% https://github.com/danieldkato/dFF/blob/master/trimExp.m


% INPUTS:
% 1) activity - path to a .csv containing an N x T activity matrix, where N
% is the number of identified cells (e.g. ROIs, units, etc.) and T is the
% number of time steps (e.g., frames, samples, etc.) in the given data
% acquisition session. 

% 2) trials - path to a .csv containing an S x 2 trial matrix, where S is
% the number of trials to be analyzed. Each row corresponds to a trial; the
% first column of each row is the frame start number of a given trial, and
% the second column is a string describing the trial condition.

% 3) conditions - a C x 1 cell array of structs, where C is the number of
% distinct trial conditions presented during throughout the course of the
% trials to be analyzed. Each struct must have at least the following three
% fields:

%   a) Name - the name of the trial condition. This must exactly match the
%   trial condition descriptions in the second column of trials, described
%   above.

%   b) Color - color code for the given condition, in any valid MATLAB
%   format for encoding color. Will be used in plotting. 

%   c) abbreviation - abbreviation for the given trial condition. Will be
%   used in creating legends for each figure.

% 4) preStim - amount of time before stimulus onset from which
% trial-by-trial data should be plotted, in seconds.

% 5) postStim - amount of time after stimulus onset from which
% trial-by-trial data should be plotted, in seconds.

% 6) outputDirectory - directory where all created figures should be saved.

% OUTPUTS:
% This function has no formal return, but saves the plots described above
% for the given data acquisition session.

%%
function summarizeGrab(activityPath, trialsPath, preStim, postStim, outputDirectory, grabMetadata, conditions)
    
    % Load activity data:
    activity = csvread(activityPath);
    numROIs = size(activity,1);
    
    % Load trial data and strip header:
    [n, t, trials] = xlsread(trialsPath);
    trials = trials(7:end, :);
    
    %{
    % Load condition settings:
    conditions = importdata(conditions);
    %}
    
    % Trim flanking NaN columns from boxcar-averaged dF/F data, discard any
    % trials that occur during frames corresponding to NaN columns, apply
    % offset to trial start frames so that they are relative to start of
    % trimmed data: 
    [activity, trials] = trimExp(activity, trials);
    
    %overviewPath = plotOverview(activity, trials, conditions, outputDirectory);
    [meanPaths, rawPaths] = plotPerCell(activity, trials, preStim, postStim, outputDirectory, grabMetadata, conditions);    
    
    %% Write metadata:
    
    inputs = {{'activty traces', activityPath};
              {'trial matrix', trialsPath};
              {'conditions structure', conditions}};
    params = {{'pre-stim period', preStim};
              {'post-stim period', postStim}};
    outputs = cell(1 + numROIs*2,1);
    outputs{1} = {'overview figure', overviewPath};
    
    roiIndices = (1:1:numROIs)';
    outputs(2:numROIs+1) = arrayfun(@(a) {strcat(['ROI ', num2str(a), ' mean traces figure']), meanPaths{a}} , roiIndices, 'UniformOutput', 0);
    outputs(numROIs+2:numROIs*2 + 1) = arrayfun(@(a) {strcat(['ROI ', num2str(a), ' raw traces figure']), rawPaths{a}} , roiIndices, 'UniformOutput', 0);
    %disp(outputs);
    
    writeMetadata('summarize_grab', 'summarize_grab', inputs, outputs, params);
end