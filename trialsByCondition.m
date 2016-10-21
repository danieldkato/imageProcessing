% Last updated DDK 2016-10-21

% OVERVIEW:
% This function takes activity traces from over the course of an entire
% acquisition session, along with trial information from that session, to
% parse activity into individual trials organized by trial condition. See
% OUTUPUTS below for a detailed description of the output format. 


% INPUTS:
% 1) activity - N x T activity matrix, where N is the number of activity
% channels (e.g. ROIs, units, etc.) and T is the number of time steps (e.g.
% frames, samples) from a given acquisition session. 

% 2) trials - S x 2 cell array, where S is the number of trials delivered
% during the acquisition session. The first column of each row is an
% integer sample number during which a given trial began, and the second
% column is a string describing the trial condition. 

% 3) conditionNames - C x 1 cell array, where C is the number of trial
% conditions presented throughout the course of the acquisition. Each entry
% of C is a string containing the name of a single trial condition. These
% should match the trial condition descriptions contained in the second
% column of trials. 

% 4) preStim - integer number of samples before trial onset that should be
% included in the plots for each trial.

% 5) postStim - integer number of samples after trial onset that should be
% included in the plots for each trial. 


% OUTPUTS:
% 1) T - C x 1 cell array, where C is the number of trial conditions
% presented throughout the course of the acquisition. 

% Each element of T is itself a 1 x 2 cell array. The first element is a
% string containing the name of a trial condition. The second element is an
% N x P x D matrix, where N is the number of activity channels, P is the
% number of samples around each trial onset to be plotted for each trial,
% and D is the number of trials for the given condition. 


%%
function [T] = trialsByCondition(activity, trials, conditionNames, preStim, postStim)

    numConditions = size(conditionNames, 1);
    T = cell(numConditions, 1);

    % For each condition...
    for i = 1:numConditions
        
        conditionDat = cell(1,2);
        
        indices = trials(find(cellfun(@(c) strcmp(c, conditionNames{i}), trials(2,:))),1); % ... get the trial start index for each trial of that condition... 
        
        % ... and for each trial start index, get a range of samples around
        % it (specified by preStim and postStim) for each channel (e.g.,
        % ROI, unit, etc.). This means that for every index in indices,
        % arrayfun will return an N x P matrix, where N is the number of
        % channels, and P is the number of samples used to plot each trial.
        % These will be assembled into planes, a 1 x I cell array where I
        % is the number of indices (i.e. trials of the current condition),
        % and each element of planes is an N x P matrix.
        planes = arrayfun(@(d) activity(:,d-preStim:d+postStim), indices, 'UniformOutput', 0); 
        cube = reshape(planes, 1, 1, size(planes,2)); % Reshape planes as an N x T x I cube
        cube = cell2mat(cube); % Convert cube from a cell array to a matrix
        
        conditionDat{1} = conditionNames{i};
        conditionDat{2} = cube;
        
        T{i} = conditionDat;
    end
end