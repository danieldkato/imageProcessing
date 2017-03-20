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
% conditions presented throughout the course of the acquisition. Each
% element of C is a struct with two fields: 'Name' and 'Color', where
% 'Name' is the name of the trial condition and 'Color' is its color code
% in relevant figures.

% 4) preStim - integer number of samples before trial onset that should be
% included in the plots for each trial.

% 5) postStim - integer number of samples after trial onset that should be
% included in the plots for each trial. 


% OUTPUTS:
% 1) T - C x 2 cell array, where C is the number of trial conditions
% presented throughout the course of the acquisition. 

% Each row of T corresponds to a trial condition. The first column of each
% row contains a string describing the name for a given condition. The
% second column of each row contains an N x P x D matrix, where N is the
% number of activity channels, P is the number of samples around each trial
% onset to be plotted for each trial, and D is the number of trials for the
% given condition.


%%
function TBC = trialsByCondition(activity, trials, preStim, postStim, Conditions)

    numConditions = length(Conditions);
    disp('trials(:,2)');
    disp(trials(:,2));
    
    numSamples = size(activity, 2);
    
    TBC = cell(numConditions, 2);

    % For each condition...
    for i = 1:numConditions
        
        indices = cell2mat(trials(find(cellfun(@(c) strcmp(c, Conditions{i}.Name), trials(:,2))),1)); % ... get the trial start index for each trial of that condition... 
        indices = indices(find( indices>preStim & indices<=numSamples-postStim )); % Omit trials that begin within preStim samples of the the beginning of the acquisition or within postStim samples of the end of the acquisition.  
        
        % ... and for each trial start index, get a range of samples around
        % it (specified by preStim and postStim) for each channel (e.g.,
        % ROI, unit, etc.). This means that for every index in indices,
        % arrayfun will return an N x P matrix, where N is the number of
        % channels, and P is the number of samples used to plot each trial.
        % These will be assembled into planes, an I x 1 cell array where I
        % is the number of indices (i.e. trials of the current condition),
        % and each element of planes is an N x P matrix.
        planes = arrayfun(@(d) activity(:,d-preStim:d+postStim), indices, 'UniformOutput', 0); 
        cube = reshape(planes, 1, 1, size(planes,1)); % Reshape planes as an N x T x I cube
        cube = cell2mat(cube); % Convert cube from a cell array to a matrix
        
        TBC{i,1} = Conditions{i}.Name;
        TBC{i,2} = cube;
    end
end