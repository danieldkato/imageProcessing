% Last updated DDK 2016-10-18

% OVERVIEW:
% This function takes raw intensity traces put out by a segmentation step
% and returns dF/F traces in which the baseline F is based on a moxing
% boxcar average of frames around any given frame. Processed data is
% returned in 2 formats:

% 1) a grand k x (t - m -n) activity matrix, where k is the number of ROIs
% analyzed from the current grab, t is the number of frames analyzed from
% the current grab, m is the number of frames before any given frame to
% include in the baseline average for that frame, and n is the number of
% frames after any given frame to include in the baseline average for that
% frame.

% 2) a c x 1 cell array containing activity parsed into individual trials
% and sorted by condition, where c is the number of trial conditions to be
% analyzed. Each element of the cell array is a k x p x u slab of data
% corresponding to a single trial condition, where k is the number of ROIs
% to be analyzed, p is the number of samples in the peri-stimulus period to
% be analyzed for each trial, and u is the number of trials of the given
% trial condition.


% REQUIREMENTS:
% 1) The MATLAB function trialsByCondition.m, avaialable at 
% https://github.com/danieldkato/trial_registration/blob/master/trialsByCondition.m

% Note that this function treats the peri-stimulus period as a fixed number
% of frames, which is really only sensible if the stimulus period of every
% trial is the same.


% INPUTS:
% 1) input - path to a .csv containing a K x T matrix of raw intensity
% values from a given grab, where K is the number of ROIs extracted from
% the grab and T is the number of frames in the grab. 

% 2) m - the number of frames before the current frame included in the
% moving boxcar average.

% 3) n - the number of frames after the current frame included in the
% moving boxcar average.

% 4) Trials - an s x 3 cell array, where s is the number of trials to be
% analyzed. The first column of each row gives the sample start number of a
% given trial, the second column gives the name of the trial condition, and
% the third column gives the trial duration. Note that this fucntion
% assumes that Trials includes a 7-row header. 

% 5) preStim - the number of samples before stimulus onset to include in
% trial-by-trial plots.

% 6) postStim - the number of samples after stimulus onset to include in
% trial-by-trial plots. 

% 7) outputDir - directory to which to save outputs.

% 8) condSettingsPath - path to a .txt file contaning information about
% different trial conditions used in the current experiment. This text file
% should contain a series of MATLAB statements defining a c x 1 cell array
% of structs called Conditions, where c is the number of trial conditions
% to be analyzed from the current grab. Each struct should minimally
% include fields corresponding to whatever trial parameters are necessary
% for defining a condition in the current experiment, as well as a concise
% condition name that will be assigned to each trial and that can be used
% later on to easily parse experiments by condition. Example contents of a
% conditionSettings.txt file might include:

% Conditions{1}.STPRIDX = 1;
% Conditions{1}.SPKRIDX = 0;
% Conditions{1].Name = ' stepper only ';

% Conditions{2}.STPRIDX = 0;
% Conditions{2}.SPKRIDX = 1;
% Conditions{2}.Name = ' speaker only ';


% OUTPUTS:
% 1) dFF - a grand k x (t - m -n) activity matrix, where k is the number of ROIs
% analyzed from the current grab, t is the number of frames analyzed from
% the current grab, m is the number of frames before any given frame to
% include in the baseline average for that frame, and n is the number of
% frames after any given frame to include in the baseline average for that
% frame.

% 2) TBC - a c x 1 cell array containing activity parsed into individual trials
% and sorted by condition, where c is the number of trial conditions to be
% analyzed. Each element of the cell array is a k x p x u slab of data
% corresponding to a single trial condition, where k is the number of ROIs
% to be analyzed, p is the number of samples in the peri-stimulus period to
% be analyzed for each trial, and u is the number of trials of the given
% trial condition.

% In addition to returning dFF and TBC formally, this function also saves
% dFF to storage as a .csv and an HDF5, and saves TBC as an HDF5.

%%
function [dFF, TBC] = movingBoxcarDFF(input, m, n, Trials, preStim, postStim, outputDir, condSettingsPath)
    dat = csvread(input); % Load the data
    numROIs = size(dat, 1);
    numFrames = size(dat, 2);
    
    % Load trial data:
    [num, txt, trials] = xlsread(Trials);
    
    %% Compute whole-grab dFF trace:
    
    % Create an F x 1 vector of frames, where F is the number of frames for
    % which dF/F will be computed. Recall that for a moving boxcar average,
    % this must be less than the total number of frames in the grab K;
    % since the baseline for any given frame is based on m frames before
    % and n frames after that frame, there must be at least m frames before
    % and n frames after every frame for which we are to compute dF/F;
    % thus, we omit the first m and last n frames.
    includedFrames = (m+1:1:numFrames-n); 
    
    % Call arrayfun to apply a boxcar averaging function to every frame in
    % includedFrames. We can define the boxcar averaging function inline as
    % an anonymous function @(i) dat(:,i)/mean(dat(:,i-m:i+2)). arrayfun will
    % sequentially iterate through includedFrames and substitute every
    % element for i.
    dFF = arrayfun(@(i) (dat(:,i)-mean(dat(:,i-m:i+n),2)) ./ mean(dat(:,i-m:i+n),2), includedFrames, 'UniformOutput', 0);
    
    % Convert dFF from a cell array to a matrix (each iteration of @(i) by
    % arrayfun results in a K x 1 vector, and the output of each iteration
    % goes into its own cell, resulting in a 1 x F cell array where each
    % cell is a K x 1 vector).
    dFF = cell2mat(dFF);
    
    % Future processing steps may need to know which frames were omitted in
    % order to perform the moving boxcar averaging. We can do this by
    % padding dFF with m columns of NaNs before the dFF data and n columns
    % of NaNs after the dFF data.
    dFFpadded = [NaN(numROIs, m), dFF, NaN(numROIs, n)];
    
    
    %% Parse the data into individual trials and organize by condition:
    
    trials(:,1) = cellfun(@(c) c-m, trials(:,1), 'UniformOutput', 0); % <-- REMEMBER!! we've removed the first m frames of the data, so trial start indices must be offset appropriately
    TBC = trialsByCondition(dFF, trials, preStim, postStim, condSettingsPath); 
    
    
    %% Save the outputs to disk:
    status = exist(outputDir, 'dir');
    if status == 0
        mkdir(outputDir);
    end 

    old = cd(outputDir);
    
    % Create and save files containing full trace dataset:
    csvwrite('dFF.csv', dFF); % create a .csv for easy manual inspection
    h5create('dFF_full_traces.h5', '/fullActivityTraces', size(dFFpadded)); % create an HDF5 in case this is useful later on
    h5write('dFF_full_traces.h5', '/fullActivityTraces', dFFpadded);

    disp(TBC(:,1));
    
    % Create a separate HDF5 file for data parsed into individual trials, organized into one dataset per condition:
    for c = 1:length(TBC)
        dSetName = strcat(['/', TBC{c,1}]);
        h5create('dFF_parsed_trials.h5', dSetName, size(TBC{c,2}));
        h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_pre_stim', preStim);
        h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_post_stim', postStim);
        h5write('dFF_parsed_trials.h5', dSetName, TBC{c,2});
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Color', Conditions{c}.Color);
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Abbreviation', Conditions{c}.Abbreviation);
    end

    
    %% Write metadata:
    
    inputs = {{'raw intensity traces', input}};
    outputs = {{'dF/F traces', outputDir}};
    params = {{'pre', m};
              {'post', n}};
    
    writeMetadata('compute_dF/F', 'moving boxcar average', inputs, outputs, params);
    
    cd(old);
end