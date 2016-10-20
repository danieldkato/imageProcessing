% Last updated DDK 2016-10-18

% OVERVIEW:
% This function takes raw intensity traces put out by a segmentation step
% and returns dF/F traces in which the baseline F is based on a moxing
% boxcar average of frames around any given frame.


% INPUTS:
% 1) input - path to a .csv containing a K x T matrix of raw intensity
% values from a given grab, where K is the number of ROIs extracted from
% the grab and T is the number of frames in the grab. 

% 2) m - the number of frames before the current frame included in the
% moving boxcar average.

% 3) n - the number of frames after the current frame included in the
% moving boxcar average.


% OUTPUTS:
% 1) dFF - K x T matrix of dF/F matrix, where K is the number of ROIs and T is
% the number of frames in the input grab. Note that because the baseline
% for each frame is computed from m frames before and n frames after that
% frame, we cannot compute dF/F for the first m or last n frames. The
% columns corresponding to these frames are therefore replaced with NaNs.

% In addition to returning dFF formally, this function also saves dFF to
% disk as a K x (T-pre-post) .csv for use in subsequent processing steps.

%%
function dFF = movingBoxcarDFF(input, m, n, output)
    dat = csvread(input); % Load the data
    numROIs = size(dat, 1);
    numFrames = size(dat, 2);
    
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
    % an anonymous function @(i) mean(dat(:,i-m:i+2)). arrayfun will
    % sequentially iterate through includedFrames and substitute every
    % element for i.
    dFF = arrayfun(@(i) mean(dat(:,i-m:i+n),2), includedFrames, 'UniformOutput', 0);
    
    % Convert dFF from a cell array to a matrix (each iteration of @(i) by
    % arrayfun results in a K x 1 vector, and the output of each iteration
    % goes into its own cell, resulting in a 1 x F cell array where each
    % cell is a K x 1 vector).
    dFF = cell2mat(dFF);
    
    % Future processing steps may need to know which frames were omitted in
    % order to perform the moving boxcar averaging. We can do this by
    % padding dFF with m columns of NaNs before the dFF data and n columns
    % of NaNs after the dFF data.
    dFF = [NaN(numROIs, m), dFF, NaN(numROIs, n)];
    
    % Save the outputs to disk:
    status = exist(output, 'dir');
    if status == 0
        mkdir(output);
    end 
    
    old = cd(output);
    csvwrite('dFF.csv', dFF);
    
    %% Write metadata:
    
    inputs = {{'raw intensity traces', input}};
    outputs = {{'dF/F traces', output}};
    params = {{'pre', m};
              {'post', n}};
    
    writeMetadata('compute_dF/F', 'moving boxcar average', inputs, outputs, params);
    
    cd(old);
end