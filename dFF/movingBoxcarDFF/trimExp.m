% Last updated DDK 2016-20-21

% OVERVIEW: 
% This function trims any flanking NaN columns in a K x T matrix of
% activity traces from a given acquisition, where K is the number of
% channels (e.g. ROIs, units, etc.) and T is the number of time steps (e.g.
% frames, samples, etc.). In addition, it also makes the appropriate
% changes to the corresponding trial start sample numbers for any trials
% that have been registered to the original (i.e. untrimmed) activity
% traces.

% For example, using a boxcar average to compute the baseline F in dF/F
% necessitates discarding the first m and last n frames of the grab, where
% m is the number of frames before any given frame to be included in the
% baseline average for that frame, and n is the number of frames after that
% frame to be inlcuded in the baseline average for that frame.

% Moreover, since the initial trial registration step registers trial start
% times to frame numbers of the raw data, the trial start frames computed
% by the initial trial registration must be adjusted to be relative to the
% beginning of the trimmed, rather than the raw, data. This may necessitate
% discarding some trials, as some may have been delivered before the first
% usable frame.

% While this function is primarily intended for preprocessing
% boxcar-averaged dF/F data and corresponding trials, it can be used for
% any matrix of activity traces where the first and last few sample columns
% are NaNs.


% INPUTS:
% 1) activityIn - K x S activity matrix, where K is the number of activity
% channels (ROIs, units, etc.) and S is the number of time steps (samples,
% frames, etc.) in the given acquisition (grab, recording, etc).

% 2) trialsIn - T x 2 cell array, where T is the number of trials delivered
% in the current acquisition. The first column of each row is an integer
% sample number of the raw data at which a given trial starts, and the
% second column is a string describing the trial type.


% OUTPUTS: 
% 1) activityOut - K x R activity matrix, where R is the number of samples
% between flanking NaN sample columns. Identical to activityIn, except all
% flanking NaN columns have been discarded.

% 2) trialsOut - U x 2 cell array, where U is the number of trials
% delivered over the course of usable samples. Identical to trialsIn except
% trials delivered outside of usable samples have been omitted, and sample
% start numbers of retained trials have are expressed relative to the
% beginning of the trimmed, rather than the raw, data. 

%%
function [activityOut, trialsOut] = trimExp(activityIn, trialsIn)
    
    [activityOut, startOffset, stopOffset] = trimNaNs(activityIn); % Trim NaNs from activity matrix
    trialsIn(:,1) = num2cell(cellfun(@(c) c-startOffset, trialsIn(:,1))); % Substract startOffset from all trials start frames to adjust for the fact that we may be chopping off the first few frames
    
    % Omit any trials that occur before the first usable frame or after the last usable frame:
    trials2keep = cellfun(@(c) c>0 & c<size(activityOut, 2), trialsIn(:,1));
    trialsOut = trialsIn(find(trials2keep),:);
    
end