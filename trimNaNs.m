% Last updated DDK 2016-10-20

% OVERVIEW:
% This function trims any columns of NaNs from a K x T activity matrix for
% a given grab, where K is the number of ROIs in the grab and T is the
% number of frames in the grab. 

% This trimming is sometimes necessary because certain methods for
% computing dF/F (e.g., using a moving boxcar average to compute the
% baseline F) necessitate discarding a number of frames at the beginning
% and/or end of the grab. In order to preserve informaton about which
% frames were discarded, I've written these dF/F functions so that the
% columns to be discarded are substituted with NaNs rather than omitted
% wholesale. This makes it possible for downstream functions (like this
% one) to pass along which frames of the raw data are actually being
% exlcuded. This is useful for registering dF/F traces to other
% experimental data (like trial type information) that have already been
% registered to raw activity traces.

% Note that this function is only equipped to deal with data where NaNs are
% confined to the first m or last n consecutive columns, where m and n are
% integers >=0. In other words, this function does not deal with columns of
% NaNs in the middle of the data.


% INPUTS: 
% 1) tracesIn - K x T activity matrix for a given grab, where K is the
% number of ROIs in the grab and T is the number of frames in the original
% grab. The first m and last n columns of tracesIn may be populated with
% NaNs, where m and n are integers >=0.


% OUTPUTS:
% 1) tracesOut - K x F activity matrix where K is the number of ROIs in the
% grab and F is the number of non-NaN columns in tracesIn.

% 2) startOffset - index of the first non-NaN column of tracesIn.

% 3) stopOffset - index of the last non-NaN column of tracesIn. 

%%
function [tracesOut, startOffset, stopOffset] = trimNaNs(tracesIn)
   
    whereNans = (isnan(tracesIn));
    tracesOut = tracesIn(~whereNans);

    row = sum(whereNans,1); 
    shiftRowRight = [0, row(1:end-1)];
    shiftRowLeft = [row(2:end), 0];
    
    startOffset = find(shiftRowRight>row);
    stopOffset = find(shiftRowLeft>row);
end