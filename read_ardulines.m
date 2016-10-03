% Last updated DDK 2016-09-27    

% OVERVIEW: 
% This function takes a .txt file containing serial output read from an
% Arduino running a given ArduFSM protocol for a single behavior session,
% and returns a T x 1 cell array of the trial type of each trial delivered
% during the session, where T is the number of trials delivered during the
% session.


% INPUTS: 
% 1) input - a path to a .txt file containing serial output from an Arduino
% running an ArduFSM protocol for a single behavior session. In accordance
% with the general ArduFSM framework, each line of output either
% acknowledges the receipt of instructions from the host PC, asserts
% upcoming trial parameters, reports recorded behavior parameters, or
% signals the start of a trial. More information about the ArduFSM
% framework can be found at: 

% https://github.com/cxrodgers/ArduFSM


% OUTPUTS:
% 1) trialTypes - a T x 1 cell array of strings describing the trial type
% of each trial delivered over the course of a behavior session, where T is
% the number of trials. Trial types are listed in the order they are
% delivered.


% TODO: 
% This function is currently hard-coded to look for trial parameters and
% return trial types specific to the multiSens protocol. It may be
% worthwhile building more flexibility into this function. s

%%
function trialTypes = read_ardulines(input)
    fileID = fopen(input);
    C = textscan(fileID, '%s', 'Delimiter', '\r\n');
    C = C{1,1};
    fclose(fileID);

    last = @(x) x(end);

    trialStartRegex = 'TRL_[1-9]*_START';
    trialStarting = ~cellfun(@isempty, regexp(C,trialStartRegex));
    trialStartLines = find(trialStarting);

    stprExp = 'TRLP STPRIDX';
    isStprLine = ~cellfun(@isempty, regexp(C, stprExp));
    stprLines = C(find(isStprLine));
    stprVals = str2num(cellfun(last, stprLines));

    spkrExp = 'TRLP SPKRIDX';
    isSpkrLine = ~cellfun(@isempty, regexp(C, spkrExp));
    spkrLines = C(find(isSpkrLine));
    spkrVals = str2num(cellfun(last, spkrLines));

    trialTypes = cell(min([length(stprVals), length(spkrVals)]), 1);

    for i = 1:length(trialTypes)
        if stprVals(i) == 0 && spkrVals(i) == 0
            trialTypes{i} = 'none';
        elseif stprVals(i) == 1 && spkrVals(i) == 0
            trialTypes{i} = 'stepper only';
        elseif stprVals(i) == 0 && spkrVals(i) == 1
            trialTypes{i} = 'speaker only';
        elseif stprVals(i) == 1 && spkrVals(i) == 1
            trialTypes{i} = 'stepper and speaker';
        end
    end
end