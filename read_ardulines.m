    
function trialTypes = read_ardulines(input)
    %inputPath = '\\10.112.43.46\Public\dank\multiSens\raw\1146-1\2P\160830\site6\grab001\ardulines.20160830130611';
    %inputPath = '\\10.112.43.46\Public\dank\multiSens\raw\1146-1\2P\160901\site6\grab001\ardulines - Copy.20160901134327';
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