%vargin = {inputs, outputs, parameters}
%inputs = {{'inKey1','inVal1'},{'inKey2','inVal2'}}
%outputs = {{'outKey1','outValue2'},{'outKey2','outVal2'}}
%parameters = {{'paramKey1','paramVal1'},{'paramKey2','paramVal2'}}

%number of lines:
%metadata  =          1+
%{ input line 1
%  inputs line 2 <--  n inputs +
%  output line 1 
%  output line 2 <-- n outputs +
%  param line 1
% param line 2  <-- n params +
% close bracket      1
function writeMetadata(step, pipeline, varargin)
    
    disp(size(varargin));

    % Parse vargin into information about inputs, outputs and parameters:
    inputs = varargin{1};
    outputs = varargin{2};
    parameters = varargin{3};
    
    % (should probably do some validation here)
    
    numLines = length(inputs)+length(outputs)+length(parameters)+2;
    
    %% Reformat all cells of the form {'key', 'value'} as strings of the form " 'key':'value' " 
    
    % Make individual dictionary entries for each input, i.e., strings of
    % the form 
    
    % " 'key1':'value1' "
    
    inputEntries = cell(1, length(inputs));
    for i = 1:length(inputs)
        if isempty((strfind(inputs{i}{2}, '\')))
            splitted = strsplit(inputs{i}{2}, '\');
            inputs{i}{2} = strjoin(splitted, '\\\\');
        end
        inputEntries{i} = strcat(['''', inputs{i}{1},''':''', inputs{i}{2},'''']);
    end 
    
    
    % Make indiviual dictionary entries for each output:
    outputEntries = cell(1, length(inputs));
    for j = 1:length(outputs)
        if ~isempty(strfind(outputs{i}{2}, '\'))
            splitted = strsplit(outputs{i}{2}, '\');
            outputs{i}{2} = strjoin(splitted, '\\\\');
        end
        outputEntries{j} = strcat(['''', outputs{j}{1},''':''', outputs{j}{2},'''']);
    end 
    
    
    % Make indiviual dictionary entries for each parameter:
    paramEntries = cell(1, length(inputs));
    for k = 1:length(parameters)
        if isa(parameters{k}{2}, 'char')
            paramEntries{k} = strcat(['''', parameters{k}{1},''':''', parameters{k}{2},'''']);
        elseif isnumeric(parameters{k}{2})
            paramEntries{k} = strcat(['''', parameters{k}{1},''':', num2str(parameters{k}{2})]);
        end
    end
    
    
    
    %% Get the dependencies of the calling function, and, where possible, their versions
    
    ST = dbstack(1, '-completenames');
    [fList, pList] = matlab.codetools.requiredFilesAndProducts(ST.file);
    fListMax = 50;
    fList = fList(1:fListMax); % we can truncate this list; in practice, it tends to turn out to be several hundred
    
    calledFunctions = cell(1, length(fListMax));
    for i = 1:length(fList)
        [fullPath, commit] = getVersion(fList{i});
        if ~isempty((strfind(fullPath, '\')))
            splitted = strsplit(fullPath, '\');
            fullPath = strjoin(splitted, '\\\\');
        end
        calledFunctions{i} = strcat([fullPath, ' ' commit]);
    end

    disp(strcat(['length(calledFunctions) = ', num2str(length(calledFunctions)) ]));
    
    %% Write all " 'key':'value' " pairs into a dictionary with the format:
    
    %
    %
    %
    
    % Create the metaData cell array; each entry in this array will be a
    % line of the metadata file
    line = 1;
    metaData = cell(numLines, 1);
    metaData{line} = strcat([step, '_metadata = \r\n']); line = line + 1;
    
    
    % Write the input dictionary to metaData
    metaData{line} = strcat(['{''inputs'':{', inputEntries{1}, ', \r\n']); line = line + 1; %opening bracket of the input dictionary
    if length(inputEntries) > 2
        for m = 2:length(inputEntries)-1
            metaData{line} = strcat(['           ', inputEntries{m}, ', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['           ', inputEntries{length(inputEntries)}, '}, \r\n']); line = line + 1; %closing bracket of the input dictionary
    
    
    % Write the output dictionary to metaData
    metaData{line} = strcat([' ''outputs'':{', outputEntries{1}, ', \r\n']); line = line + 1; %opening bracket of the output dictionary
    if length(outputEntries) > 2
        for n = 2:length(outputEntries)-1
            metaData{line} = strcat(['            ', outputEntries{n}, ', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['            ', outputEntries{length(outputEntries)}, '}, \r\n']); line = line + 1; %closing bracket of output dictionary
    
    
    metaData{line} = strcat([' ''pipeline'':''', pipeline, ''', \r\n']); line = line + 1;
    
    
    % Write the parameters dictionary to metaData
    metaData{line} = strcat([' ''parameters'':{', paramEntries{1}, ', \r\n']); line = line + 1; %opening bracket of parameters dictionary
    if length(paramEntries) > 2
        for n = 2:length(paramEntries)-1
            metaData{line} = strcat(['               ', paramEntries{n}, ', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['               ', paramEntries{length(paramEntries)}, '}, \r\n']); line = line + 1; %closing bracket of parameters dictionary
    
    %Write the calling function's dependencies to metaData
    metaData{line} = strcat([' ''dependencies'':[''', calledFunctions{1}, ''', \r\n']); line = line + 1;
    if length(calledFunctions) > 2
        for q = 2:length(calledFunctions)-1
            metaData{line} = strcat(['                 ''', calledFunctions{q}, ''', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['                 ''', calledFunctions{length(calledFunctions)}, '''] \r\n']); line = line + 1;
    metaData{line} = '} \r\n';
    
    %% Save as output
    
    % Create the file into which metaData will be saved:
    fileID = fopen('meta.txt', 'w');
    
    % Write the metadata into the file:
    for p = 1:length(metaData)
        fprintf(fileID, metaData{p});
    end
    
    fclose(fileID);
end