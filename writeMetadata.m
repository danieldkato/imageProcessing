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
    
    
    % For each input, make individual dictionary entries, i.e., strings of
    % the form 
    
    % " 'key1':'value1' "
    
    inputEntries = cell(1, length(inputs));
    for i = 1:length(inputs)
        inputEntries{i} = strcat(['''', inputs{i}{1},''':''', inputs{i}{2},'''']);
    end 
    
    
    % Make indiviual dictionary entries for each output:
    outputEntries = cell(1, length(inputs));
    for j = 1:length(outputs)
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
    
    
    %Write the output dictionary to metaData
    metaData{line} = strcat([' ''outputs'':{', outputEntries{1}, ', \r\n']); line = line + 1; %opening bracket of the output dictionary
    if length(outputEntries) > 2
        for n = 2:length(outputEntries)-1
            metaData{line} = strcat(['            ', outputEntries{n}, ', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['            ', outputEntries{length(outputEntries)}, '}, \r\n']); line = line + 1; %closing bracket of output dictionary
    
    
    metaData{line} = strcat([' ''pipeline'':', pipeline, ', \r\n']); line = line + 1;
    
    
    %Write the parameters dictionary to metaData
    metaData{line} = strcat([' ''parameters'':{', paramEntries{1}, ', \r\n']); line = line + 1; %opening bracket of parameters dictionary
    if length(paramEntries) > 2
        for n = 2:length(paramEntries)-1
            metaData{line} = strcat(['               ', paramEntries{n}, ', \r\n']); line = line + 1;
        end
    end
    metaData{line} = strcat(['               ', paramEntries{length(paramEntries)}, '}, \r\n']); line = line + 1; %closing bracket of parameters dictionary
    metaData{line} = '} \r\n';
    
    %create the file into which metaData will be saved:
    fileID = fopen('meta.txt', 'w');
    
    for p = 1:length(metaData)
        fprintf(fileID, metaData{p});
    end
    
    fclose(fileID);
end