% Last updated DDK 2016-9-29

% OVERVIEW:
% This function accepts a number of cell arrays containing analysis-related
% metadata and writes them to a .txt file formatted as a dictionary-like
% structure (ideally saved alongside the outputs of the analysis itself).
% It should be called from every MATLAB function responsible for saving any
% kind of output.

% The idea is that all saved, processed data should be associated with some
% metadata file that allows the user to reconstruct exactly how that output
% was derived, all the way back to the raw data. This is described in more
% detail at:

% \\10.112.43.45\Public\dank\multiSens\analysis\README.txt



% INPUTS:
% 1) step - a string that succinctly describes what step in the data
% processing workflow is performed by the calling, output-saving function.
% Examples include 'motion_correction', 'image_segmentation',
% 'compute_dFF'.

% 2) pipeline - a string that succinctly describes what method was used to
% perform the current data processing step. For example, if the calling
% function were a way of doing image segmentation, possible values for
% pipeline might include 'NMF', or 'STICA'. This level of abstraction is
% described more at the address listed above.

% 3) inputCells - a cell array containing paths to the data processed by the
% calling function. inputs should be formatted as a cell array of 1 x 2
% cell arrays of field-value pairs, as follows:

% inputCells = { {'input field 1', 'path\to\input1'};
%            {'input field 2', 'path\to\input2'}; 
%            {'input field 3', 'path\to\input3'}};

% The second element of each nested pair should be a path to something
% processed by the calling function - e.g., a raw TIFF to be
% motion-corrected, or raw gavlo and timer data to be registered to each
% other.

% 4) outputCells - a cell array containing paths to whatever files or
% directories are created and saved by the calling function, e.g., a
% motion-corrected TIFF. This should have the same format as inputs.

% 5) paramCells - a cell array containing important parameters used to
% perform the current processing step. Should be formatted the same as
% inputs and outputs. For example:

% paramCells = { {'K', 100};
%                {'tau', 2};
%                {'p', 1}};


% OUTPUTS:
% This function does not have any formal return, but creates a file called
% meta.txt containing information about what the calling function just did.
% Note that meta.txt will just be saved in the current working directory at
% the time writeMetadata is called, so it its the calling function's
% responsibility to ensure that this function is called from meta.txt's
% desired destination (which should ideally be the lowest directory
% containing all of the outputs of the calling function).

% writeMetadata will also find whatever other functions were invoked (even
% indirectly) during execution of the calling function, and try to retrieve
% their most recent git commits to give a clear picture of all the code
% involved in the data processing step described by the metadata. 

% The final product of this function will have the format:

% step_metadata = 
% {'pipeline':'one of perhaps several alternative methods of doing the present step',
%  'inputs':{'input1':'path\to\input1',
%            'input2':'path\to\input2'},
%  'outputs':{'output1':'path\to\output1',
%             'output2':'path\to\output2'},
%  'dependencies':['path\to\dependency1 (commit)',
%                  'path\to\dependency2 (commit)',
%                  'path\to\dependency3 (commit)],
%  'parameters':{'param1':'value1',
%                'param2':'value2',
%                'param3':'value3'}
% }


% TODO: 
% Need to think about how to copy over metadata associated with inputs.
% This is critical to the metadata forwarding scheme described in: 

% \\10.112.43.45\Public\dank\multiSens\analysis\README.txt

% Should probably include some way of dealing with situations when a
% meta.txt file already exists. 

% Also, should think about whether the idea of the 'pipeline' is a useful
% level of abstraction or if it will just confuse people; on the one hand,
% in situations when there exist several alternative possible strategies
% for accomplishing the current data-processing step, I find it a succinct
% way of getting across which of several alternative strategies was used to
% accomplish the present data-processing step; on the other hand, it's not
% strictly necessary, as meta.txt already contains finer-grained
% information about what scripts, etc., were used to create the output.


%%
function writeMetadata(step, pipeline, inputCells, outputCells, paramCells)    
    % (should probably do some validation here)
    
    numLines = length(inputCells)+length(outputCells)+length(paramCells)+2;
    
    %% Reformat all cells of the form {'key', 'value'} as strings of the form " 'key':'value' " :
    
    % Make individual dictionary entries for each input, i.e., strings of
    % the form 
    
    % " 'key1':'value1' "
    
    inputDictEntries = cell(1, length(inputCells));
    for i = 1:length(inputCells)
        inputCells{i}{2} = strrep(inputCells{1}{2}, '\', '\\');
        inputDictEntries{i} = strcat(['''', inputCells{i}{1},''':''', inputCells{i}{2},'''']);
    end 
    
    % Make indiviual dictionary entries for each output:
    outputDictEntries = cell(1, length(outputCells));
    for j = 1:length(outputCells)
        outputCells{j}{2} = strrep(outputCells{j}{2}, '\', '\\');
        outputDictEntries{j} = strcat(['''', outputCells{j}{1},''':''', outputCells{j}{2},'''']);
    end 

    % Make indiviual dictionary entries for each parameter:
    paramDictEntries = cell(1, length(paramCells));
    disp('length parameters = ');
    disp(length(paramCells));
    for k = 1:length(paramCells)
        if isa(paramCells{k}{2}, 'char')
            paramDictEntries{k} = strcat(['''', paramCells{k}{1},''':''', paramCells{k}{2},'''']);
        elseif isnumeric(paramCells{k}{2})
            paramDictEntries{k} = strcat(['''', paramCells{k}{1},''':', num2str(paramCells{k}{2})]);
        end
    end
    disp('length paramEntries = ');
    disp(length(paramDictEntries));
    
    
    
    %% Get the dependencies of the calling function, and, where possible, their versions:
    
    ST = dbstack(1, '-completenames');
    [fList, pList] = matlab.codetools.requiredFilesAndProducts(ST.file);
    fListMax = 50;
    fList = fList(1:min([length(fList), fListMax])); % we can truncate this list; in practice, it tends to turn out to be several hundred
    
    calledFunctions = cell(1, length(fListMax));
    for i = 1:length(fList)
        [fullPath, commit] = getVersion(fList{i});
        fullPath = strrep(fullPath, '\', '\\');
        calledFunctions{i} = strcat([fullPath, ' ' commit]);
    end
    
    %% Write all " 'key':'value' " pairs into a dictionary:
    
    % Create the metaData cell array; each entry in this array will be a
    % line of the metadata file
    line = 1;
    metaData = cell(numLines, 1);
    metaData{line} = strcat([step, '_metadata = \r\n']); line = line + 1;
    
    metaData{line} = strcat(['{''pipeline'':''', pipeline, ''', \r\n']); line = line + 1;
    
    % Write the input dictionary to metaData
    metaData{line} = strcat([' ''inputs'':{', inputDictEntries{1}]);  %opening bracket of the input dictionary
    if length(inputDictEntries) == 1
        metaData{line} = strcat([metaData{line}, '} \r\n']); line = line + 1;
    else
        metaData{line} = strcat([metaData{line}, ', \r\n']);
        line = line + 1;
        for m = 2:length(inputDictEntries)-1
            metaData{line} = strcat(['           ', inputDictEntries{m}, ', \r\n']); line = line + 1;
        end
        metaData{line} = strcat(['           ', inputDictEntries{length(inputDictEntries)}, '}, \r\n']); line = line + 1; %closing bracket of the input dictionary
    end
    
    
    % Write the output dictionary to metaData
    metaData{line} = strcat([' ''outputs'':{', outputDictEntries{1}]); %opening bracket of the output dictionary
    if length(outputDictEntries) < 2
        metaData{line} = strcat([metaData{line}, '} \r\n']); line = line + 1;
    else
        metaData{line} = strcat([metaData{line}, ', \r\n']);
        line = line + 1;
        for n = 2:length(outputDictEntries)-1
            metaData{line} = strcat(['            ', outputDictEntries{n}, ', \r\n']); line = line + 1;
        end    
        metaData{line} = strcat(['            ', outputDictEntries{length(outputDictEntries)}, '}, \r\n']); line = line + 1; %closing bracket of output dictionary
    end
    
    %Write the calling function's dependencies to metaData
    metaData{line} = strcat([' ''dependencies'':[''', calledFunctions{1}, '''']); 
    if length(calledFunctions) < 2
        metaData{line} = strcat([metaData{line}, '} \r\n']); line = line + 1;
    else
        metaData{line} = strcat([metaData{line}, ', \r\n']);
        line = line + 1;
        for q = 2:length(calledFunctions)-1
            metaData{line} = strcat(['                 ''', calledFunctions{q}, ''', \r\n']); line = line + 1;
        end
        metaData{line} = strcat(['                 ''', calledFunctions{length(calledFunctions)}, '''] \r\n']); line = line + 1;
    end
    
    % Write the parameters dictionary to metaData
    if length(paramDictEntries) > 0
        metaData{line} = strcat([' ''parameters'':{', paramDictEntries{1}]); %opening bracket of parameters dictionary
        if length(paramDictEntries) < 2
            metaData{line} = strcat([metaData{line}, '} \r\n']); line = line + 1;
        else
            metaData{line} = strcat([metaData{line}, ', \r\n']);
            line = line + 1;
            for n = 2:length(paramDictEntries)-1
                metaData{line} = strcat(['               ', paramDictEntries{n}, ', \r\n']); line = line + 1;
            end
            metaData{line} = strcat(['               ', paramDictEntries{length(paramDictEntries)}, '}, \r\n']); line = line + 1; %closing bracket of parameters dictionary
        end
    elseif length(paramDictEntries) == 0
        metaData{line} = ' ''parameters'':{}, \r\n'; line = line + 1;%opening bracket of parameters dictionary
    end
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