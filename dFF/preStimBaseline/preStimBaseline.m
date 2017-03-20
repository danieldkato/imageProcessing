
function Conditions = preStimBaseline(rawPath, trialPath, preStim, postStim, outputDir, condSettingsPath)

% Function to compute dF/F with baseline F computed as mean pre-stimulus
% period activity. This requires parsing full grabs into trials and
% organizing them by condition. 

    % Load raw image data: 
    Raw = csvread(rawPath);

    % Load trial data from .csv:
    [num, txt, Trials] = xlsread(trialPath);
    Trials = Trials(7:end, :);
    
    % Parse raw F traces into individual trials, and sort by condition:
    [Conditions, duration] = trialsByCondition(Raw, Trials, preStim, postStim, condSettingsPath);
    numConditions = length(Conditions);
    
    for c = 1:numConditions
        disp('Condition data:');
        disp(size(Conditions(c).Data));
    end
    
    % Compute dFF for every trial in every condition:
    for c = 1:numConditions
        
        % Get the number of trials in the current condition:
        numTrials = size(Conditions(c).Data, 3);
        
        % For each trial in the current condtion, compute an n x p dF/F matrix, where n is the number of ROIs and p is the number of frames in the peri-stimulus period: 
        for t = 1:numTrials
            disp('size F');
            disp(size(Conditions(c).Data(:,1,t)));
            disp('size F0');
            disp(size(mean(Conditions(c).Data(:,1:preStim,t),2)));
            Conditions(c).Data(:,:,t) = cell2mat(arrayfun(@(a) (Conditions(c).Data(:,a,t)-mean(Conditions(c).Data(:,1:preStim,t),2)) ./ mean(Conditions(c).Data(:,1:preStim,t),2), (1:1:preStim+postStim+1), 'UniformOutput', 0));
        end
        
    end
    
    %% Save output data: 
    
    %Save data into an HDF5 with one dataset per condition:
    status = exist(outputDir, 'dir');
    if status == 0
        mkdir(outputDir);
    end 

    old = cd(outputDir);
    
    for c = 1:numConditions
        dSetName = strcat(['/', Conditions(c).Name]);
        h5create('dFF_parsed_trials.h5', dSetName, size(Conditions(c).Data));
        h5write('dFF_parsed_trials.h5', dSetName, Conditions(c).Data);
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Color', Conditions(c).Color);
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Abbreviation', Conditions(c).Abbreviation);
    end
    
    h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_pre_stim', preStim);
    h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_post_stim', postStim);
    h5writeatt('dFF_parsed_trials.h5', '/', 'stim_duration', duration);
    
    % For easy manual inspection, save output data as a .xsl file, with one sheet per trial:
    for c = 1:numConditions
        numTrials = size(Conditions(c).Data,3);
        
        for t = 1:numTrials
            xlswrite('dFF_parsed_trials.xls', Conditions(c).Data(:,:,t), strcat([Conditions(c).Name, ', trial ', num2str(t)]));
        end 
    end
 
    

end