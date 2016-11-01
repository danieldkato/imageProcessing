
function preStimBaseline(Raw, Trials, preStim, postStim, condSettingsPath)

% Function to compute dF/F with baseline F computed as mean pre-stimulus
% period activity. This requires parsing full grabs into trials and
% organizing them by condition. 

    % Need to trim any trials that begin within preStim frames of the grab
    % start or postStim frames of grab end:
    
    
    
    % Parse raw F traces into individual trials, and sort by condition:
    [Conditions, duration] = trialsByCondition(Raw, Trials, preStim, postStim, condSettingsPath);
    
    % Compute dFF for every trial in every condition:
    for c = 1:length(Conditions)
        
        % Get the number of trials in the current condition:
        numTrials = size(Conditions(c).Data, 3);
        
        % For each trial in the current condtion, compute an n x p dF/F matrix, where n is the number of ROIs and p is the number of frames in the peri-stimulus period: 
        for t = 1:numTrials
            Conditions(c).Data(:,:,t) = arrayfun(@(a) (Conditions(c).Data(:,a,t)-mean(Conditions(c).Data(1:preStim,a,t),2)) ./ mean(Conditions(c).Data(1:preStim,a,t),2), (1:1:preStim+postStim+1), 'UniformOutput', 0);
        end
        
    end
    
    % Save data into an HDF5 with one dataset per condition:
    for c = 1:length(Conditions)
        dSetName = strcat(['/', Conditions(c).Name]);
        h5create('dFF_parsed_trials.h5', dSetName, size(Conditions(c).Data));
        h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_pre_stim', preStim);
        h5writeatt('dFF_parsed_trials.h5', '/', 'num_samples_post_stim', postStim);
        h5writeatt('dFF_parsed_trials.h5', '/', 'stim_duration', duration);
        h5write('dFF_parsed_trials.h5', dSetName, Conditions(c).Data);
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Color', Conditions(c).Color);
        h5writeatt('dFF_parsed_trials.h5', dSetName, 'Abbreviation', Conditions(c).Abbreviation);
    end
    
    

end