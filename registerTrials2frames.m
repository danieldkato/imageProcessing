% Last updated DDK 2016-09-27

% OVERVIEW:
% This function returns a T x 2 matrix containing the start time and trial
% type for every trial delivered during a given grab, where T is the number
% of trials delivered during the grab.


% REQUIREMENTS:
% This function requires the following data:
% 1) A trace of the analog voltage signal used to drive the slow
% scan-mirror galvanomter during the grab under analysis
% 2) A trace of an analog trial timer signal recorded during the grab 
% 3) A text file containing Arduino feedback received over a serial port
% during the grab. The script currently assumes that all analog data has
% been saved in a LabView .dat format.

% This function requires the following software:
% 1) The MATLAB function readContinuousDAT, available at https://github.com/gpierce5/BehaviorAnalysis/blob/master/readContinuousDAT.m (commit 71b3a3c)
% 2) The MATLAB function LocalMinima, available at \\hsbruno05\Users\dan\Documents\MATLAB\clay\LocalMinima.m
% 3) The MATLAB function read_ardulines, available at https://github.com/danieldkato/trial_registration/blob/master/read_ardulines.m


% INPUTS:
% 1) galvoPath - path to a trace of the analog voltage signal used to drive
% the slow scan-mirror galvanometer during the grab. Should be saved as a
% LabView .dat file. Contains information about frame start times.

% 2) timerPath - path to a trace of the analog trial timer signal recorded
% during the grab. In the current protocol, trial onset is immediately
% preceded by a 10 ms, 5 V TTL pulse sent from the Arduino responsible for
% controlling stimulus hardware. Should be saved as a LabView .dat file.
% Contains information about trial start times.

% 3) ardulines - path to a .txt file containing serial output received from
% an Arduino over the course of the grab. Contains information about trial type.

% 4) outputPath - path to the directory where the output matrix should be
% saved.

% 5) grabPath - optional string argument containing path to the raw TIFF of
% the grab being analyzed.

% 6) showInflectionPoints - optional boolean argument controlling whether
% or not to plot the galvo and timer traces along with identified frame and
% trial start times. 


% OUTPUTS: 
% This function returns a T x 2 matrix containing the start time and trial
% type of every trial delivered during the grab, where T is the number of
% complete trials delivered during the grab. Each row represents a trial;
% the first column gives the frame number during which the trial starts,
% and the second column gives the trial type. 

% This function also saves the output matrix as a .csv for future
% processing, along with metadata about the registration itself.


%%
function [trialMatrix] = registerTrials2frames(galvoPath, timerPath, ardulines, grabPath, conditionSettings, outputPath, showInflectionPoints)
    %% Read image grab and galvo metadata files to get the necessary data acquisition parameters:
    
    % Read image grab metadata to get framerate:
    [grabDir, nm, ext] = fileparts(grabPath);
    list = dir(grabDir);
    grabMeta = list(arrayfun(@(a) strcmp(a.name, 'meta.txt'), list)); % look for a file called 'meta.txt' in the same directory as the raw grab file
    
    % Raise an error if meta.txt is not found:
    if length(grabMeta) == 0
        error('Metadata file for grab not found; make sure that meta.txt is located in the same directory as raw multi-page TIFF.');
    else
        disp('grabMeta(1).name');
        disp(grabMeta(1).name);
        grabMetaFid = fopen(fullfile(grabDir, grabMeta(1).name), 'r');
        disp('grabMetaFid');
        disp(grabMetaFid);
        content = fscanf(grabMetaFid, '%c');
        eval(content);
        
        % Raise an error if meta.txt does not contain a variable called frame_rate:
        if exist('frame_rate', 'var') == 0
            error('Frame rate not found; make sure that grab metadata file includes line ''frame_rate = <f>'', where <f> stands for frame rate in Hz.');
        end
        
    end
    
    % Read galvo header to get sample rate:
    galvoFid = fopen(galvoPath, 'r', 'b');
    [headerSize, header] = SkipHeader(galvoFid);
    sampleRate = str2double(header{7}(18:end));
    
    
    %% Set output display parameters:
    if nargin< 6
        showInflectionPoints = 0;
    end
    
    if nargin < 5
        outputPath = cd;
    end
    
    
    %% Load galvo, timer and Arduino data:
    galvoTrace = readContinuousDAT(galvoPath); % Load the galvanometer data from the raw .dat file into an s x 1 vector, where s is number of samples taken during grab 
    timerTrace = readContinuousDAT(timerPath); % Load the trial timer data from the raw .dat file into an s x 1 vector, where s is the number of samples taken during a grab
    trialTypes = read_ardulines(ardulines, conditionSettings); %% Get an ordered list of trial types from arudlines
    
    
    %% Get the start time of every frame, in terms an index into the galvanometer trace:
    
    % Compute some input parameters for LocalMinima, called below:
    framePeriod = 1/frame_rate;
    minDistanceGalvo = framePeriod * sampleRate; % The function LocalMinima will include only the lowest of any local minima found within this many samples of each other
    minDistanceGalvo = minDistanceGalvo * .9; % Fudge factor; the true number of samples between certain pairs of frame-start times is slightly less than the theoretical value
    galvoThreshold = -1.6; % Whatever units gavloTrace is expressed in (Volts, I think); the function LocalMinima will exclude any local minima higher than this value; for the time being, I just got this from eyeballing a sample galvo trace, but I may ultimately need more sophisticated ways of getting this if there's any variability
    
    % Get a vector of every galvo signal sample at which a frame begins:
    frameOnsetSamples = LocalMinima(galvoTrace, minDistanceGalvo, galvoThreshold);
    
    
    %% Get the start time of every trial, in terms of an index into the trial timer signal trace:
    
    % Get every sample index in timerTrace corresponding to the onset of a
    % new trial; trial onsets are indicated by local maxima, so run
    % LocalMinima on -trialTrace:
    
    % Compute some parameters for LocalMinima, called below:
    minITI = 3; % Seconds; again, it would be better if there were a way to do this dynamically
    minDistanceTimer = minITI * sampleRate;
    timerThreshold = -4; % timerTrace units (Volts, I think); I just eyeballed this for now, but I should probably find a way to get this dynamically.
    
    % Get a vector of every timer signal sample at which a trial begins: 
    trialOnsetSamples = LocalMinima(-timerTrace, minDistanceTimer, timerThreshold);
    
    
    %% Show traces if requested by user:
    
    % Plot the local minima on top the galvo trace if desired; this can be
    % handy just to double check that reasonable parameters for LocalMinima
    % have been chosen, but may be cumbersome if processing large batches
    % of data.
    
    if showInflectionPoints == 1
        figure;
        hold on;
        t = (1:1:length(galvoTrace));
        plot(galvoTrace);
        plot(t(frameOnsetSamples), galvoTrace(frameOnsetSamples), 'r.');
        plot(timerTrace);
        plot(t(trialOnsetSamples), timerTrace(trialOnsetSamples), 'r.'); 
    end
    
    
    %% Omit any trials delivered before the first frame or after the last frame 
    trialOnsetSamples = trialOnsetSamples( trialOnsetSamples>=min(frameOnsetSamples) & trialOnsetSamples<=max(frameOnsetSamples) );
    
    
    %% Match every trial to the frame within which it started
    trialStartFrames = cell(length(trialOnsetSamples), 1);
    
    % For each sample number corresponding to a trial start time, find the highest sample number corresponding to a frame start time below it:
    for i = 1:length(trialStartFrames)
        [M, I] = max(frameOnsetSamples( frameOnsetSamples <= trialOnsetSamples(i) ));
        trialStartFrames{i} = I + 1; % have to add 1 because there's one frame that completes before the first local minimum
    end

    
    %% Merge the trial start time and trial type information:
    trialMatrix = cell(length(trialOnsetSamples), 3);
    trialMatrix(:, 1) = trialStartFrames;
    disp(size(trialOnsetSamples));
    disp(size(trialTypes));
    trialMatrix(:, 2:3) = trialTypes; 
    
    
    %% Write trialMatrix to a .csv: 
    
    % Check that the output path exists:
    status = exist(outputPath);
    if status == 0
        mkdir(outputPath);
    end
    
    % Create and open the output file for writing:
    filename = fullfile(outputPath, 'trialMatrix.csv');
    fileID = fopen(filename, 'w');
    
    % Write header:
    if exist('grabPath', 'var')
        fprintf(fileID, strcat(['Grab,', strrep(grabPath,'\','\\'), '\n']));
    else
        fprintf(fileID, strcat(['Grab, none specified \n']));
    end
    fprintf(fileID, strcat(['Galvanometer trace, ', strrep(galvoPath,'\','\\'), '\n']));
    fprintf(fileID, strcat(['Trial timer signal trace, ', strrep(timerPath,'\','\\'), '\n']));
    fprintf(fileID, strcat(['Arduino output, ', strrep(ardulines,'\','\\'), '\n']));
    fprintf(fileID, '\n');
    fprintf(fileID, 'Trial start frame number, Trial type, Trial duration (ms) \n');
    
    % Write body:
    formatSpec = '%d, %s, %d \n';
    for i = 1:size(trialMatrix, 1)
        fprintf(fileID, formatSpec, trialMatrix{i,:});
    end
    
    fclose(fileID);
    
    
    %% Write metadata:
    
    inputs = {{'galvanometer trace', galvoPath};
              {'timer signal trace', timerPath};
              {'arduino feedback', ardulines}
        };
    
    outputs = {{'trial matrix', strcat([outputPath, 'triaMatrix.csv'])}};
    parameters = {};
    
    old = cd(outputPath);
    writeMetadata('trial_registration', 'sampleDomain', inputs, outputs, parameters);
    cd(old);
    
end