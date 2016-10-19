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

function [trialMatrix] = registerTrials2frames(galvoPath, timerPath, ardulines, outputPath, grabPath, showInflectionPoints)
    
    if nargin<6
        showInflectionPoints = 0;
    end
    
    
    %% Load data:
    galvoTrace = readContinuousDAT(galvoPath); % Load the galvanometer data from the raw .dat file into an s x 1 vector, where s is number of samples taken during grab 
    timerTrace = readContinuousDAT(timerPath); % Load the trial timer data from the raw .dat file into an s x 1 vector, where s is the number of samples taken during a grab
    trialTypes = read_ardulines(ardulines); %% Get an ordered list of trial types from arudlines
    
    
    %% Get the start time of every frame (in terms an index into the galvanometer trace)
    
    % Frame start times correspond to local minima in the sawtooth pattern
    % of the galvanometer trace; these can be found using the LocalMinima
    % function.
    
    frameRate = 3.37; % Frames per second; this is a constant for now, but I should think about how to make this get the frame rate at runtime
    framePeriod = 1/frameRate;
    sampleRate = 16000; % Samples per second; Also a constant for now, but I should think about how to make this get the sample rate at runtime
    minDistanceGalvo = framePeriod * sampleRate; % The function LocalMinima will include only the lowest of any local minima found within this many samples of each other
    minDistanceGalvo = minDistanceGalvo * .9;
    galvoThreshold = -1.6; % Whatever units gavloTrace is expressed in (Volts, I think); the function LocalMinima will exclude any local minima higher than this value; for the time being, I just got this from eyeballing a sample galvo trace, but I may ultimately need more sophisticated ways of getting this if there's any variability
    
    frameOnsetSamples = LocalMinima(galvoTrace, minDistanceGalvo, galvoThreshold); %returns a vector of indices into the input trace
    
    % Plot the local minima on top the galvo trace if desired; this can be
    % handy just to double check that reasonable parameters for LocalMinima
    % have been chosen, but may be cumbersome if processing large batches
    % of data
    
    if showInflectionPoints == 1
        figure;
        plot(galvoTrace);
        hold on;
        t = (1:1:length(galvoTrace));
        plot(t(frameOnsetSamples), galvoTrace(frameOnsetSamples), 'r.');
    end 
    
    
    %% Get the start time of every trial (in terms of an index into the trial timer signal trace)
    
    % Get every sample index in timerTrace corresponding to the onset of a
    % new trial; trial onsets are indicated by local maxima, so run
    % LocalMinima on -trialTrace:
    
    minITI = 3; % Seconds; again, it would be better if there were a way to do this dynamically
    minDistanceTimer = minITI * sampleRate;
    timerThreshold = -4; % timerTrace units (Volts, I think); again, should think of a way to get this dynamically
    trialOnsetSamples = LocalMinima(-timerTrace, minDistanceTimer, timerThreshold);
    
    if showInflectionPoints == 1
        plot(timerTrace);
        plot(t(trialOnsetSamples), timerTrace(trialOnsetSamples), 'r.'); 
    end
    
    
    %% Omit any trials delivered before the first frame or after the last frame 
    trialOnsetSamples = trialOnsetSamples( trialOnsetSamples>=min(frameOnsetSamples) & trialOnsetSamples<=max(frameOnsetSamples) );
    
    
    %% Match every trial to the frame within which it started
    trialStartFrames = cell(length(trialOnsetSamples), 1);
    
    % For each sample number corresponding to a trial start time, find the
    % highest sample number corresponding to a frame start time below it
    for i = 1:length(trialStartFrames)
        [M, I] = max(frameOnsetSamples( frameOnsetSamples <= trialOnsetSamples(i) ));
        trialStartFrames{i} = I;
    end

    
    %% Merge the trial start time and trial type information
    trialMatrix = cell(length(trialOnsetSamples), 2);
    trialMatrix(:, 1) = trialStartFrames;
    disp(size(trialOnsetSamples));
    disp(size(trialTypes));
    trialMatrix(:, 2) = trialTypes; 
    
    
    
    %% Write trialMatrix to a .csv 
    filename = fullfile(outputPath, 'trialMatrix.csv');
    fileID = fopen(filename, 'w');
    %fileID = fopen('trialMatrix.csv', 'w');
    
    % write header:
    if exist('grabPath', 'var')
        fprintf(fileID, strcat(['Grab,', strrep(grabPath,'\','\\'), '\n']));
    else
        fprintf(fileID, strcat(['Grab, none specified \n']));
    end
    fprintf(fileID, strcat(['Galvanometer trace, ', strrep(galvoPath,'\','\\'), '\n']));
    fprintf(fileID, strcat(['Trial timer signal trace, ', strrep(timerPath,'\','\\'), '\n']));
    fprintf(fileID, strcat(['Arduino output, ', strrep(ardulines,'\','\\'), '\n']));
    fprintf(fileID, '\n');
    fprintf(fileID, 'Trial start frame number, Trial type \n');
    
    % write trialMatrix:
    formatSpec = '%d, %s \n';
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