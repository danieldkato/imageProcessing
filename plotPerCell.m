% Last updated DDK 2016-10-23

% OVERVIEW:
% This function produces summary plots for every ROI recorded from a given
% acquisition session during which multiple different trial or stimulus
% conditions are presented. For every ROI, this function produces 2 plots:

% 1) A plot with mean dF/F traces plus SEM bars for each condition, and
% 2) A plot with dF/F traces for every individual trial, color-coded by condition


% REQUIREMENTS:
% 1) The MATLAB function trialsByCondition, available at:
% https://github.com/danieldkato/trial_registration/blob/master/trialsByCondition.m

% This function should only be used on datasets where the stimulus period
% for all trials is the same. This function draws a shaded rectangle over
% the area corresponding to the stimulus period on every plot, so the
% figures are really only sensible if the stimulus period for all trials is
% the same.


% INPUTS:
% 1) activity - N x T activity matrix, where N is the number of ROIs and T
% is the number of frames to be analyzed.

% 2) trials - S x 2 cell array, where S is the number of trials to be
% analyzed. Each row corresponds to a trial; the first column of each row
% is the frame start number of a given trial, and the second column is a
% string describing the trial condition.

% 3) conditions - a C x 1 cell array of structs, where C is the number of
% distinct trial conditions presented during throughout the course of the
% trials to be analyzed. Each struct must have at least the following three
% fields:

%   a) Name - the name of the trial condition. This must exactly match the
%   trial condition descriptions in the second column of trials, described
%   above.

%   b) Color - color code for the given condition, in any valid MATLAB
%   format for encoding color. Will be used in plotting. 

%   c) abbreviation - abbreviation for the given trial condition. Will be
%   used in creating legends for each figure.

% 4) preStimTime - amount of time before stimulus onset from which to plot data, in seconds.

% 5) postStimTime - amount of time after stimulus onset from which to plot data, in seconds. 

% 6) outputDirectory - directory where all created figures should be saved.


% OUTPUTS:
% This function saves to disk 2 plots for each ROI from the analyzed
% dataset:

% 1) A plot with mean dF/F traces plus SEM bars for each condition, and
% 2) A plot with dF/F traces for every individual trial, color-coded by condition

% In addition, this function formally returns:
% 1) meanPaths - an N x 1 cell array of strings containing full paths to
% the saved mean dF/F figures, and 
% 2) rawPaths - an N x 1 cell array of strings containing full paths to the
% saved individual trial dF/F figures. 


% TODO:
% 1) The framerate is currently hard-coded into this script. Framerate
% should be read in dynamically from some parameters file.

% 2) Stimulus duration is currently hard-coded into the script. This should
% also be read in dynamically from somewhere else, like trials (although
% trials currently only includes trial start frame and condition, not
% duration, so this would entail changes to how the trial matrix is
% created). In getting the trial duration from trials, one could also
% verify that all trial durations are the same (which is the only situation
% in which the shaded stimulus period rectangle makes sense).

% 3) This function requires that the condition names in trials exactly
% match the condition names in conditions. Should probably throw up an
% error if this doesn't hold. 

% 4) Maybe think about making the conditions input optional, as it's
% probably really not necessary; it could just read out the trial
% conditions from trials, then automatically assign colors and
% abbreviations to each if none are specified in conditions.

% 5) Perhaps make the outputDirectory argument optional, and have it
% default to the current working directory.


%%
function [meanPaths, rawPaths] = plotPerCell(activity, trials, conditions, preStimTime, postStimTime, outputDirectory)
    %% Load data, spell out some basic parameters:
    
    % Imaging/trial parameters: currently hard-coded, need to get this dynamically
    frameRate = 3.37; % Hz 
    stimDur = 2;  % seconds
    
    numROIs = size(activity, 1);
    numConditions = length(conditions);
    
    preStimSamples = ceil(preStimTime * frameRate);
    postStimSamples = ceil(postStimTime * frameRate);
    periStimPeriod = preStimSamples + postStimSamples + 1;
    
    % Parse trials into TBC (for 'trials by condition'), a C x 2 matrix,
    % where C is the number of condition types. The first element of each
    % row is a condition name, and the second element of each row is an N x
    % P x D matrix, where N is the number of cells recorded in the current
    % acquisition, P is the number of samples to plot around each trial
    % onset, and D is the number of trials for the given condition.
    TBC = trialsByCondition(activity, trials, conditions, preStimSamples, postStimSamples); 
    
    % Create a cell array that stores how many trials of each condition
    % were delivered (will be necessary for plot legends):
    numTrialsPerC = cellfun(@(c) size(c,3), TBC(:,2), 'UniformOutput', 0)';
    
    
    %% Compute means and SEM for each condition for each ROI:
    
    % For each condition, get an N x P means matrix
    means = cellfun(@(c) mean(c, 3), TBC(:,2), 'UniformOutput', 0);
    
    % For each condition, get an N x P standard error of the mean matrix
    SEM = cellfun(@(c) std(c,0,3)/sqrt(size(c,3)), TBC(:,2), 'UniformOutput', 0);
    
    
    %% Create plots:
    
    % Make sure the output directory exists:
    status = exist(outputDirectory, 'dir');
    if status == 0
        mkdir(outputDirectory);
    end
    old = cd(outputDirectory);
    
    % Create cell arrays that will contain full paths to created figures;
    % will be returned to calling function
    meanPaths = cell(numROIs, 1);
    rawPaths = cell(numROIs, 1);
    
    % Define some vectors that will be used for plotting shaded SEM areas:
    domain = (1:1:preStimSamples+postStimSamples+1);
    semX = [domain fliplr(domain)]; 
    
    
    % Define vectors that will be used for plotting shaded stim period
    % rectangle (recall, stim begins at index preStim + 1):
    stimPeriod = [preStimSamples+1 preStimSamples+1+stimDur*frameRate];
    
    
    % Compute appropriate x-tick coordinates (in data units) so that there
    % is one tick at stimulus onset and ticks marks every 2 sec
    secPerTick = 2;
    samplesPerTick = frameRate*secPerTick;
    
    nTicksBelowTrialStart = floor((preStimSamples)/samplesPerTick);
    ticksBelowTrialStart = (preStimSamples+1) - fliplr(samplesPerTick*(0:1:nTicksBelowTrialStart));

    nTicksAboveTrialStart = floor((postStimSamples)/samplesPerTick);
    ticksAboveTrialStart = (preStimSamples+1) + samplesPerTick*(1:1:nTicksAboveTrialStart);

    xticks = [ticksBelowTrialStart ticksAboveTrialStart];
    
    % Create cell array with appropriate labels for x-ticks:
    labels = (xticks - (preStimSamples+1))/frameRate;
    labels = arrayfun(@(a) num2str(a), labels, 'UniformOutput', 0);
    
    
    % Create legend text:
    legText = cellfun(@(c, d) strcat([c.abbreviation, ', n=', num2str(d)]), conditions, numTrialsPerC, 'UniformOutput', 0);

    
    % Create figure windows (these will be erased then reused between ROIs):
    meanFig = figure(); % one for mean traces
    rawFig = figure(); % one for raw traces
    figures = {meanFig, rawFig};
    titles = {'mean', 'individual'};
    outputs = {meanPaths, rawPaths};
        
    % For each ROI, create 2 figures: 
    % 1) a figure plotting mean traces (with SEM bars) for each condition
    % 2) a figure plotting traces for all individual trials (color coded by condition)
    for r = 1:numROIs
        
        % Initialize empty vectors for storing numeric handles to plots; these will be necessary for creating legends
        meanPlotHandles = NaN(1,numConditions); 
        rawPlotHandles = NaN(1,numConditions); 
        
        % For current ROI, plot data from each condition:
        for c = 1:numConditions
            
            % Plot mean of current condition for current ROI:
            figure(meanFig);
            hold on;
            meanPlotHandles(c) = plot(means{c}(r,:), 'Color', conditions{c}.Color, 'LineWidth', 1.25, 'DisplayName', conditions{c}.Name);
            
            % Plot SEM bars of current condition for current ROI:
            semOver = [means{c}(r,:) fliplr( means{c}(r,:) + SEM{c}(r,:) )];
            semUnder = [means{c}(r,:) fliplr( means{c}(r,:) - SEM{c}(r,:) )];
            semOverPatch = patch(semX, semOver, conditions{c}.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
            semUnderPatch = patch(semX, semUnder, conditions{c}.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
            
            % Plot traces for every individual trial of current condition for current ROI:
            figure(rawFig);
            hold on;
            currROIallTraces = reshape(TBC{c,2}(r,:,:), numTrialsPerC{c}, periStimPeriod);
            h = plot(currROIallTraces', 'Color', conditions{c}.Color);
            rawPlotHandles(c) = h(1);
        end
        
        % Create legends:
        figure(meanFig);
        legend(meanPlotHandles, legText);
        legend('boxoff');
        
        figure(rawFig);
        legend(rawPlotHandles, legText);
        legend('boxoff');
        
        % For both figures for the current ROI, format properly and save:
        for f = 1:length(figures)
            
            % Draw a rectangle covering the stimulus period:
            figure(figures{f});
            hold on;
            yl = ylim;
            recY = [yl fliplr(yl)];
            h = fill([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, 'b', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
            uistack(h, 'bottom');
            ylim(yl);
            
            % Plot y = 0 for reference:
            plot(zeros(1,periStimPeriod), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.1);
            
            % Position and label x-ticks appropriately:
            set(gca, 'XTick', xticks, 'XTickLabel', labels);            
            
            % Label x- and y- axes:
            ylabel('dF/F (a.u.)');
            xlabel('Time rel. to stim onset (s)');
            
            % Print title:
            title(strcat(['ROI #', num2str(r), ' ', titles{f}, ' dF/F traces by condition']));
            
            xl = xlim;
            xlim([1 periStimPeriod]);
            
            % Save figure:
            if r < 10
                numStr = strcat(['00', num2str(r)]);
            elseif r >=10 && r <100
                numStr = strcat(['0', num2str(r)]);
            elseif r>100
                numStr = num2str(r);
            end
            
            titleStr = strcat(['ROI_', numStr, '_', titles{f}, '_traces']);
            saveas(gcf, titleStr);
            %disp(fullfile(outputDirectory, titleStr));
            %outputs{f}{r} = fullfile(outputDirectory, titleStr);
            
            if mod(f,2) == 0
                rawPaths{r} = fullfile(outputDirectory, titleStr);
            else
                meanPaths{r} = fullfile(outputDirectory, titleStr);
            end
            
            % Clear figure for the next ROI:
            clf(gcf)
        end
        
        %{
        % Print title, label axes:
        figure(meanFig);
        title(strcat(['ROI #', num2str(r), ' mean dF/F traces by condition']));
        
        figure(rawFig);
        title(strcat(['ROI #', num2str(r), ' individual dF/F traces by condition']));

        
        % Save figure:
        figure(meanFig);
        titleStr = strcat(['ROI', num2str(r), '_mean_activity']);
        saveas(gcf, titleStr, 'fig');
        
        figure(rawFig);
        titleStr = strcat(['ROI', num2str(r), '_individual_traces']);
        saveas(gcf, titleStr, 'fig');
        
       
        clf(meanFig);
        clf(rawFig);
        %}
    end
    
end