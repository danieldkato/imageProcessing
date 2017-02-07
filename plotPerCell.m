% Last updated DDK 2017-02-07

% OVERVIEW:
% This function creates peri-stimulus activity plots for every unit recorded
% from a given acquisition session. For every unit, this function creates 2
% plots:

% 1) A plot with mean activity traces plus SEM bars, and
% 2) A plot with activity traces for every individual trial

% All plots are saved to storage in an output directory specified in the
% function call.

% If the data include responses to multiple stimulus or trial conditions,
% then individual activity traces will be color-coded according to the
% mapping defined either in the data file itself or in a condition settings
% file specified in the function call (see below for detail).

% Note that this function itself is ONLY responsible for taking the mean
% and SEM of the data passed as input and plotting it. It dos NOT perform
% any other computation, e.g. compute dF/F from raw fluorescence values.
% Thus, if the data passed as input are raw fluorescence values, then
% the values plotted will be raw fluorescence values. 


% REQUIREMENTS:
% This function should only be used on datasets where the stimulus period
% for all trials is the same. This function draws a shaded rectangle over
% the area corresponding to the stimulus period on every plot, so the
% figures are really only sensible if the stimulus period for all trials is
% the same.


% INPUTS:
% 1) activity - path to an HDF5 file containing the data to be plotted.
% Data must be parsed as follows: the HDF5 must include one dataset for
% each trial or stimulus condition to be analyzed. Each dataset should be
% an N x T x P activity matrix, where N is the number of ROIs to be
% analyzed, T is the number of frames in the peri-stimulus period to be
% plotted, and P is the number of repetitions of the corresponding trial or
% stimulus condition. N and and T must be the same for all datasets, but P
% may be different for each.

% In order to properly lay out figure objects like shaded rectangles for
% stimulus epochs, the HDF5 root must have the following attributes:
    
%   a) num_samples_pre_stim: the number of frames before stimulus onset
%   included in the trace for each trial

%   b) num_samples_post_stim: the number of frames after stimulus onset
%   included in the trace for each trial

% In addition to the required attributes described above, the HDF5 may also
% have the following optional attributes for specifying the appearance of
% figures:


% 2) outputDirectory - directory where all created figures should be saved.

% 3) grabMetadata - path to a MATLAB-evaluable .txt file containing
% information about the acquisition session. This must include the
% data acquisition rate in samples per second. 

% 4) conditionSettings - path to a MATLAB-evaluable .txt file defining a 
% C x 1 cell array of structs, where C is the number of distinct trial
% conditions presented during throughout the course of the trials to be
% analyzed. Each struct must have at least the following three fields:

%   a) Name - the name of the trial condition. This must exactly match the
%   trial condition descriptions in the second column of trials, described
%   above.

%   b) Color - color code for the given condition, in any valid MATLAB
%   format for encoding color. Will be used in plotting. 

%   c) abbreviation - abbreviation for the given trial condition. Will be
%   used in creating legends for each figure.


% OUTPUTS:
% This function saves to disk 2 plots for each ROI from the analyzed
% dataset:

% 1) A plot with mean peri-stimulus activity traces plus SEM bars for each condition, and
% 2) A plot with peri-stimulus traces for every individual trial, color-coded by condition

% In addition, this function formally returns:
% 1) meanPaths - an N x 1 cell array of strings containing full paths to
% the saved mean activity figures, and 
% 2) rawPaths - an N x 1 cell array of strings containing full paths to the
% saved individual trial activity figures. 


% TODO:
% 1) Stimulus duration is currently hard-coded into the script. This should
% also be read in dynamically from somewhere else, like trials (although
% trials currently only includes trial start frame and condition, not
% duration, so this would entail changes to how the trial matrix is
% created). In getting the trial duration from trials, one could also
% verify that all trial durations are the same (which is the only situation
% in which the shaded stimulus period rectangle makes sense).

% 2) This function requires that the condition names in trials exactly
% match the condition names in conditions. Should probably throw up an
% error if this doesn't hold. 

% 3) Maybe think about making the conditions input optional, as it's
% probably really not necessary; it could just read out the trial
% conditions from trials, then automatically assign colors and
% abbreviations to each if none are specified in conditions.

% 4) Perhaps make the outputDirectory argument optional, and have it
% default to the current working directory.


%%
function [meanPaths, rawPaths] = plotPerCell(activity, outputDirectory, grabMetadata, conditionSettings)
    %% Load data, spell out some basic parameters:
       
    stimDur = 2; % this needs to be found dynamically
    
    %% Load grab metadata if available:
    fid2 = fopen(grabMetadata);
    grabMeta = fscanf(fid2, '%c');
    eval(grabMeta);
    
    
    %% Load condition metadata:
    info = h5info(activity);
    numConditions = length(info.Datasets);
    
    if nargin < 4 % if no conditionSettings.txt is specified, get condition metadata from HDF5 attributes:
        for c = 1:numConditions
            Conditions(c).Name = info.Datasets(c).Name;
            Conditions(c).Abbreviation = h5readatt(activity, strcat(['/', info.Datasets(c).Name]), 'Abbreviation');
            Conditions(c).Color = h5readatt(activity, strcat('/', info.Datasets(c).Color), 'Color');
        end
    else % if conditionSettings.txt is specified, get condition metadata from there:
       fid = fopen(conditionSettings);
       fidContent = fscanf(fid, '%c');
       eval(fidContent);
       
       % (should include some code here to check that eval(fidContent) actually creates a struct called Conditions)
       % (also might want to check whether it has appropriate attributes)
       
    end
    
    % (should auto-assign abbreviation and color if condition settings are specified in neither conditionSettings.txt or HDF5 attributes)
    
    
    %% Load activity data:
    
    % Recall that if conditionSettings.txt is specified, steps must be
    % taken to ensure that data and condition metadata match up; e.g., the
    % datasets in the HDF5 might be ordered 'stepper only', 'speaker only',
    % then 'stepper and speaker', while in conditionSettings.txt they may
    % be ordered 'speaker only', 'stepper and speaker', then 'stepper
    % only'. Thus, conditions need to be explicitly matched with their
    % metadata by name, and not just by index.
    conditionNames = extractfield(Conditions, 'Name');
    for c = 1:numConditions
        cIdx = find(arrayfun(@(a) strcmp(a, info.Datasets(c).Name), conditionNames));
        Conditions(cIdx).Data = h5read(activity, strcat(['/', info.Datasets(c).Name]));
    end
    
    % Convert peri-stimulus period from seconds to samples:
    preStimSamples = h5readatt(activity, '/', 'num_samples_pre_stim');
    postStimSamples = h5readatt(activity, '/', 'num_samples_post_stim');
    periStimPeriod = preStimSamples + postStimSamples + 1;
    
    % (should set default behavior for if these attributes aren't found)
    
    % For convenience and human readability:
    numROIs = size(Conditions(1).Data,1); % should probably have something to confirm that everything has the same number of ROIs
       
    
    %% Compute means and SEM for each condition for each ROI:
    
    for c = 1:numConditions
        Conditions(c).Means = mean(Conditions(c).Data,3);
        Conditions(c).SEM = std(Conditions(c).Data, 0, 3)/sqrt(size(Conditions(c).Data, 3));
    end
    
    
    %% Compute some useful variables for plotting:
    
    % Compute appropriate x-tick coordinates (in data units) so that there is one tick at stimulus onset and ticks marks every 2 sec:
    secPerTick = 2;     samplesPerTick = frame_rate*secPerTick;
    nTicksBelowTrialStart = floor((preStimSamples)/samplesPerTick);     ticksBelowTrialStart = (preStimSamples+1) - fliplr(samplesPerTick*(0:1:nTicksBelowTrialStart));
    nTicksAboveTrialStart = floor((postStimSamples)/samplesPerTick);    ticksAboveTrialStart = (preStimSamples+1) + samplesPerTick*(1:1:nTicksAboveTrialStart);
    xticks = [ticksBelowTrialStart ticksAboveTrialStart];
    
    % Write x-tick labels into a cell array:
    labels = (xticks - (preStimSamples+1))/frame_rate;
    labels = arrayfun(@(a) num2str(a), labels, 'UniformOutput', 0);
    
    % Define some vectors that will be used for plotting shaded SEM areas:
    domain = (1:1:preStimSamples+postStimSamples+1);
    semX = [domain fliplr(domain)]; 
    
    % Define vectors that will be used for plotting shaded stim period rectangle (recall, stim begins at index preStim + 1):
    stimPeriod = [preStimSamples+1 preStimSamples+1+stimDur*frame_rate];
    
    % Create a cell array that stores how many trials of each condition were delivered:
    numTrialsPerC = arrayfun(@(c) size(c.Data,3), Conditions, 'UniformOutput', 0)';
    
    % Write legend text into a cell array:
    legText = arrayfun(@(c) strcat([c.Abbreviation, ', n=', num2str(size(c.Data,3))]), Conditions, 'UniformOutput', 0);
    
    
    %% Prepare figure windows and file paths for plotting: 
    
    % Create cell arrays that will contain full paths to created figures; these will be returned to calling function:
    meanPaths = cell(numROIs, 1);
    rawPaths = cell(numROIs, 1);
    
    % Create figure windows (these will be erased then reused between ROIs):
    meanFig = figure(); % one for mean traces
    rawFig = figure(); % one for raw traces
    figures = {meanFig, rawFig};
    titles = {'mean', 'individual'};
    outputs = {meanPaths, rawPaths};
        
    % Make sure the output directory exists, create it if it doesn't then cd into it:
    status = exist(outputDirectory, 'dir');
    if status == 0
        mkdir(outputDirectory);
    end
    old = cd(outputDirectory);
    
    
    %% For each ROI, create 2 figures: 
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
            meanPlotHandles(c) = plot(Conditions(c).Means(r,:), 'Color', Conditions(c).Color, 'LineWidth', 1.25);
            
            % Plot SEM bars of current condition for current ROI:
            semOver = [Conditions(c).Means(r,:) fliplr( Conditions(c).Means(r,:) + Conditions(c).SEM(r,:) )];
            disp('size(semOver)');
            disp(size(semOver));
            disp('size(semX)');
            disp(size(semX));
            semUnder = [Conditions(c).Means(r,:) fliplr( Conditions(c).Means(r,:) - Conditions(c).SEM(r,:) )];
            semOverPatch = patch(semX, semOver, Conditions(c).Color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
            semUnderPatch = patch(semX, semUnder, Conditions(c).Color, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
            
            % Plot traces for every individual trial of current condition for current ROI:
            figure(rawFig);
            hold on;
            currROIallTraces = reshape(Conditions(c).Data(r,:,:), size(Conditions(c).Data,3), periStimPeriod);
            h = plot(currROIallTraces', 'Color', Conditions(c).Color);
            rawPlotHandles(c) = h(1);
        end
        
        % Create legends:
        figure(meanFig);
        legend(meanPlotHandles, legText);
        legend('boxoff');
        
        figure(rawFig);
        legend(rawPlotHandles, legText);
        legend('boxoff');
        
        % Format and save mean and individual trace figures for the current ROI:
        for f = 1:length(figures)
            
            % Draw a rectangle covering the stimulus period:
            figure(figures{f});
            hold on;
            yl = ylim;
            recY = [yl fliplr(yl)];
            %h = fill([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, [0.01 0.01 0.01], 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
            p = patch([stimPeriod(1) stimPeriod(1) stimPeriod(2) stimPeriod(2)], recY, [0.75 0.75 0.75], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
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
            
            if mod(f,2) == 0
                rawPaths{r} = fullfile(outputDirectory, titleStr);
            else
                meanPaths{r} = fullfile(outputDirectory, titleStr);
            end
            
            % Clear figure for the next ROI:
            clf(gcf)
        end
    end
    
end