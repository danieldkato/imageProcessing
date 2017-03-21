% Last updated DDK 2016-10-20

% OVERVIEW:
% This function produces a summary plot for a given grab, inlcuding full
% activity traces for all ROIs and trial timing and type information. 

% The plot includes a full trace for each ROI throughout the entire course
% of the grab, with frame number on the x-axis and dF/F on the y-axis; note
% that traces are stacked on top of each other along the y-axis for
% legibility, so dF/F values are not absolute but rather are just there to
% illustrate relative change over time.

% In addition, the plot includes shaded rectangles denoting the times of
% stimulus presentations. Rectangles are color-coded according to trial
% type.


% INPUTS: 
% 1) inputTraces - path to a .csv containing a K x F dF/F matrix from the
% given grab, where K is the number of ROIs in the given grab and F is the
% number of frames.

% 2) inputTrials - path to a .csv containing a T x 2 trial matrix, where T
% is the number of trials in the given grab. The first column should of
% each row should contain the starting frame number of a given trial, and
% the second column should state the trial type. 


% OUTPUTS:
% 1) figPath - full path to the created figure.


% TODO:
% 1) It turns out that overlaying a rectangle across the vertical extent of
% the figure for each stimulus delivery significantly reduces legibility
% given how many and how closely-spaced they are. Should think about
% alternative ways of representing trial information.

% 2) Frame rate and trial duration are currently hard-coded in this
% function. Should find a way around this. 

% 3) Think about making conditions optional in case someone just wants
% traces


%%
function figPath = plotOverview(activity, trials, conditions, outputDirectory)
    %% Plot dF/F traces:
    numROIs = size(activity,1);
    numFrames = size(activity,2);
    
    % Sort activity traces by ascending maximum value:
    maxes = max(activity, [], 2);
    [maxes, ind] = sortrows(maxes); 
    activity = activity(ind, :); 
    
    % To ensure legibility, make sure that the vertical space between any
    % two traces is at least the maximum height of the lower trace:
    rowCenters = cumsum(maxes);
    dat2plot = repmat(rowCenters, 1, numFrames);
    dat2plot = dat2plot + activity;
 
    % By convention, when MATLAB's plot() function is given a matrix as an
    % input, it creates a separate plot for each column, rather than each
    % row. Since we want one plot for each ROI, we need one column per ROI.
    % The original data has one row per ROI, so transpose to get one column
    % per ROI.
    dat2plot = dat2plot';  
    
    % Draw the plots:
    figure;
    hold on; 
    frameDomain = (1:1:numFrames);
    plot(frameDomain, dat2plot, 'color', [0.5 0.5 0.5]);
    
    
    %% For each trial, draw a rectangle corresponding to the trial epoch over the activity traces and color-code it according to trial type:
 
    numTrials = size(trials, 1);
    trialDur = 2; % in seconds; currently hard-coded, should probably get this dynamically
    frameRate = 3.37; % in Hz; also currently hard-coded, should be found dynamically
    
    yl = ylim;
    ymin = yl(1);
    ymax = yl(2);
    
    width = trialDur * frameRate;
    height = ymax-ymin;
    
    numConditions = length(conditions);
    
    for i = 1:numTrials
        
        for j = 1:numConditions
            if strcmp(trials{i,2}, conditions{j}.Name)
                color = conditions{j}.Color;
                break
            end
        end
        
        r = rectangle('Position', [trials{i,1} ymin width height], 'FaceColor', color, 'EdgeColor', 'none');
        uistack(r, 'bottom');
    end
    
    set(gca, 'YTick', {}, 'YTickLabel', {});
    xlabel('Frame #');
    ylabel('dF/F');
    title('Grab overview');
    %close;
    
    status = exist(outputDirectory, 'dir');
    if status == 0
        mkdir(outputDirectory);
    end
    
    old = cd(outputDirectory);
    
    figPath = fullfile(outputDirectory, 'grab_overview.fig');
    saveas(gcf, figPath, 'fig');
    close;
    
    cd(old);
    
end