
% Create a plot for each ROI, including a mean activity trace (with
% shaded SEM areas) for each condition, for each session; perhaps
% traces from the same condition but different sessions can be
% different shades of the same color:

% For each session, create an axes object and plot one histogram for
% each condition:

% For each condition, create an axes object and plot one histogram for
% each session:

% Create an figure with s x c subplots, where s is the number of
% sessions and c is the number of pairwise combinations of conditions;
% each subplot will be a scatterplot with responses to one condition on
% one axis and responses to another conditions on another axis

% Load pre data

% Load post data

% For post and pre:
% make a vector of the mean/max/whatever across time for each condition
% make a scatterplot for each pair of vectors

% the only trick here is to ensure that the axes match

function preVpostCorr(pre, post, plotOrder, outputPath, condSettingsPath)
    %% Define a few variables that will be useful for organizing session-specific data later on:
    
    h5names = {pre, post}; % cell array containing paths to the HDF5 files to be loaded (one for each session)
    infos = cellfun(@(c) h5info(c), h5names); % struct array containing HDF5 metadata, including dataset names
    tables = cell(2,1); % cell array that will contain one table for each session
    sessionLabels = {'Pre', 'Post'}; 
    

    %% Load plot configuration:

    % plotOrder will be used to specify the order and orientation of the
    % plot axes. It can be either a cell array or loaded from a .txt. Check
    % which one of these plotOrder is:
    if ~iscell(plotOrder)
        fid = fopen(plotOrder);
        content = fscanf(fid, '%c');
        eval(content);
    end

    
    %% Validate data:
    
    % Get every unique abbreviation specified in the config file passed to the function:
    c = cat(1, plotOrder{:});
    disp('c');
    disp(c);
    tgtConditions = unique(c);
    disp('tgtConditions');
    disp(tgtConditions);
    numConditions = length(tgtConditions);
    
    % Make sure that each HDF5 file contains all of the conditions specified in plotOrder:
    for i = 1:2
        
        % Check that both pre- and post- HDF5s contain all of the conditions specified by the abbreviations in the config file:
        currAbbrevs = arrayfun(@(a) h5readatt(h5names{i}, strcat(['/',a.Name]), 'Abbreviation'), infos(i).Datasets, 'UniformOutput', 0);
        disp('currAbbrevs');
        disp(currAbbrevs);
        hits = cellfun(@(c) sum(strcmp(currAbbrevs, tgtConditions)), currAbbrevs);
        
        % If the number of matches between currAbbrevs and tgtAbbrevs is less than the number of tgtAbbrevs, raise an error:
        if sum(hits) < length(tgtConditions)
            error('Requested trial condition not found in input dataset. Make sure that requested trial condition abbreviations exactly match those of input datasets.'); 
        end
    end
    
    % Get the pre-stim period, post-stim period and stimulus duration from both HF5 files, and make sure that they are the same:
    timing = cell(4,3);
    timing{1,1} = 'Num pre-stimulus samples' ;
    timing{2,1} = 'Num post-stimulus samples';
    timing{3,1} = 'Num stim samples';
    timing{4,1} = 'Stimulus duration';
    
    for s = 1:2
        timing{1,s+1} = h5readatt(h5names{s},'/','num_samples_pre_stim');
        timing{2,s+1} = h5readatt(h5names{s},'/','num_samples_post_stim');
        timing{3,s+1} = h5readatt(h5names{s},'/','num_stim_samples');
        timing{4,s+1} = h5readatt(h5names{s},'/','stim_duration');
    end
    
    disp(timing);
    
    match = cellfun(@(a,b) isequal(a,b), timing(:,2), timing(:,3));
    if sum(match) < length(match)
        mismatches = find( match == 0 );
        for m = 1:length(mismatches)
            error(strcat([timing{mismatches(m),1}, ' does not match between requested sessions']));
        end
    end
    
    preStimSamples = timing{1,2};
    stimSamples = timing{3,2};
    stimRange = (preStimSamples+1:1:preStimSamples+stimSamples);
    
    disp('stimRange');
    disp(stimRange);
    
    
    %% For each session, copy the data from the HDF5 file over  into a table of the form:
    %{
                 Data             Means           Maxes           Color
          __________________  ______________  ______________  _____________
    'W'   [n x t x p double]  [n x 1 double]  [n x 1 double]  [1 x 3 double]
    'T'   [n x t x q double]  [n x 1 double]  [n x 1 double]  [1 x 3 double] 
    'W+T' [n x t x r double]  [n x 1 double]  [n x 1 double]  [1 x 3 double]
    
    where n is the number of ROIs being analyzed in the session, t is the
    number of frames in the peri-stimulus period, and p, q and r are the
    number of stimulus presentations in different conditions. 
    
    Tables are nice in this case because by making the row names the
    condition names, we can search the table by the condition abbreviations
    specified in plotOrder.
    %}
        
    for s = 1:2
        
        % Get the index of the HDF5 dataset corresponding to each condition named in tgtConditions:
        tgtConds2datasets = arrayfun(@(a) find(strcmp(tgtConditions, h5readatt(h5names{s}, strcat(['/',a.Name]),'Abbreviation'))), infos(s).Datasets)';
        
        
        % make sure that these indices are correct:
        disp('tgtConditions:');
        disp(tgtConditions);
        
        disp('dataset names:');
        infos(s).Datasets(:).Name
        
        disp('tgtConds2datasets');
        disp(tgtConds2datasets);
        
        % Gather the data from each condition into  c x 1 cell array, where c is the number of conditions; each element of c will be an n x t x p data slab:
        Data = arrayfun(@(a) h5read(h5names{s}, strcat(['/',infos(s).Datasets(a).Name])), tgtConds2datasets, 'UniformOutput', 0);
        
        % Compute some measure of the mean and max response of each ROI in each condition:   
        Means = cellfun(@(c) mean( mean(c, 3), 2), Data, 'UniformOutput', 0);
        %Maxes = cellfun(@(c) max( mean(c, 3), 1), Data, 'UniformOutput', 0);        
        
        % For each condition, get some measure of the max response of each cell:
        Maxes = cell(1,numConditions);
        for c = 1:numConditions
            
            meanTraces = mean(Data{c}, 3); % Get the mean activity trace over the peri-stimulus period for every ROI (an n x t matrix)
            condMaxes = zeros(size(meanTraces,1),1); % Initialize an n x 1 vector; each element will be the maximum excursion from baseline for a given cell in a given condition
            
            for n = 1:size(meanTraces,1)
                extremes = [min(meanTraces(n, stimRange)) max(meanTraces(n, stimRange))];
                [m, i] = max(abs(extremes));
                condMaxes(n) = extremes(i);
            end
            
            disp('condition:');
            disp(tgtConditions{c});
            
            disp('num ROIs:');
            disp(length(condMaxes));
            
            disp('condMaxes');
            disp(condMaxes);
            
            Maxes{c} = condMaxes;
        end
        
        % If condSettingsPath has been provided, get the color code from the condition settings file:
        if nargin > 4
            
            fid = fopen(condSettingsPath);
            content = fscanf(fid, '%c');
            eval(content);
            
            tgtConds2condSettings = arrayfun(@(a) find(strcmp(tgtConditions, a.Abbreviation)), Conditions, 'UniformOutput', 0); % Get the index of the Conditions entry corresponding to each condition named in tgtConditions:
            disp(cat(1,tgtConds2condSettings{:}));
            Colors = arrayfun(@(b) Conditions(b).Color, cat(2,tgtConds2condSettings{:}), 'UniformOutput', 0);
        
        % If condSettingsPath has not been provided, get the color code from the HDF5 attributes:
        else
            
            Colors = arrayfun(@(a) h5readatt(h5names{i}, strcat(['/',infos(s).Datasets(a).Name]), 'Color'), tgtConds2datasets, 'UniformOutput', 0);
        
        end
        
        % Assemble everything into a table:
        tables{s} = table(tgtConds2datasets', Data', Means', Maxes', Colors', 'RowNames', tgtConditions, 'VariableNames', {'Idx', 'Data', 'Means', 'Maxes', 'Color'});
        
    end
        
    
    %% Compute pairwise correlations specified by config file and store the results into a cell array of structures:
    
    % Create a cell array where correlation info will be stored: each row represents a session, and each column represents a pair of variables specified in plotOrder:
    corrs = cell(2, length(plotOrder));
    
    for s = 1:2
        for p = 1:length(plotOrder)
            
            pair = plotOrder{p};
            S.Xdat = tables{s}(pair(1),:).Maxes{:};
            S.Ydat = tables{s}(pair(2),:).Maxes{:};
            mdl = fitlm(S.Xdat,S.Ydat);
            
            S.Session = h5names{s};
            S.Xname = pair{1};
            S.Yname = pair{2};
            S.mdl = mdl;
            
            corrs{s,p} = S;
        end
    end
    
    % Convert the cell array corrs into a struct array (I find the latter easier to deal with in some ways):
    corrs = cell2mat(corrs);
    
    
    %% cd to appropriate directory:
    status = exist(outputPath, 'dir');
    if status == 0
        mkdir(outputPath);
    end
    old = cd(outputPath);
    
    
    %% For each session, create a histogram with distributions for each condition:
    histsBySession = figure;
    histsBySessionHandles = gobjects(2,numConditions);
    histsBySessionAxes = gobjects(2,1);
    %title('Pre- and post-pairing histograms of mean responses for each condition');
    
    for i = 1:2
        subplot(2,1,i);
        hold on;
        for c = 1:numConditions
            histsBySessionHandles(i,c) = histogram(tables{i}(tgtConditions(c),:).Maxes{:}, 'FaceColor', tables{i}(tgtConditions(c),:).Color{:});
        end
        histsBySessionAxes(i) = gca;
        legend(histsBySessionHandles(i,:), tgtConditions);
        xlabel('Mean dF/F');
        ylabel('Probability');
        v = text(0,0, sessionLabels{i});
        set(v, 'Units', 'normalized');
        set(v, 'Position', [-0.1 0.5]);
    end
    
    % Set normalization to probability rather than count:
    arrayfun(@(c) set(c, 'Normalization', 'probability'), histsBySessionHandles);
    
    disp(histsBySessionHandles);
    disp(length(histsBySessionHandles));
    
    % Standardize number of bins and minimum bin width across all histograms in figure (can't find an easy way to get this out of fig1histHandles in one line, will use for loop):
    binCounts = zeros(size(histsBySessionHandles,1)*size(histsBySessionHandles,2));
    binWidths = zeros(size(histsBySessionHandles,1)*size(histsBySessionHandles,2));
    
    for h = 1:size(histsBySessionHandles,1)*size(histsBySessionHandles,2)
        binCounts(h) = histsBySessionHandles(h).NumBins;
        binWidths(h) = histsBySessionHandles(h).BinWidth;
    end
    
    maxBins = max(max(binCounts));
    minWidth = min(min(binWidths(binWidths>0)));
    
    arrayfun(@(b) set(b, 'NumBins', maxBins), histsBySessionHandles, 'UniformOutput', 0);
    arrayfun(@(c) set(c, 'BinWidth', minWidth), histsBySessionHandles, 'UniformOutput', 0);
    
    % Standardize X and Y limits, bin counts and bin widths:
    maxVals = zeros(size(histsBySessionHandles,1)*size(histsBySessionHandles,2),1);
    
    for h = 1:size(histsBySessionHandles,1)*size(histsBySessionHandles,2)
        maxVals(h) = max(histsBySessionHandles(h).Values);
    end
    
    xMax = max(max(histsBySessionAxes.XLim));
    xMin = min(min(histsBySessionAxes.XLim));
    yMax = max(max(maxVals)) * 1.1;
    
    arrayfun(@(a) set(a, 'XLim', [xMin xMax], 'YLim', [0 yMax]), histsBySessionAxes, 'UniformOutput', 0);
    
    savefig('sessionHistograms.fig');
    
    
    %% For each condition, plot a histogram of pre- vs. post:
    condHistsFig = figure;
    hold on;
    
    for c = 1:numConditions
        
        subplot(1, numConditions, c);
        hold on;
        
        % Plot initial histograms:
        condHistHandles = cellfun(@(t) histogram(t(tgtConditions(c),:).Maxes{:}), tables, 'UniformOutput', 0);
        
        % Standardize axes, bin counts and bin widths:
        set(gca, 'XLim', [xMin xMax], 'YLim', [0 yMax]);
        cellfun(@(a) set(a, 'Normalization', 'probability'), condHistHandles);
        cellfun(@(b) set(b, 'NumBins', maxBins), condHistHandles);
        cellfun(@(c) set(c, 'BinWidth', minWidth), condHistHandles);
        
        % Write axis labels:
        ylabel('Probability');
        xlabel(tgtConditions{c});
        
        % Add legend:
        legend('Pre', 'Post'); % this assumes that pre is plotted first, followed by post; ok for now, but might want to think about including safeguard
    end
    
    savefig('conditionHistograms.fig');
    
    
    %% For each condition pair specified in plotOrder, plot a histogram of their ratio before pairing and a histogram of their ratio after pairing:
    
    % NOTE: this looks terrible
    
    figure;
    hold on;
    
    fig2histHandles = gobjects(1,2*length(plotOrder));
    fig2axesHandles = gobjects(1,length(plotOrder));
    
    colors = {[255/255 228/255 122/255],[255/255 204/255 0/255]};
    
    % Do initial plotting of histograms:
    for p = 1:length(plotOrder)
        subplot(1,length(plotOrder), p);
        hold on;
        
        for s = 1:2
            S = corrs(s,p);
            ratio = S.Ydat ./ S.Xdat;
            fig2histHandles( (p-1)*2 + s ) = histogram(ratio);
            fig2histHandles( (p-1)*2 + s ).Normalization = 'probability';
        end
        
        fig2axesHandles(p) = gca;
        xlabel(strcat([ '(', S.Yname, ')/', S.Xname]));
        ylabel('Probablility');
        
        legend('Pre','Post');
    end
    
    % Equalize scale of all axes objects for easy comparison:
    xMax = max(max([fig2axesHandles.XLim]));
    xMin = min(min([fig2axesHandles.XLim]));
    yMax = max(max([fig2axesHandles.YLim]));
    yMin = max(max([fig2axesHandles.YLim]));
    
    binCounts = zeros(size(fig2histHandles,1)*size(fig2histHandles,2));
    binWidths = zeros(size(fig2histHandles,1)*size(fig2histHandles,2));
    for h = 1:size(fig2histHandles,1)*size(fig2histHandles,2)
        binCounts(h) = fig2histHandles(h).NumBins;
        binWidths(h) = fig2histHandles(h).BinWidth;
    end
    
    maxBins = max(max(binCounts));
    minWidth = min(min(binWidths(binWidths>0)));
    
    %arrayfun(@(a) set(a, 'XLim', [xMin xMax]), fig2axesHandles, 'UniformOutput', 0);
    arrayfun(@(b) set(b, 'NumBins', maxBins), fig2histHandles, 'UniformOutput', 0);
    arrayfun(@(c) set(c, 'BinWidth', minWidth), fig2histHandles, 'UniformOutput', 0);
    
    savefig('ratioHistograms.fig');
    
    
    %% Compute and plot stats:
    
    corrFig = figure;
    hold on;
    
    
    % Create arrays to store axes handles:
    corrHandles = gobjects(2, length(plotOrder));
    corrAxes = gobjects(2, length(plotOrder));
    
    % 
    
    
    for p = 1:2
        
        %{
        % Find the indices of datasets corresponding to the conditions listed in tgtConditions:
        Indices = arrayfun(@(a) find(strcmp(tgtConditions, h5readatt(h5names{p}, strcat(['/',a.Name]),'Abbreviation'))), infos{p}.Datasets)'; 
        
        % Create a lookup table that maps a condition abbreviation to relevant data and metadata:
        Colors = arrayfun(@(a) h5readatt(h5names{p}, strcat(['/',infos{p}.Datasets(a).Name]), 'Color'), Indices, 'UniformOutput', 0);
        Data = arrayfun(@(a) h5read(h5names{p}, strcat(['/',infos{p}.Datasets(a).Name])), Indices, 'UniformOutput', 0);
        disp('Colors');
        disp(Colors);
        disp('Indices');
        disp(Indices);
        disp('Data');
        disp(Data);
        disp('tgtConditions:');
        disp(tgtConditions');
        LUT = table(Indices', Colors', Data', 'RowNames', tgtConditions, 'VariableNames', {'Idx', 'Color', 'Data'});

        disp(LUT);
        
        
        
        % For each condition, we want an n x 1 vector of mean responses,
        % where n is the number of ROIs. In other words, we want one data
        % point per ROI. Recall that for each ROI in a given condition, the
        % raw data has p x t data points, where p is the number of frames
        % in the peri-stimulus period and t is the number of trials of that
        % condition.
        
        % To condense the raw data for each condition to one point per ROI,
        % 1) take the mean activity trace across trials, then 2) take the
        % mean value over time of that mean activity trace.
        
        % Initialize two empty cell arrays that will become table columns:
        means = cell(numConditions,1);
        maxes = cell(numConditions,1);
        
        % Compute stats:
        for c = 1:numConditions
            % Get all the available data for the condition:
            dat = LUT(c,:).Data{:};
            means{c} = mean( mean(dat,3), 2 ); % compute the time-averaged mean of the mean activity trace
            maxes{c} = max( mean(dat,3), 2 ); % get the max of the mean activity trace
        end
        
        LUT.Means = means;
        LUT.Maxes = maxes;
        
        %annotation('textbox', [0, 0.25 .1 .1], 'String', 'test annotation');
        %}
        
        % Plot:
        for d = 1:length(plotOrder)
            
            % Select subplot:
            plotIdx = length(plotOrder)*(p-1)+d;
            subplot(2,3,plotIdx);
            hold on;
            
            % Retrieve the structure storing the correlation information for the current pairing:
            S = corrs(p,d);
            
            % Get the abbreviations for the variables to plot along X and Y for the current plot
            pair = plotOrder{d}; 
            
            % Use the abbreviations to retrieve data from the specified conditions for the current experiment:
            %X = tables{p}(pair(1),:).Maxes{:};  
            %Y = tables{p}(pair(2),:).Maxes{:};
            corrHandles(p,d) = scatter(S.Xdat,S.Ydat);
            corrAxes(p,d) = gca;
            
            set(gca, 'XLim', [min(S.Xdat)*1.1, max(S.Xdat)*1.1]);
            set(gca, 'YLim', [min(S.Ydat)*1.1, max(S.Ydat)*1.1]);
            
            xlabel(pair{1}, 'Color', tables{p}(pair(1),:).Color{:});
            ylabel(pair{2}, 'Color', tables{p}(pair(2),:).Color{:});
            
            %{



            
            mdl = fitlm(X,Y);
            B1 = ;
            B0 = ;
            R2 = ;
            pVal = ;
            [rho, pRho] = corr(X,Y);
            
            


            disp(LUT({pair{1}},:).Color);
            xlabel(pair{1}, 'Color', LUT({pair{1}},:).Color{:});
            xl = xlim;
            yl = ylim;
            disp(xl);
            ylabel(pair{2}, 'Color', LUT({pair{2}},:).Color{:});
            domain = (xl(1):.01:xl(2));
            disp(domain);
            bestFit = domain*B1 + B0;
            plot(domain, bestFit);
            xlim(xl);
            %}
            
            %{
            if d == 1
                text(xl(1)-.2, mean(ylim), 'test annotation');
            end
            
            text(xl(2)*.9, yl(2)*.9, sprintf(strcat(['m=', num2str(B1), '\n', 'R^2=', num2str(R2), '\n', 'p=', num2str(pVal)])));
            %}
            
        end
        
    end
    
    % Make everything the same scale:
    Xs = zeros(2, length(plotOrder));
    Ys = zeros(2, length(plotOrder));
    for s = 1:2
        for p = 1:length(corrAxes)
            Xs(:,(s-1)*length(plotOrder)+p) = corrAxes(s,p).XLim;
            Ys(:,(s-1)*length(plotOrder)+p) = corrAxes(s,p).YLim;
        end 
    end
    corrFigMinX = min(min(Xs));
    corrFigMaxX = max(max(Xs));
    corrFigMinY = min(min(Ys));
    corrFigMaxY = max(max(Ys));
    
    grandMin = min([corrFigMinX corrFigMinY]);
    grandMax = max([corrFigMaxX corrFigMaxY]);
    
    arrayfun(@(a) set(a, 'XLim', [grandMin grandMax]), corrAxes, 'UniformOutput', 0);
    arrayfun(@(b) set(b, 'YLim', [grandMin grandMax]), corrAxes, 'UniformOutput', 0);
    
    
    % Overlay axes, best fit line and correlation stats over each plot:
    domain = (grandMin:.01:grandMax);
    codomain = (grandMin:.01:grandMax);
    for s = 1:2
        for p = 1:length(plotOrder)
            
            axes(corrAxes(s,p));
            S = corrs(s,p);
            
            % Plot axes;
            x0 = plot(domain, zeros(1,length(domain)), 'Color', [0.85 0.85 0.85], 'LineWidth', 0.125);
            y0 = plot(zeros(1,length(codomain)), codomain, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.125);
            uistack(x0, 'bottom');
            uistack(y0, 'bottom');
            
            % Plot best fit line:
            plot(domain, domain*S.mdl.Coefficients{2,1} + mdl.Coefficients{1,1});
            
            % Print regression stats on plot:
            t = text(corrFigMaxX*.9, corrFigMaxY*.9, sprintf(strcat(['m = ', num2str(S.mdl.Coefficients{2,1}, 2), '\n', 'R^2 = ', num2str(S.mdl.Rsquared.ordinary, 2), '\n', 'p = ', num2str(S.mdl.Coefficients{2,4}, 2)])));
            set(t, 'Units', 'normalized')
            set(t, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'Position', [0.95 1]);        
        end
        
        % Label rows of plots:
        axes(corrAxes(s,1));
        u = text(0,0,strcat([sessionLabels{s}, ', n = ', num2str(length(corrs(s,1).Xdat))]));
        set(u, 'Units', 'normalized');
        set(u, 'Position', [-0.5 0.5], 'FontSize', 10)
    end
    
    savefig('correlations.fig');
    
    %{
    preDatasetNames = extractfield(preInfo.Datasets, 'Name');
    postDatasetNames = extractfield(postInfo.Datasets, 'Name');
    
    if ~isequal(preDatasetNames, postDatasetNames)
        
        % If they contain the same elements, but are just ordered differently, then re-order them so that they are the same:
        matches = cellfun(@(c) sum(strcmp(c, preDatasetNames)), postDatasetNames);
        if matches == length(info.preDatasetNames) && matches == length(info.postDatasetNames)
            
        
        elseif % If they don't even contain the same elements, then raise an error:
        end
    end
    %}
    
    %{
    {pre, post};
    
    for i = 1:2
        Conditions(1).Data = h5read();
    end
    %}
    
end