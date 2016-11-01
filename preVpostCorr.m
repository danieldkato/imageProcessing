
% Load pre data

% Load post data

% For post and pre:
% make a vector of the mean/max/whatever across time for each condition
% make a scatterplot for each pair of vectors

% the only trick here is to ensure that the axes match

function preVpostCorr(pre, post, config)
    %% Load plot configuration:

    % config will be used to specify the order and orientation of the
    % plot axes. It can be either a cell array or loaded from a .txt. Check
    % which one of these config is:
    if ~iscell(config)
        fid = fopen(config);
        content = fscanf(fid, '%c');
        eval(content);
    end

    
    %% Validate and load data:
    
    % Get every unique abbreviation specified in the config file passed to the function:
    c = cat(1, plotOrder{:});
    disp('c');
    disp(c);
    tgtConditions = unique(c);
    disp('tgtConditions');
    disp(tgtConditions);
    numConditions = length(tgtConditions);
    
    % Validate the data from each HDF5 file then load it into a table that can be searched by condition abbreviation:
    h5names = {pre, post};
    infos = cell(2,1);
    tables = cell(2,1);
    for i = 1:2
        
        % Check that both pre- and post- HDF5s contain all of the conditions specified by the abbreviations in the config file:
        infos{i} = h5info(h5names{i});
        currAbbrevs = arrayfun(@(a) h5readatt(h5names{i}, strcat(['/',a.Name]), 'Abbreviation'), infos{i}.Datasets, 'UniformOutput', 0);
        disp('currAbbrevs');
        disp(currAbbrevs);
        hits = cellfun(@(c) sum(strcmp(currAbbrevs, tgtConditions)), currAbbrevs);
        
        % If the number of matches between currAbbrevs and tgtAbbrevs is less than the number of tgtAbbrevs, raise an error:
        if sum(hits) < length(tgtConditions)
            error('Requested trial condition not found in input dataset. Make sure that requested trial condition abbreviations exactly match those of input datasets.'); 
        end
        
        % If the data are successfully validated, load the data from the HDF5 into a table:
        Indices = arrayfun(@(a) find(strcmp(tgtConditions, h5readatt(h5names{p}, strcat(['/',a.Name]),'Abbreviation'))), infos{p}.Datasets)'; 
        Colors = arrayfun(@(a) h5readatt(h5names{p}, strcat(['/',infos{p}.Datasets(a).Name]), 'Color'), Indices, 'UniformOutput', 0);
        
        tbl = ;
    end
    
    %% Load contents of each HDF5 into a table:
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Compute and plot stats:
    
    figure;
    hold on;
    
    
    % Create arrays to store axes handles:
    handles = NaN(numConditions*2,1);
    
    % 
    
    
    for p = 1:2
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
        
        % Plot:
        for d = 1:length(plotOrder)
            pair = plotOrder{d};
            X = LUT({pair{1}},:).Means{:};
            Y = LUT({pair{2}},:).Means{:};
            
            mdl = fitlm(X,Y);
            B1 = mdl.Coefficients{2,1};
            B0 = mdl.Coefficients{1,1};
            R2 = mdl.Rsquared.ordinary;
            pVal = mdl.Coefficients{2,4};
            [rho, pRho] = corr(X,Y);
            
            
            plotIdx = length(plotOrder)*(p-1)+d;
            subplot(2,3,plotIdx);
            hold on;
            handles(plotIdx) = scatter(X,Y);
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
            
            if d == 1
                text(xl(1)-.2, mean(ylim), 'test annotation');
            end
            
            text(xl(2)*.9, yl(2)*.9, sprintf(strcat(['m=', num2str(B1), '\n', 'R^2=', num2str(R2), '\n', 'p=', num2str(pVal)])));
            
        end
        
    end
    
    %minX = min(handles.XLim);
    %maxX = max(handles.X);
    
    
    
    
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
    
    {pre, post};
    
    for i = 1:2
        Conditions(1).Data = h5read();
    end
    
    
end