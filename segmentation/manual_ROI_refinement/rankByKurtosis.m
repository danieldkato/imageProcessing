% Last updated DDK 2016-10-12

% OVERVIEW:
% This function takes activity traces of ROIs segmented from a given grab
% and sorts them by the kurtosis of their pixel values in descending order.
% It is intended as a heuristic in aiding the manual refinement of
% automatically-indetified ROIs, the idea being that
% automatically-identified ROIs with low kurtosis values may be more likely
% to be false positives.

% This function takes a K x T matrix of activity from a given grab, where K
% is the number of ROIs segmented from the grab and T is the number of
% frames in the grab, and returns a list of ROIs ordered by the kurtosis of
% their pixel values from highest to lowest. 


% INPUTS:
% 1) input - path to a K x T matrix of activity measured from K
% automatically-segmented ROIs, where T is the number of frames in the
% grab, saved as a .csv (e.g., the rawTraces.csv output by
% ca_source_extract.m).


% OUTPUTS: 
% 1) Kurtoses - a K x 3 matrix, where K is the number of ROIs automatically
% segmented from the input grab. The first column gives the ROI number, the
% second column gives the kurtosis of the ROI's pixel values, and the third
% column gives the ROI's maximum pixel value. 

% In addition to outputting Kurtoses as a formal return, this script also
% saves Kurtoses as a .csv in the directory from which it's run. The .csv
% appends two additional columns to Kurtoses: a blank column where a human
% annotator can indicate whether to accept or reject the ROI (for example,
% with a 1 or a 0, respectively), and another for any additional
% miscellaneous notes. The script also prepends a header row to the saved
% matrix for human readability.

%%
function Kurtoses = rankByKurtosis(input)
    
    dat = csvread(input)'; % dat will be t x n, where t is the number of frames and n is the number of ROIs
    disp(size(dat));
    numROIs = size(dat, 2);
    
    pxDistributions = cell(numROIs,1); % n x 1 cell array, where n is the number of ROIs; each cell will store a histogram of pixel values for an ROI
    Kurtoses = zeros(numROIs, 3); % n x 2 matrix, where n is the number of ROIs; the first column stores ROI indices, the second stores kurtosis values for the sample distribution of each ROI
    Kurtoses(:,1) = (1:1:numROIs);

    % Populate pxDistributions with histograms of pixel values for each ROI;
    % populate kurtoses with kurtosis values for the sample distribution for
    % each ROI
    for i = 1:numROIs
        pxDistributions{i} = histcounts(dat(:,i)); 
        Kurtoses(i, 2) = kurtosis(dat(:,i));
        Kurtoses(i, 3) = max(dat(:,i));
    end

    % Sort ROIs in order of descending kurtosis
    [Kurtoses, index] = sortrows(Kurtoses, -2);
    pxDistributions = pxDistributions(index);

    % Plot pixel value histograms of each ROI, and overlay its kurtosis
    numRows = floor(sqrt(numROIs));
    numCols = ceil(numROIs/numRows);
    for k = 1:numROIs
        hold on;
        subplot(numRows, numCols, k);
        plot(pxDistributions{k});
        set(gca, 'xtick', [], 'xticklabel', {});
        set(gca, 'ytick', [], 'yticklabel', {});
        title(strcat(['ROI ', num2str(Kurtoses(k, 1))]),'FontSize' , 8);
        text(0.1, 0.80, strcat(['k = ', num2str(Kurtoses(k, 2))]), 'Units', 'normalized', 'FontSize', 6);
    end

    % Plot histogram of kurtosis values:
    figure;
    histogram(Kurtoses(:,2));
    
    % Write header row:
    header = cell(1, 5);
    header{1, 1} = 'ROI #';
    header{1, 2} = 'kurtosis';
    header{1, 3} = 'max px value';
    header{1, 4} = 'annotation';
    header{1, 5} = 'notes';
    
    % Add header row to kurtosis-ordered ROI matrix and save to disk:
    fid = fopen('kurtoses.csv', 'w');
    fprintf(fid, '%s,', header{1,1:end-1});
    fprintf(fid, '%s\n', header{1,end});
    fclose(fid);
    
    dlmwrite('kurtoses.csv', Kurtoses(2:end,:), '-append');
end


