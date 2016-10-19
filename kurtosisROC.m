% Last updated DDK 2016-10-17

% OVERVIEW: 
% This function plots a receiver operating characteristic (ROC) curve that
% charcterizes the performance of using kurtosis as a criterion for
% accepting or rejecting automatically-segmented ROIs from a given grab.
% Its purpose is to help assess whether rejecting all ROIs whose pixel
% value distributions have a kurtosis beneath a certain threshold is an
% effective heuristic for filtering out false positives from
% automatically-identified ROIs.


% REQUIREMENTS:
% This function requires a set of automatically-segmented ROIs, each of
% which has been subsequently hand-annotated as either accepted as a real
% cell soma or rejected as note corresponding to a real cell soma.


% INPUTS:
% 1) rois - path to a .csv containing a K x 4 matrix of
% automatically-segmented ROIs from a given grab, sorted in descending
% order by kurtosis and hand-annotated as either accepted or rejected. This
% file should be of the same format returned by rankByKurtosis.m, i.e.:

% column 1: ROI number 
% column 2: kurtosis
% column 3: max value
% column 4: hand-annotated value: 1 for accepted, 0 for rejected

% Also note that the file is expected to have one header row, so the first
% row is skipped when read into this function.


% OUTPUTS: 
% This function has no formal return, but creates 2 plots: 1) an ROC curve
% obtained by using different kurtosis values as a threshold for deciding
% whether to accept or reject an automatically-segmented ROI, and 2) a plot
% of number of false negatives vs. kurtosis threshold. 

%%
function kurtosisROC(rois)
    kurtoses = xlsread(rois);
    
    kurtoses = flipud(kurtoses); %flip this so that kurtosis increases as we iterate through it
    FP = zeros(length(kurtoses),1); % vector of number of false positives for each possible kurtosis cutoff value
    TP = zeros(length(kurtoses),1); % vector of number of true positives for each possible kurtosis cutoff value

    % For each kurtosis value, get the number of true and false positives above
    % that kurtosis value:
    for i = 1:length(kurtoses)
        TP(i) = sum(kurtoses(i:end,4)); % number of true positives above i-th kurtosis level
        FP(i) = ((length(kurtoses)-i+1) - TP(i)); % number of false positives above i-th kurtosis level
    end

    % Get the total number of automatically-identified ROIs that correspond
    % to true cell soma (as judged by a human annotator) and the number
    % that don't
    totalPos = sum(kurtoses(:,4)); % total number of ROIs that actually correspond to a cell
    totalNeg = length(kurtoses) - totalPos; % total number of ROIs that don't actually correcpond to a cell

    FPR = FP./(totalNeg*ones(length(kurtoses),1)); % false positive rate (sum(false positive)/sum(condition negative))
    TPR = TP./(totalPos*ones(length(kurtoses),1)); % true positive rate (sum(true positive)/sum(condition positive))

    % Plot ROC curve:
    figure;
    plot(FPR, TPR);
    xl = xlim;
    yl = ylim;
    axMax = max(xl(2), yl(2));
    unity = (0:.1:axMax);
    hold on;
    plot(unity, unity);
    ylabel('True positive rate (a.u.)');
    xlabel('False negative rate (a.u.)');
    title('Receiver operating characteristic curve');

    % Plot false positives vs. kurtosis:
    figure;
    plot(kurtoses(:,2), FP);
    ylabel('False positives');
    xlabel('Kurtosis');
    title('False positives vs. kurtosis cutoff');

end