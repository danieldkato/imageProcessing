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
% This function requires a set of automatically-segmented ROIs that have
% been hand-annotated as either true positives or false positives.


% INPUTS:
% 1) rois - path to a .csv containing a K x 4 matrix of
% automatically-segmented ROIs from a given grab, sorted in descending
% order by kurtosis and hand-annotated as either true positives or false
% positives. This file should be of the same format returned by
% rankByKurtosis.m, i.e.:

% column 1: ROI number 
% column 2: kurtosis
% column 3: max value
% column 4: hand-annotated value: 1 for true positive, 0 for false positive

% Also note that the file is expected to have one header row, so the first
% row is skipped when read into this function.


% OUTPUTS: 
% This function has no formal return, but creates 2 plots: 1) an ROC curve
% obtained by using different kurtosis values as a threshold for deciding
% whether to accept or reject an automatically-segmented ROI, and 2) a plot
% of number of false negatives vs. kurtosis threshold. 

%%
function kurtosisROC(rois)
    kurtoses = csvread(rois, 1, 0);

    kurtoses = flipud(kurtoses); %flip this so that kurtosis increases as we iterate through it
    FP = zeros(length(kurtoses),1); % vector of number of false positives for each possible kurtosis cutoff value
    TP = zeros(length(kurtoses),1); % vector of number of true positives for each possible kurtosis cutoff value

    % For each kurtosis value, get the number of true and false positives above
    % that kurtosis value:
    for i = 1:length(kurtoses)
        TP(i) = sum(kurtoses(i:end,4));
        FP(i) = ((length(kurtoses)-i+1) - TP(i));
    end

    totalTP = sum(kurtoses(:,4));
    totalTN = length(kurtoses) - totalTP;

    FPR = FP./(totalTN*ones(length(kurtoses),1)); % false positive rate
    TPR = TP./(totalTP*ones(length(kurtoses),1)); % true positive rate

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