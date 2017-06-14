% load data
fName = 'dFF_parsed_trials.h5';
info = h5info(fName);
names = arrayfun(@(a) a.Name, [info.Datasets], 'UniformOutput', 0); % condition names
dSets = cellfun(@(c) h5read(fName, strcat('/', c)), names, 'UniformOutput',0); % load datasets into cell array

nConds = size(dSets, 1); % num conditions
n = size(dSets{1}, 1); % get num neurons
disp(strcat(['num neurons = ', num2str(n)]));

% calculate total num trials
t = 0;
trialNums = zeros(nConds,1);


%% prepare analysis where each condition is a different level of the same categorical factor; this CAN'T be used to ascertain intercative effects
condLabels = [];
for c = 1:length(names)
    names{c} = strrep(names{c}, ' ', []);
    nCondTrials = size(dSets{c},3);
    trialNums = nCondTrials;
    t = t + size(dSets{c}, 3);
    labels = zeros(nCondTrials, nConds);
    labels(:,c) = 1;
    condLabels = vertcat(condLabels, labels);
end
%disp(condLabels)

T1 = array2table(condLabels, 'VariableNames', names);
%disp(T);

% maxProj = squeeze(maxProj)

pVals1 = zeros(n, 3); % column 1: speaker, column 2: stepper and speaker, column 3: stepper only


%% preapre analysis that includes interactive effects
condLabels2 = zeros( size(condLabels,1), 2); % column 1: whether speaker is present, column 2: whether stepper is present
condLabels2(:,1) = (condLabels(:,1) | condLabels(:,2));
condLabels2(:,2) = (condLabels(:,3) | condLabels(:,2));
disp(condLabels2);
T2 = array2table(condLabels2, 'VariableNames', {'speakeronly','stepperonly'});
pVals2 = zeros(n, 3); % column 1: speaker, column 2: stepper, column 3: interaction


%% Do both analyses for each neuron
for i = 1:n
    dFF = [];
    for d = 1:length(dSets)
        mat = squeeze(dSets{d}(i, :, :))'; % f x s, where f is number of frames in peri-stim period and and s is number of stim presentations; since I'm doing max projection, filtering out pre-stim period probably isn't necessary
        %proj = max(mat, [], 2); 
        proj = mean(mat,2);
        %disp(size(maxProj));
        dFF = vertcat(dFF, proj);
    end
    %disp(dFF);

    % for the first analysis
    T1.dFF = dFF;
    %disp(T);
    out1 = fitlm(T1, 'dFF~-1+speakeronly + stepperandspeaker + stepperonly');
    %disp(out1);    

    pVals1(i, 1) = out1.Coefficients{1,4};
    pVals1(i, 2) = out1.Coefficients{2,4};
    pVals1(i, 3) = out1.Coefficients{3,4};
    
    % for the second analysis
    T2.dFF = dFF;
    out2 = fitlm(T1, 'dFF~-1+speakeronly*stepperonly');
    disp(out2);

end

%sum(pVals(:,3)>0);
