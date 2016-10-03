% Last updated DDK 2016-09-28

% OVERVIEW: 
% This function is a wrapper for doing image segmentation using Eftykios
% Pnevmatikakis' ca_source_extraction package. For documentation, see:

% https://github.com/epnev/ca_source_extraction/blob/master/documentation.pdf.

% The purpose of this wrapper is to: 
% 1) save the outputs in a standardized format that will be compatible with
% subsequent processing stages in the modular workflow described here:

% 10.112.43.46\Public\dank\multiSens\analysis\README.txt.

% 2) to forward and append similarly standardized metadata, also described
% in the link above.


% REQUIREMENTS:
% 1) ca_source_extraction, available at https://github.com/epnev/ca_source_extraction
% 2) coor2txt.m, available at https://github.com/danieldkato/image_segmentation/blob/master/ca_source_extraction/coor2txt.m
% 3) writeMetadata.m, available at https://github.com/danieldkato/metadata


% INPUTS:
% 1) inputPath - path to a motion-corrected multi-page TIFF stack to be segmented
% 2) outputPath - path where the output of the segmentation should be saved
% 3) K - number of components to be found
% 4) tau - std of gaussian kernel (size of neuron)
% 5) order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
% 6) merging threshold


% OUTPUTS:
% This script does not formally return anything, but it creates and saves
% four outputs in the directory specified by outputPath:

% 1) rawTraces.csv - a K x T matrix of demixed activity traces for each
% ROI, where K is the number of ROIs and T is the number of frames in the
% grab.
% 
% 2) Coor.mat - a K x 1 cell array containing the coordinates of every
% identified ROI. Each cell is a 2 x n matrix of coordinates defining
% the contours of the ROI, where n is the number of points used to
% define its contours. 

% 3) Coor2txt - a directory containing one text file for each ROI
% identified by ca_source_extraction. Each text file contains an n x 2
% matrix where n is the number of points that define the ROI. This can be
% used to load and visualize the ROIs in ImageJ, either for comparison to
% to other segmentation methods or for manual refinement.
%
% 4) meta.txt - a metadata file.

function ca_source_extract(inputPath, outputPath, K, tau, p, merge_thr)
    %% load file 
    
    addpath(genpath('../../ca_source_extraction'));

    nam = inputPath;                % insert path to tiff stack here
    sframe=1;						% user input: first frame to read (optional, default 1)
    num2read=2000;					% user input: how many frames to read   (optional, default until the end)

    Y = bigread2(nam,sframe,num2read);
    Y = Y - min(Y(:)); 
    if ~isa(Y,'single');    Y = single(Y);  end         % convert to single

    [d1,d2,T] = size(Y);                                % dimensions of dataset
    d = d1*d2;                                          % total number of pixels

    %% Set parameters
    
    options = CNMFSetParms(...                      
        'd1',d1,'d2',d2,...                         % dimensions of datasets
        'search_method','ellipse','dist',1,...      % search locations when updating spatial components %DK 160914 changed from 3 to 2
        'deconv_method','constrained_foopsi',...    % activity deconvolution method
        'temporal_iter',2,...                       % number of block-coordinate descent steps 
        'fudge_factor',0.98,...                     % bias correction for AR coefficients
        'merge_thr',merge_thr,...                    % merging threshold 
        'gSig',tau...
        );
    %% Data pre-processing

    [P,Y] = preprocess_data(Y,p);
    
    %% fast initialization of spatial components using greedyROI and HALS

    [Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize
    
    % display centers of found components
    Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
    figure;imagesc(Cn);
        axis equal; axis tight; hold all;
        scatter(center(:,2),center(:,1),'mo');
        title('Center of ROIs found from initialization algorithm');
        drawnow;

    %% manually refine components (optional)
    refine_components = false;  % flag for manual refinement
    if refine_components
        [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
    end

    %% update spatial components
    Yr = reshape(Y,d,T);
    clear Y;
    [A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

    %% update temporal components
    P.p = 0;    % set AR temporarily to zero for speed
    [C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    %% merge found components
    [Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

    %%
    display_merging = 1; % flag for displaying merging example
    if and(display_merging, ~isempty(merged_ROIs))
        i = 1; %randi(length(merged_ROIs));
        ln = length(merged_ROIs{i});
        figure;
            set(gcf,'Position',[300,300,(ln+2)*300,300]);
            for j = 1:ln
                subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
                    title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
            end
            subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                    title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
            subplot(1,ln+2,ln+2);
                plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
                hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
                title('Temporal Components','fontsize',16,'fontweight','bold')
            drawnow;
    end

    %% repeat
    P.p = p;    % restore AR value
    [A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
    [C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);

    %% do some plotting

    [A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
    K_m = size(C_or,1);
    [C_df,~] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],K_m+1); % extract DF/F values (optional)

    %contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
    figure;
    [Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),options,1); % contour plot of spatial footprints
    %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
    %% display components

    plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options)

    %% make movie

    make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)

    %% save ROI coordinates
    
    status = exist(outputPath, 'dir');
    if status == 0
        mkdir(outputPath);
    end 

    cd(outputPath);
    csvwrite('rawTraces.csv', C);
    save('Coor.mat', 'Coor');
    coor2txt(Coor);
    
    %% print metadata
    
    inputs = {{'motion-corrected TIFF', inputPath}};
    outputs = {{'segmented raw traces', strcat([outputPath, '\rawTraces.mat'])};
               {'ROI coordinates (.mat)', strcat([outputPath, '\Coor.mat'])};
               {'ROI coordinates (.txt)', strcat([outputPath, '\Coor2txt\'])}
              };
    params = {{'K', K};
              {'tau', tau};
              {'p', p};
              {'merge_thr', merge_thr}
             };
   
    writeMetadata('segmentation', 'ca_source_extraction', inputs, outputs, params);
end