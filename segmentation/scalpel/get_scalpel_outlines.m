% get_scalpel_outlines.m
%%
% DOCUMENTATION TABLE OF CONTENTS
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated ddk 2017-10-25


%% I. OVERVIEW:
% This script converts the output of the image segmentation software
% SCALPEL into a set of .txt files that can be imported into ImageJ to
% easily view ROIs superimposed over the source movie.


%% II. REQUIREMENTS: 
% 1) JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files


%% III. INPUTS:
% This script is not (yet) a function and thus takes no formal input
% arguments. Instead, users should specify in the first line of code below
% the path to a .mat file containing the 'A' output of SCALPEL step2 or
% above.This .mat file should contain a single (p x q)-by-n binary matrix
% called `A`, where p is the video width, q is the video height, and n is
% the number of ROIs.Each element of the matrix is 1 if and only if the
% corresponding pixel is part of the corresponding ROI, 0 if otherwise.


%% IV. OUTPUTS:
% This script is not (yet) a function and thus has no formal return.
% However, it does save to secondary storage a directory of .txt files, one
% .txt file per ROI. Each .txt file contains a list of coordinates that
% defines the corresponding ROI, formatted into b rows and 2 columns, where
% b is the number of points that defines the boundary of the ROI, the first
% column represents the corresponding point's x-coordinate, and the second
% column represents the corresponding point's y-coordinate. 


%% Load data:
input_path = '/mnt/nas2/homes/dan/MultiSens/data/test_movies/5036-2_short/segmentation/Step1_21743/Step2_omega_0.2_cutoff_0.25/A.mat'; % DEFINE INPUT PATH HERE
[input_dir, name, ext] = fileparts(input_path);
cd(input_dir)
load(input_path); 
vid_height = 512;
vid_width = 512;
num_rois = size(A,2);


%% Compute binary mask of of ROI boundaries:

% We want to make a (vid_height-1)-by-(vid_width-1) matrix where the value is
% 0 if and only if the corresponding pixel in A is surrounded by 1's:
A = reshape(A, vid_height,vid_width,num_rois); % reshape as h x w x n slab, where n is the number of ROIs
A_clipped = A(2:vid_height-1,2:vid_width-1,:);

row_indices = (2:vid_height-1);
col_indices = (2:vid_width-1);

above_filled = A(row_indices-1,col_indices,:) == 1; % for each pixel in A (not along the border of the frame), is the pixel ABOVE filled in? 1 if yes, 0 if no
below_filled = A(row_indices+1,col_indices,:) == 1; % for each pixel in A (not along the border of the frame), is the pixel BELOW filled in? 1 if yes, 0 if no
left_filled = A(row_indices,col_indices-1,:) == 1; % for each pixel in A (not along the border of the frame), is the pixel TO THE LEFT filled in? 1 if yes, 0 if no
right_filled = A(row_indices,col_indices+1,:) == 1; % for each pixel in A (not along the border of the frame), is the pixel TO THE RIGHT filled in? 1 if yes, 0 if no
surrounded = above_filled & below_filled & left_filled & right_filled; % for each pixel in A (not along the border of the frame), is the pixel TOTALLY SURROUNDED by filled pixels (i.e., is it NOT part of the outline of an ROI)? 1 if yes, 0 if no

% Need to pad 'surrounded' matrix so it has the same dimensions as A (so they can be logically added):
col_zeros = zeros(size(surrounded,1),1,num_rois);
surrounded = [col_zeros surrounded col_zeros]; % Pad with left and right column of zeros:
row_zeros = zeros(1,size(surrounded,2),num_rois);
surrounded = [row_zeros; surrounded; row_zeros]; % Pad with top and bottom rows of zeros:

% Create 'outlines' matrix: 1 if the corresponding pixel in A is filled and NOT surrounded by other filled pixels, 0 otherwise
outlines = A & ~surrounded;


%% Write coordinates of ROI boundary pixels to .txt files:

% These will be added to outlines to get x and y coordinates of each outline pixel
cols = repmat((1:vid_width),vid_height,1); % matrix where every element is that element's column number (i.e., x-coordinate)
rows = repmat((1:vid_height)',1,vid_width); % matrix where every element is that element's row number (i.e., y-coordinate)

% Now we need to go through every plane and get and save the x,y coordinates of every non-zero in 'outlines':
mkdir('ROI_coords_txt');
cd ROI_coords_txt
for r = 1:num_rois
    
    % Get X and Y coordinates of each point that constitutes part of an ROI outline:
    X = cols(outlines(:,:,r)==1);
    Y = rows(outlines(:,:,r)==1);
    
    % Order them so that they'll render properly as a boundary in ImageJ:
    order_as_boundary = boundary(X,Y);
    X = X(order_as_boundary);
    Y = Y(order_as_boundary);
    
    % Writre the coordinates as a text file:
    fname = ['ROI_' num2str(r) '_coords.txt'];
    disp(fname)
    disp(class(fname))
    id = fopen(fname, 'wt');
    for p = 1:length(X)
        fprintf(id, [num2str(X(p)), ' ']);
        fprintf(id, [num2str(Y(p)), '\r\n']);
    end
    fclose(id);
    
    M.outputs(r).path = [input_dir filesep 'ROI_coords_txt' filesep fname];
    [status, out] = system(['sha1sum ' fname]);
    sha1 = out(1:40);
    M.outputs(r).sha1 = sha1;
end


%% Get and save metadata:

[status, input_sha1] = system(['sha1sum ' input_path]);
input_sha1 = input_sha1(1:40);

M.input.path = input_path;
M.input.sha1 = input_sha1;

savejson('',M,[input_dir filesep 'ROI_coords_txt' filesep 'metadata.json']);

%{
indices = ; 
above = indices-1;
below = indices+1;
left = indices-512;
right = indices+512;
%}