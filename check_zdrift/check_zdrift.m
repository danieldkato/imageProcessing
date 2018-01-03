 function [r, sse] = check_zdrift(path)

% DOCUMENTATION TOC:
% I. OVERVIEW
% II. REQUIREMENTS
% III. USAGE
% IV. INPUTS
% V. OUTPUTS

% Last updated DDK 2018-01-03


%% I. OVERVIEW
% This function provides a metric of how much focal plane drift has
% occurred over the course of a movie. 

% It does this by taking the 2-D correlation coefficient between the
% average of the first 1000 frames and the aveerage of the last 1000
% frames. In order to control for any differences that might arise from X-Y
% translations, it first registers the images with the optimal 2D affine
% transformation before computing the correlation coefficient.


%% II. REQUIREMENTS
% 1) The MATLAB function write_metadata.m, available at
% https://github.com/danieldkato/utilities/blob/master/metadata/MATLAB/write_metadata.m.
% Note that this function itself has a few dependencies.


%% III. USAGE
% This function can be called form within another MATLAB script or from the
% MATLAB command window. Alternatively, it can be invoked directly from the
% command line as follows:

% matlab -nodisplay -nosplash -r "check_zdrift <path\to\movie>"

% where <path\to\movie> stands for the path to the movie to be analyzed.


%% IV. INPUTS
% 1) path - path to a movie saved as a multi-page TIFF.


%% V. OUTPUTS
% 1) r - 2D correlation coefficient between the average of the first 1000
% frames and the average of the last 1000 frames of the input movie. Note
% that the average of the last 1000 frames has been registered to the
% average of the first 1000 frames using the optimal 2D affine
% transformation to eliminate any differences due to X-Y translation.

% In addition to formally returning r, this function saves the following to
% secondary storage:

% 1) D - a MATLAB struct with the following fields:
%   a) r - 2D correlation coefficient between the average of the first
%   1000 frames and the registered average of the last 1000 frames.
%   b) beginning - m x n matrix containing the mean image of the first 1000
%   frames of the input movie, where m and n are the movie height and
%   width, respectively. The gray values have been adjusted to maximize the
%   dynamic range.
%   c) end - m x n matrix containing the mean image of the last 1000
%   frames of the input movie, where m and n are the movie height and
%   width, respectively. The gray values have been adjusted to maximize the
%   dynamic range.
%   d) diff_map - m x n matrix containing the difference between the mean
%   image of the first 1000 frames and the mean of the last 1000 frames.
%   The gray values have been adjusted to maximize the dynamic range.

% 2) A TIFF of the mean of the first 1000 frames of the input movie. The
% gray values have been adjusted to maximize the dynamic range.

% 3) A TIFF of the mean of the last 1000 frames of the input movie. The
% gray values have been adjusted to maximize the dynamic range.

% 4) A TIFF of the difference between the mean of the first 1000 frames and
% the last 1000 frames of the input movie. The gray values have been
% adjusted to maximize the dynamic range.


%% TOODO:

% 1) Provide ways of dealing with short (<2000 frames long) movies?


%% CD to directory of movie:
[directory, name, ext] = fileparts(which(path));
cd(directory);


%% Find number of frames in movie:
tic;
disp('Getting file info...');
info = imfinfo(path);
disp('... done'); toc;
num_frames = numel(info);
width = info(1).Width;
height = info(1).Height;


%% Get average images of beginning and end of movie:

image_data = Tiff(path);

% Load first thousand frames and get average:
tic;
disp('Computing mean of first 1000 frames...');
F = zeros(height, width, 1000);
for i = 201:1200
    image_data.setDirectory(i)
    F(:, :, i-200) = double(image_data.read());
end
avg_first = mean(F, 3); 
disp('... done.');
toc;

% Load last thousand frames and get average:
tic;
L = zeros(height, width, 1000);
disp('Computing mean of last 1000 frames...');
for j = num_frames-999:num_frames
    image_data.setDirectory(j)
    L(:, :, j-(num_frames-1000)) = double(image_data.read());
end
disp('... done');
avg_last = mean(L, 3);
toc;


%% Transform average of last 1000 frames to compensate for any XY drift:

% Find best 2D transform between avg_first and avg_last:
tform = imregcorr(avg_last, avg_first);

% Register avg_last to avg_first:
rfixed = imref2d(size(avg_first));
avg_last_reg = imwarp(avg_last, tform, 'OutputView', rfixed);

% Crop any pixels that don't appear in both avg_first and avg_last_reg:
[avg_first, avg_last_reg] = crop_translated(avg_first, avg_last_reg, tform);


%% Quantify overlap of mean images:

% Find 2D correlation coefficient:
r = corr2(avg_first, avg_last_reg);
disp(r);

% Sum of squared errors:
sse = sum(sum((avg_last_reg - avg_first).^2));
disp(sse);


%% Compute correlaton map:
%{
Fbar = mean(mean(avg_first));
Lbar = mean(mean(avg_last_reg));

F_error = avg_first - Fbar;
L_error = avg_last_reg - Lbar;

F_sq_error = F_error.^2;
L_sq_error = L_error.^2;

F_sse = sum(reshape(F_sq_error, [1, width*height]));
L_sse = sum(reshape(L_sq_error, [1, width*height]));

d = sqrt(F_sse * L_sse);

corr_map = F_error.*L_error;
%}

%% Visualize images:

% Adjust gray values to maximize dynamic range:
Fmax = max(max(avg_first));
avg_first_gray = avg_first/Fmax;
avg_first_adjusted = imadjust(avg_first_gray);

Lmax = max(max(avg_last_reg));
avg_last_gray = avg_last_reg/Lmax;
avg_last_adjusted = imadjust(avg_last_gray);

% Visualize difference map:
diff_map = avg_last_adjusted - avg_first_adjusted;
diff_map = diff_map - min(min(diff_map));
diff_map = diff_map/max(max(diff_map));

%{
corr_map = corr_map - min(min(corr_map));
corr_map = corr_map/max(max(corr_map));
%}

%  Display images:
figure();
imshow(avg_first_adjusted);
figure();
imshow(avg_last_adjusted);
figure();
imshow(diff_map);

%{
figure();
imshow(corr_map);
%}


%% Try to register the mean images to the z-stack:

disp('Attempting to register average images to z-stack...');

% Find if there's a z-stack directory associated with the imaging site:

% CD to the site directory ad get contents:
movie_dir = cd('../..'); % assuming a directory structure of /mouse/session/site/grab/2P/movie.tif, this should cd to the site directory
ls = dir();

% Check if directory contents include zstack directory:
names = arrayfun(@(x) x.name, ls, 'UniformOutput', false);
is_zdir = cellfun(@(c) regexp(c, 'zstack'), names, 'UniformOutput', false);
zstack_exists = ~isempty(cell2mat(is_zdir));

% If zstack exists, try to register beginning and end average images to plane of z-stack:
if zstack_exists
    site_dir = cd('zstack');
    
    % Get the name of the actual z-stack file:
    ls = dir();
    names = arrayfun(@(x) x.name, ls, 'UniformOutput', false);
    is_zstack = cellfun(@(c) ~isempty(regexp(c, 'file_[0-9]*.tif', 'ONCE')), names, 'UniformOutput', false);
    
    % TODO: deal with situations where there's more than one regexp match?
    % TODO: deal with situations where there's no regexp match?
    z_idx = find(cell2mat(is_zstack)); 
    zstack_name = ls(z_idx).name;
    
    % Get z-stack metadata:
    z_metadata = loadjson('metadata.json'); 
    % TODO: throw warning if metadata doesn't exist
    % TODO: validate metadata
    disp('Getting z-stack metadata...')
    zinfo = imfinfo(zstack_name);
    num_zframes = length(zinfo);
    frames_per_slice = z_metadata.framesPerSlice;  
    disp('... done.')
    
    % Load the z-stack image data:
    disp('Loading z-stack image data...');
    z_tiff = Tiff(zstack_name);
    Z = NaN(zinfo(1).Height, zinfo(2).Width, num_zframes);
    for slice = 1:num_zframes
        z_tiff.setDirectory(slice);
        Z(:,:,slice) = z_tiff.read(); 
    end
    disp('... done');
    
    % If the z-stack inclues more than one frame per slice, average together frames
    % from each slice
    ('Averaging frames from the same slice...');
    if frames_per_slice > 1
        if mod(num_zframes, frames_per_slice) == 0 % confirm that the recorded number of frames per slice divides evenly into the number of frames in the z-stack
            
            num_slices = num_zframes/frames_per_slice;
            Z_avg = NaN(zinfo(1).Height, zinfo(2).Width, num_slices);
            
            ssi = frames_per_slice * ones(1, num_slices);
            slice_start_indices = cumsum(ssi) - slices_per_frame + 1;
            
            for i = 1:length(slice_start_indices)
                Z_avg(:,:,i) = mean(Z(:,:,slice_start_indices(i):start_slice_indices(i)+frames_per_slice-1),3);
            end
            
            disp('... done');
            
        else
            warn('Frames per slice must divide evenly into number of frames in z-stack; check that frames per slice recorded in z-stack metadata file is correct; skipping registration of average images to z-stack');
        end
        
    end
    
    
end


%% Save outputs:

% Check if output directory exists, and if not, create it and cd into it:
ls = dir();
names = arrayfun(@(x) x.name, ls, 'UniformOutput', false);
is_cz = cellfun(@(c) regexp(c, 'check_zdrift'), names, 'UniformOutput', false);
cz_exists = ~isempty(cell2mat(is_cz));
if ~cz_exists
    mkdir('check_zdrift');
end
old = cd('check_zdrift');

% Create an output structure including the R value, mean beginning image,
% mean end image, and difference
D.r = r;
D.sse = sse;
D.beginning = avg_first_adjusted;
D.end = avg_last_adjusted;
D.diff_map = diff_map;
save('zdrift.mat', 'D');

% Save average images:
f_name = 'AVG_frames_201-1200_adjusted.tif';
l_name = ['AVG_frames_' num2str(num_frames - 999) '-' num2str(num_frames) '_adjusted.tif'];

imwrite(avg_first_adjusted, f_name);
imwrite(avg_last_adjusted, l_name);
imwrite(diff_map, 'dif_map.tif');


%% Save metadata:

Metadata.inputs(1).path = which(path);
Metadata.outputs(1).path = which(f_name);
Metadata.outputs(2).path = which(l_name);

write_metadata(Metadata,'check_zdrift_metadata.json');

cd('..');

%{
% Display images adjacent to each other;
figure();
imshowpair(avg_first, avg_last_reg, 'montage');
%}