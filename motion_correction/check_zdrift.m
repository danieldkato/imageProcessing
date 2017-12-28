 function [r, sse] = check_zdrift(path)

% DOCUMENTATION TOC:
% I. OVERVIEW
% II. REQUIREMENTS
% III. USAGE
% IV. INPUTS
% V. OUTPUTS

% Last updated DDK 2017-12-26


%% I. OVERVIEW
% This function provides a metric of how much focal plane drift has
% occurred over the course of a movie. 

% It does this by taking the 2-D correlation coefficient between the
% average of the first 1000 frames and the aveerage of the last 1000
% frames. In order to control for any differences that might arise from X-Y
% translations, it first registers the images with the optimal 2D affine
% transformation before computing the correlation coefficient.


%% II. REQUIREMENTS
% MATLAB v >= ???


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


%% Find number of frames in movie:
tic;
disp('Getting file info...');
info = imfinfo(path);
disp('... done'); toc;
t = Tiff(path);
num_frames = numel(info);
width = info(1).Width;
height = info(1).Height;


%% Get average images of beginning and end of movie:

% Get average of first thousand frames:
tic;
disp('Computing mean of first 1000 frames...');
F = zeros(height, width, 1000);
for i = 201:1200
    t.setDirectory(i)
    F(:, :, i-200) = double(t.read());
end
avg_first = mean(F, 3); 
disp('... done.');
toc;

% Get average of last thousand frames:
tic;
L = zeros(height, width, 1000);
disp('Computing mean of last 1000 frames...');
for j = num_frames-999:num_frames
    t.setDirectory(j)
    L(:, :, j-(num_frames-1000)) = double(t.read());
end
disp('... done');
avg_last = mean(L, 3);
toc;


%% Transform images to compensate for any XY drift:

% Find best 2D transform between avg_first and avg_last:
tform = imregcorr(avg_last, avg_first);

% Register avg_last to avg_first:
rfixed = imref2d(size(avg_first));
avg_last_reg = imwarp(avg_last, tform, 'OutputView', rfixed);


%% Crop padded pixels due to translation:

% Recall that any translation will result in avg_last being padded with
% 0's, which will artefactually affect the correlation between avg_last and
% avg_first. We should omit these pixels from the correlation. To find the
% translation, look to the last row of the T field of the affine2d object
% tform. 

% Specifically, tform.T(3,1) is the horizontal translation, where a
% positive value means a shift to the right. A rightward shift means we
% have to omit the first ceil(tform.T(3,1)) columns. A leftward shift means
% we have to omit the last ceil(tform.T(3,1)) columns.

% tform.T(3,2) is the vertical translation, where a positive value means a
% shift down. A downward shift means we have to omit the first
% ceil(tform.T(3,2)) rows. An upwards shift means we have to omit the last
% ceil(tform.T(3,2)) rows.

if tform.T(3,1) > 0 
    first_col = 1 + ceil(tform.T(3,1));
    last_col = width;
elseif tform.T(3,1) < 0
    first_col = 1;
    last_col = width - ceil(abs(tform.T(3,1)));
end

if tform.T(3,2) > 0 
    first_row = 1 + ceil(tform.T(3,2));
    last_row = height;
elseif tform.T(3,2) < 0
    first_row = 1;
    last_row = height - ceil(abs(tform.T(3,2)));
end


avg_first = avg_first(first_row:last_row, first_col:last_col);
avg_last_reg = avg_last_reg(first_row:last_row, first_col:last_col);


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


%% Save outputs:

% Create an output structure including the R value, mean beginning image,
% mean end image, and difference
D.r = r;
D.sse = sse;
D.beginnig = avg_first_adjusted;
D.end = avg_last_adjusted;
D.diff_map = diff_map;
save('zdrift.mat', 'D');

% Save images:
imwrite(avg_first_adjusted, 'AVG_frames_201-1200_adjusted.tif');
imwrite(avg_last_adjusted, ['AVG_frames_' num2str(num_frames - 999) '-' num2str(num_frames) '_adjusted.tif']);
imwrite(diff_map, 'dif_map.tif');

%{
% Display images adjacent to each other;
figure();
imshowpair(avg_first, avg_last_reg, 'montage');
%}

%% Save metadata: