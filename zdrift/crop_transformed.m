function [fixed, moving, rows, cols] = crop_transformed(fixed, moving, tform)
% DOCUMENTATION TOC:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-03


%% I. OVERVIEW
% This function crops any pixels that lie outside the intersection of a
% pair of images one of which has been registered to the other.


%% II. REQUIREMENTS
% 1) MATLAB v >= ???


%% IV. INPUTS
% 1) fixed - m x n matrix of the fixed image (i.e., the reference image),
% where m is the image height and n is the image width

% 2) moving - m x n matrix of the moving image (i.e., the movie that is
% registered to the reference image), where m is the image height and n is
% the image width

% 3) tform - an affine2D object describing the transformation that
% optimally registeres the moving image to the fixed image.


%% V. OUTPUTS
% 1) fixed - m x n matrix of the fixed image (i.e., the reference image),
% where m is the image height and n is the image width

% 2) moving - m x n matrix of the moving image (i.e., the movie that is
% registered to the reference image), where m is the image height and n is
% the image width

% 3) rows - 1 x 2 vector specifying first and last row of original image
% included in cropped image

% 4) cols - 1 x 2 vector specifying first and last columns of original
% image included in cropped image


%%

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

height = size(fixed, 1);
width = size(fixed, 2);

%{
if tform.T(3,1) > 0 
    first_col = 1 + ceil(tform.T(3,1));
    last_col = width;
elseif tform.T(3,1) <= 0
    first_col = 1;
    last_col = width - ceil(abs(tform.T(3,1)));
end

if tform.T(3,2) > 0 
    first_row = 1 + ceil(tform.T(3,2));
    last_row = height;
elseif tform.T(3,2) <= 0
    first_row = 1;
    last_row = height - ceil(abs(tform.T(3,2)));
end
%}

rfixed = imref2d(size(fixed));
mask = ones(height, width);
mask_reg = imwarp(mask, tform, 'OutputView', rfixed);

first_row = 0;
last_row = height;
for c = 1:width
    uncropped = find(mask_reg(:,c) == 1);
    
    if min(uncropped) > first_row
        first_row = min(uncropped);
    end
    
    if max(uncropped) < last_row
        last_row = max(uncropped);
    end
end

first_col = 0;
last_col = width;
for r = 1:height
    uncropped = find(mask_reg(r,:) == 1);
    
    if min(uncropped) > first_col
        first_col = min(uncropped);
    end
    
    if max(uncropped) < last_col
        last_col = max(uncropped);
    end
end

fixed = fixed(first_row:last_row, first_col:last_col);
moving = moving(first_row:last_row, first_col:last_col);

rows = [first_row, last_row];
cols = [first_col, last_col];