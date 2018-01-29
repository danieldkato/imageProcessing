function [fixed, moving] = blank_nonoverlapping_pixels(fixed, moving, tform)
% DOCUMENTATION TOC:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-29


%% I. OVERVIEW
% For a pair of images registered using the imref2d() function, this
% function replaces any non-overlapping pixels in each image with the pixel
% value of the image.

% When registering a moving image to a fixed image, the transformation may
% result in non-overlapping pixels; that is, there may be pixels in the
% fixed image that do not correspond to any pixel in the registered image
% (i.e., the output of the imwarp function). These non-overlapping pixels
% will have a value of 0 for all pixels in the fixed image, which will
% artificially change the result of the 2-D correlation function corr2().
% Making the non-overlapping pixels in each image the mean pixel value of
% the image will make it so that the non-overlapping pixels have no impact
% on the result of the 2-D correlation.


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


%%

% Define reference image:
height = size(fixed, 1);
width = size(fixed, 2);
rfixed = imreg2d(size(fixed));

% Define a mask of all 1's:
mask = ones(height, width); 

% Transform the mask using the same transformation that was used to
% register moving to fixed. This should reusult in a binary mask of all 0's
% in the pixels that don't overlap between moving and fixed.
mask_reg = imwarp(mask, tform, 'OutputView', rfixed); 

% Replace non-overlapping pixels of each image with the mean of each
% respective image:
fixed(~mask_reg) = mean(mean(fixed));
moving(~mask_reg) = mean(mean(moving));