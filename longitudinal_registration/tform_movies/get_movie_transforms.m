function [Movies, mask] = get_movie_transforms(Movies)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:
% This function goes through each element of the struct array Movies, each
% of which represents a movie to be concentenated by
% register_longitudinal_data, and computes the 2D affine transformation
% needed to register it to a reference image (namel,y the mean image of the
% first movie).


%% II. REQUIREMENTS:
% 1) The MATLAB parallel computing toolbox


%% III. INPUTS:
% 1) Movies - 1 x m struct array where each element represents a movie to be
%    concatenated. For more detail on how this should be formatted, see the
%    OUTPUTS section of the documentation for initialize_movie_struct.m.


%% IV. OUTPUTS:
% 1) Movies - same as the input, but now each element also includes a
%    non-empty 'transform' field, a MATLAB affine2d transformation object
%    specifying the transformation needed to register the corresponding
%    movie to the reference movie. For more detail on affine2d
%    transformation objects, see the MATLAB documentation on affined2d().

% 2) mask - m-by-n binary matrix where m is the width of the reference
%    movie in pixels and n is the height of the reference movie in pixels.
%    The value is 0 for pixels to be excluded from the final concatenated
%    movie and 1 for all others. This is to keep track of pixels that have
%    to be clipped (e.g., if a movie has to be translated down to be
%    registered to a reference movie, then the transformed movie B will
%    have a black stripe at the top). Note that this mask will ultimately
%    need to be applied to all movies, rather than just being
%    movie-specific; if a pixel is clipped from one movie, then it must be
%    excluded from ALL movies; e.g., if a transformed movie has a black
%    stripe at the top because it had to be translated down to be
%    registered to the reference movie, then we have to exclude the pixels
%    corresponding to the location of the black stripe from EVERY other
%    movie that we register.


%%
% Define the reference image; I'm arbitrarily making this the average image
% of the first movie (I don't think it makes a difference which is the
% reference):
ref_img = Movies(1).mean_img;

% Create a mask that will be used to keep track of pixels that have to be
% clipped (e.g., if a movie has to be translated down to be registered to a
% reference movie, then the transformed movie B will have a black stripe at
% the top);
mask = ones(size(ref_img));

% For each movie to be registered (i.e., every movie other than the
% reference movie), find the transformation needed to match the reference
% movie:
disp('Finding affine transformation for each movie..');
parfor m = 2:length(Movies)
    disp(['find movie ' num2str(m) ' transform']);
    Movies(m).transform = imregcorr(Movies(m).mean_img,ref_img); % estimate 2D affine transform for current movie
    Movies(m).registered_avg = imwarp(Movies(m).mean_img, Movies(m).transform, 'OutputView', rfixed); % apply transformation to average image of current movie
    
    % If a pixel is clipped from any movie, then it must be excluded from
    % ALL movies; e.g., if a transformed movie has a black stripe at the
    % top because it had to be translated down to be registered to the
    % reference movie, then we have to exclude the pixels corresponding to
    % the location of the black stripe from EVERY other movie that we
    % register
    mask = mask & Movies(m).registered_avg ~= 0;
    disp(['done finding movie ' num2str(m) ' transform' ]);
end

% In the special case of the reference movie, the transformation matrix is
% just the identity matrix:
Movies(1).transform = affine2d([1 0 0; 0 1 0; 0 0 1]);

disp('... done finding affine transformation for each movie.');

end