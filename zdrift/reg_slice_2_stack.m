function s = reg_slice_2_stack(slice, stack)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-29


%% I. OVERVIEW:
% This function takes a 2-dimensional image and determines which plane it a
% 3-dimensional stack it's most correlated with.


%% II. REQUIREMENTS:
% 1) MATLAB >= ???


%% III. INPUTS:
% 1) slice - h x w image matrix, where h is the image height in and w is
% the image width in pixels.

% 2) stack - h x w x z image matrix, where h is the image height and w is
% the image width in pixels, and z is the number of slices in the stack.


%% IV. OUTPUT:
% 1) s - index of the slice in stack most strongly correlated with input slice image.


%%
num_slices = size(stack, 3);
corrs = zeros(num_slices, 1);
for z = 1:num_slices
    corr(z) = corr2(slice, stack(:,:,z));
end

[m, s] = max(corrs);

