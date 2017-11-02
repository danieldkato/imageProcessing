function tiff_2_mat(input_path, output_path)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-02

%% I. OVERVIEW:
% This function converts a multi-page TIFF into a w x h x t matrix and
% saves it into a .mat file, where w is the image width, h is the image
% height, and t is the number of frames in the TIFF. Note that in order to
% deal with larger files that can't be read into memory all at once, it
% reads in the TIFF frame-by-frame. 


%% II. REQUIREMENTS:
% None (beyond MATLAB v >= ?)


%% III. INPUTS:
% 1) input_path - path to a multi-page TIFF.
% 2) output_path - path where output .mat file should be saved.


%% IV. OUTPUTS:
% This fucntion has no formal return, but saves to secondary storage a .mat
% file containing a single w x h x t matrix called Y, where w is the image
% width, h is the image height, and t is the number of frames in the TIFF. 


%%
% Get data dimensions:
info = imfinfo(input_path);

% Initialize output object:
disp('Initializing memory-mapped output file...');
Y = matfile(output_path, 'Writable', true);
Y.Y = uint16(NaN(info(1).Width,info(2).Height,length(info)));
disp('... done.');

% Copy data from input TIFF to memory-mapped output file:
disp('Copying image data...');
for f = 1:length(info) 
    Y.Y(:,:,f) = imread(input_path,'Index',f);
end
disp('... done.');

end