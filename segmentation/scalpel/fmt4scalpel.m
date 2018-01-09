function fmt4scalpel(path)

% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated ddk 2017-10-24


%% I. OVERVIEW:
% This function re-formats 2-dimensional, single-color fluorescence movies
% for use with the SCALPEL automated image segmentation package. This
% entails re-shaping an m-by-n-by-t matrix as an (m x n)-by-t matrix and
% saving it a as a .mat, where m is the image height and n is the image
% width (in pixels), and t is the number of frames.


%% II. REQUIREMENTS:
% 1) JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files


%% III. INPUTS:
% 1) path - absolute path to a .mat file containing NoRMCorre-corrected
% imaging data.This file should contain only one object: an m-by-n-by-t
% matrix of motion-corrected imaging data, where m is the video height, n
% is the video width, and t is the number of frames.


%% IV. OUTPUTS:
% This function has no formal return, but saves the following to secondary
% storage:

% 1) Y_1.mat - a .mat file containing a single object, dat, which consists
% of an (m x n)-by-t matrix of motion-corrected imaging data, where m is
% the video height, n is the video width, and t is the number of frames.

% 2) metadata.json - a JSON file including absolute paths and SHA1
% checksums for the inputs used by and outputs generated by this function.


%% Load inputs: 
path = char(path); % need to convert from string to char
disp('Loading input...')
S = load(path);
disp('... done loading input.');
fnames = fieldnames(S);
dat = S.(fnames{1});
dims = size(dat);
nrows = dims(1);
ncols = dims(2);
nplanes = dims(3);

% Reshape the data:
dat = reshape(dat, [nrows*ncols, nplanes]);

% Save re-formatted data:
C = strsplit(path,"/");
grab_dir = C(1:end-3);
seg_dir = char(strjoin([grab_dir "segmentation"],'/'));
mkdir(seg_dir);
output_name = [seg_dir '/Y_1.mat'];
disp('Saving output...');
save(output_name, 'dat');
disp('... done saving output.');

% Write metadata:
M.inputs(1).path = path;
disp('Computing input sha1...');
[status, input_sha1] = system(['sha1sum ' path]);
disp('... done computing input sha1.');
M.inputs(1).sha1 = input_sha1(1:40);

M.outputs(1).path = output_name;
disp('Computing output sha1...');
[status, output_sha1] = system(['sha1sum ' output_name]);
disp('... done computing output sha1.');
M.outputs(1).sha1 = output_sha1(1:40);

% Save metadata:
metadata_path = [seg_dir '/Y_1_metadata.json'];
disp(metadata_path);

disp('Saving metadata...');
savejson('',M, metadata_path);
disp('... done saving metadata.');
