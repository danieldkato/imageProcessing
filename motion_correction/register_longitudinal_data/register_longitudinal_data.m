% register_longitudinal_data
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:

% This function concatenates multiple fluorescence movies of the same
% imaging site into one long movie. It also transforms individual movies to
% compensate for slight day-to-day translations and rotations between
% movies.

% The purpose of this is to pre-process longitudinal datasets for
% activity-based segmentation alogrithms. These activity-based segmentation
% algorithms must be run on the entirety of the longitudnial datasets in
% order to yield accurate results; for example, if neuron N only becomes
% active in the second imaging session, and you run the segmentation on the
% first imaging session, then these activity-based segmentation algorithms
% will fail to include neuron N in any of the movies. 


%% II. REQUIREMENTS:

% 1) The MATLAB Parallel Computing Toolbox
% 2) The MATLAB toolbox JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files
% 3) The MATLAB function write_metadata.m, available at https://github.com/danieldkato/utilities/blob/master/write_metadata.m
% 4) The MATLAB function initialize_movie_struct.m, available at https://github.com/danieldkato/imageProcessing/tree/master/motion_correction/register_longitudinal_data
% 5) The MATLAB function get_mean_images.m, avialable at https://github.com/danieldkato/imageProcessing/blob/master/motion_correction/register_longitudinal_data/get_mean_images.m
% 6) The MATLAB function get_movie_transforms.m, available at https://github.com/danieldkato/imageProcessing/blob/master/motion_correction/register_longitudinal_data/get_movie_transforms.m
% 7) The MATLAB function apply_transforms.m, available at https://github.com/danieldkato/imageProcessing/blob/master/motion_correction/register_longitudinal_data/apply_transforms.m
% 8) The MATLAB function stitch_temp_files, available at https://github.com/danieldkato/imageProcessing/blob/master/motion_correction/register_longitudinal_data/stitch_temp_files.m


%% III. INPUTS:

% Currently none


%% IV. OUTPUTS:

% This function has no formal return, but saves the following to secondary
% storage:

% 1) a .mat file containing the data from all input movies concatenated into
%    one long movie. 

% 2) a metadata JSON file specifying input and output paths and SHA1
%    checksums and parameters.


%% TODO:

% 1) Make compatible with different input types (like TIFFs)
% 2) Support different output type options


%% Load parameters:
tic;
params_file = '/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/motion_correction/register_longitudinal_data/register_longitudinal_params.json';
json_data = loadjson(params_file);
movie_name_list = cellfun(@(x) x.path, json_data.movies_to_load, 'UniformOutput', false);
data_var_names = cellfun(@(x) x.variable_name, json_data.movies_to_load, 'UniformOutput', false);
max_chunk_size = json_data.max_chunk_size;
output_file_name = json_data.output_name; 
[output_directory, name, ext ] = fileparts(output_file_name);
old = cd(output_directory);


%% Find and apply 2D affine transforms for each movie:

% Initialize Movies struct; MATLAB (sometimes) returns an error if parfor
% is given a struct array that doesn't exist yet
Movies = initialize_movie_struct(movie_name_list, max_chunk_size);

% Compute mean image for each movie:
Movies = get_mean_images(Movies);

% For each movie, estimate the affine transformation needed to match the movie to a reference movie:
[Movies, mask] = get_movie_transforms(Movies);

% For every chunk, find the corresponding movie, apply the transformation
% for that movie to the chunk, and save the transformed chunk to a
% temporary .mat file:
Chunks = apply_transforms(Movies, mask);

% Stitch the transformed data together into one output file:
stitch_temp_files(Movies, Chunks, output_file_name);


%% Write metadata:

% define inputs:
for m = 1:length(Movies)
    Metadata.inputs(m).path = movie_name_list{m};
end

Metadata.params.max_chunk_size = max_chunk_size; % define parameters:
Metadata.outputs(1).path = output_file_name; % define outputs:
Metadata.duration = [num2str(toc) ' seconds']; % define miscellaneous other metadata:

% write metadata:
metadata_path = [output_directory filesep 'register_longitudinal_movies_metadata.json'];
write_metadata(Metadata, metadata_path);

% For each chunk, get a filename, chunk size, corresponding movie number,
% and frame start number within the movie
chunk_offsets = cumsum(n_chunks);