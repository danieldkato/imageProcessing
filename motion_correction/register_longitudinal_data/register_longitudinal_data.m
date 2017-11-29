% register_longitudinal_data

% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-28


%% I. OVERVIEW:
% This function concatenates longitudinal series of fluorescence movies
% into one long movie and transforms each constituent movie to compensate
% for slight day-to-day translations and rotations in the field of view.


%% II. REQUIREMENTS:
% 1) The MATLAB toolbox, JSONlab, available at https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files
% 2) The MATLAB function write_metadata.m, available at https://github.com/danieldkato/utilities/blob/master/write_metadata.m


%% III. INPUTS:



%% IV. OUTPUTS:
% This function has no formal return, but save to secondary storage a
% registered and concatenated 


%% TODO:


%% Load names of all motion-corrected data files - .mat? .tiff? Flexible?

tic;

params_file = '/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/motion_correction/register_longitudinal_data/register_longitudinal_params.json';
json_data = loadjson(params_file);
%files = json_data.movies_to_load; % load these as a cell array or something
files = cellfun(@(x) x.path, json_data.movies_to_load, 'UniformOutput', false);
data_var_names = cellfun(@(x) x.variable_name, json_data.movies_to_load, 'UniformOutput', false);
max_chunk_size = json_data.max_chunk_size; % load this from params file or something
output_file_name = json_data.output_name; % load this from some params file;

[output_directory, name, ext ] = fileparts(output_file_name);
old = cd(output_directory);


%% Get basic info and mean image for each movie:

n_movies = length(files);

% Get basic info for each movie; need to do define struct Movies, into
% which we will store mean image for each movie, BEFORE running parfor loop
% to actually compute mean image for each movie; MATLAB (somtimes) returns
% an error if parfor is given a struct array that doesn't exist yet
disp("Getting movie metadata...");
for m = 1:n_movies
    Movies(m).path = files{m};
    Movies(m).matfile = matfile(Movies(m).path);
    Movies(m).n_frames = size(Movies(m).matfile.(data_var_names{m}),3);
    Movies(m).n_chunks = ceil(Movies(m).n_frames/max_chunk_size);
    
    Movies(m).chunk_sizes = zeros(1, Movies(m).n_chunks);
    Movies(m).chunk_sizes(:) = max_chunk_size * ones(1, Movies(m).n_chunks);
    Movies(m).chunk_sizes(end) = mod(Movies(m).n_frames, max_chunk_size);
    
    Movies(m).chunk_start_frames = cumsum(Movies(m).chunk_sizes) - Movies(m).chunk_sizes +1;
    
    % pre-allocate some things here to prevent parfor error messages later
    % on:
    Movies(m).mean_img = []; 
    Movies(m).transform = [];
    Movies(m).registered_avg = [];
end
disp("... done getting movie metadata.");


% Compute mean image for each movie:
disp('Finding average image for each movie...');
parfor m = 1:n_movies
    disp(['finding average image for movie ' num2str(m)]);
    Movies(m).mean_img = mean(Movies(m).matfile.(data_var_names{m}),3);
    disp(['done finding average image for movie ' num2str(m)]);
end
disp('... done finding average image for each movie.');


%% For each movie, estimate the affine transformation needed to match the movie to a reference movie:

% Define the reference image; I'm arbitrarily making this the average image
% of the first movie (I don't think it makes a difference which is the
% reference):
ref_img = Movies(1).mean_img;
rfixed = imref2d(size(ref_img));

% Create a mask that will be used to keep track of pixels that have to be
% clipped (e.g., if a movie has to be translated down to be registered to a
% reference movie, then the transformed movie B will have a black stripe at
% the top);
master_mask = ones(size(ref_img));

% For each movie to be registered (i.e., every movie other than the
% reference movie), find the transformation needed to match the reference
% movie:
disp('Finding affine transformation for each movie..');
parfor m = 2:n_movies
    disp(['find movie ' num2str(m) ' transform']);
    Movies(m).transform = imregcorr(Movies(m).mean_img,ref_img); % estimate 2D affine transform for current movie
    Movies(m).registered_avg = imwarp(Movies(m).mean_img, Movies(m).transform, 'OutputView', rfixed); % apply transformation to average image of current movie
    
    % If a pixel is clipped from any movie, then it must be excluded from
    % ALL movies; e.g., if a transformed movie has a black stripe at the
    % top because it had to be translated down to be registered to the
    % reference movie, then we have to exclude the pixels corresponding to
    % the location of the black stripe from EVERY other movie that we
    % register
    master_mask = master_mask & Movies(m).registered_avg ~= 0;
    disp(['done finding movie ' num2str(m) ' transform' ]);
end

% In the special case of the reference movie, the transformation matrix is
% just the identity matrix:
Movies(1).transform = affine2d([1 0 0; 0 1 0; 0 0 1]);

disp('... done finding affine transformation for each movie.');


%% For every chunk, find the corresponding movie, apply the transformation for that movie to the chunk:

n_chunks = sum([Movies.n_chunks]);

% 1 x m vector where m is the number of movies, and each element is the
% first chunk number corresponding to that movie:
chunk_offsets = cumsum([Movies.n_chunks]) - [Movies.n_chunks] +1; 

% Generate chunk names:
chunk_indices = (1:1:n_chunks);
chunk_names = arrayfun(@(x) pad(num2str(x),5,'left','0'), chunk_indices, 'UniformOutput', false);

% Initialize struct array `Chunks` here (before passing it to parfor, or
% else it will raise an error):
for c = 1:n_chunks
    % Get the corresponding movie, start frame number, and stop frame number for the current chunk:
    movie_number = find(chunk_offsets<=c, 1, 'last'); % get the movie corresponding to the current chunk
    chunk_within_movie = c - chunk_offsets(movie_number) + 1; % get the current chunk's order within the movie; e.g., return n if the current chunk is the nth chunk within its corresponding movie
    chunk_size = Movies(movie_number).chunk_sizes(chunk_within_movie); 
    start_frame = Movies(movie_number).chunk_start_frames(chunk_within_movie); % get the current chunk's start frame within its corresponding movie
    end_frame = start_frame + chunk_size - 1; % get the current chunk's end frame within its corresponding movie
    n_frames_in_chunk = end_frame - start_frame + 1;
    
    Chunks(c).movie_number = movie_number; 
    Chunks(c).chunk_within_movie = chunk_within_movie;  
    Chunks(c).start_frame = start_frame; 
    Chunks(c).end_frame = end_frame; 
    Chunks(c).n_frames_in_chunk = n_frames_in_chunk;    
end

% Apply transform to each chunk:
disp('Applying transformation to each chunk...');
parfor c = 1:n_chunks
    
    disp(['apply transform to chunk ' num2str(c)]);
    
    % Create a matfile object for the current chunk:
    chunk_name = chunk_names{c}
    mfile = matfile(chunk_name, 'Writable', true);
    
    % Load the image data for the chunk:
    movie_number = Chunks(c).movie_number;
    data_var_name = data_var_names{movie_number};
    img_dat = Movies(movie_number).matfile.(data_var_name)(:,:,Chunks(c).start_frame:Chunks(c).end_frame); 
    
    % Retrieve the transformation for the current chunk's corresponding
    % movie:
    transform = Movies(movie_number).transform; 
    
    % Initialize array where transformed image data will go:
    mfile.M = uint8(zeros(size(ref_img,1), size(ref_img,2), n_frames_in_chunk)); 
    
    % Apply the transformation to every frame in the chunk:
    for frame = start_frame:end_frame
        frame_dat = Movies(movie_number).matfile.(data_var_name)(:,:,frame) % load the frame data
        disp(['class(frame_dat) = ' class(frame_dat)]);
        tform_frame_dat = imwarp(frame_dat, transform, 'OutputView', rfixed); % apply the transformation
        tform_frame_dat(master_mask == 0) = 0; % apply clipping mask
        disp(['class(mfile.M) = ' class(mfile.M)]);
        disp(['class(tform_frame_dat) = ' class(tform_frame_dat)]);
        mfile.M(:,:,frame - start_frame + 1) = tform_frame_dat; 
    end 
    
    disp(['done applying transform to chunk ' num2str(c)]);
   
end
disp(' ... done applying transform to each chunk.');


%% Stitch transformed data together into output file:

output = matfile(output_file_name, 'Writable', true); % create master output file
n_frames_total = sum([Movies.n_frames]);
output.M = uint8(NaN(size(ref_img,1), size(ref_img,2), n_frames_total)); % initialize variable within master output file

% Create a 1 x m vector where m is the number of movies, and each element
% is the corresponding movie's start frame within the composite movie
movieoffsets = cumsum([Movies.n_frames]) - [Movies.n_frames] + 1;

for c = 1:n_chunks
    
    movie_number = Chunks(c).movie_number; % get movie number
    mfile = matfile(chunk_names{c});
    
    % get absoulte start and end frames of current chunk:
    absolute_start_frame = movieoffsets(movie_number);
    absolute_end_frame = absolute_start_frame + Chunks(c).n_frames_in_chunk -1;
    
    % Write data from the temporary, chunk-specific matfile to the overall
    % output file:
    output.M(:,:,absolute_start_frame:absolute_end_frame) = mfile.M(:,:,:);
end


%% Write metadata:

duration = toc;
metadata_path = [output_directory filesep 'register_longitudinal_movies_metadata.json'];

% define inputs:
for m = 1:n_movies
    Metadata.inputs(m).path = files{m};
end

% define parameters:
Metadata.params.max_chunk_size = max_chunk_size;

% define outputs:
Metadata.outputs(1).path = output_file_name;

% define miscellaneous other metadata:
Metadata.duration = [num2str(duration) ' seconds'];

% write metadata:
write_metadata(Metadata, metadata_path);

% For each chunk, get a filename, chunk size, corresponding movie number,
% and frame start number within the movie
chunk_offsets = cumsum(n_chunks);