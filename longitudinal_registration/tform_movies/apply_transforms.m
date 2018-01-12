function Chunks = apply_transforms(Movies, mask)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:
% This function goes through each element of the struct array Movies,
% splits up the corresponding movie into chunks, then applies the
% transformation for the corresponding movie to each chunk in parallel.


%% II. REQUIREMENTS:
% 1) The MATLAB Parallel Computing toolbox


%% III. INPUTS:
% 1) Movies - 1 x m struct array where each element represents a movie to be
%    concatenated. For more detail on how this should be formatted, see the
%    OUTPUTS section of the documentation for initialize_movie_struct.m.

% 2) mask - m x n binary matrix specifying which pixels to clip, where m is
%    the width of the reference movie in pixels and n is the height of the
%    reference movie in pixels. For more detail, see the OUTPUTS section of
%    the documentation for get_movie_transformations.m.


%% IV. OUTPUTS:
% 1) Chunks - 1 x c array of structs where each element represents a single
%    chunk of the concatenated data. Each element includes the following
%    fields:

%   a) movie_number - int specifying the number of the movie corresponding to the chunk
%   
%   b) chunk_within_movie - int specifying the order of the chunk within
%   its corresponding movie (i.e. Chunks(c).chunk_within_movie = n iff
%   chunk c is the n-th chunk in its corresponding movie)

%   c) start_frame - int specifying the frame at which the chunk starts
%   relative to the first frame of its corresponding movie

%   d) end_frame - int specifying the frame at which the chunk ends
%   relative to the first frame of its corresponding movie

%   e) n_frames_in_chunk - int specifying the number of frames in the chunk

%   f) mfile - matfile object linked to a temporary .mat file where the
%      registered image data from the corresponding chunk will be stored


%% Initialize struct array `Chunks` here (before passing it to parfor, or else it will raise an error):

n_chunks = sum([Movies.n_chunks]);
ref_img = Movies(1).mean_img;
rfixed = imref2d(size(ref_img));

% 1 x m vector where m is the number of movies, and each element is the
% first chunk number corresponding to that movie:
chunk_offsets = cumsum([Movies.n_chunks]) - [Movies.n_chunks] +1; 

% Generate chunk names:
chunk_indices = (1:1:n_chunks);
chunk_names = arrayfun(@(x) pad(num2str(x),5,'left','0'), chunk_indices, 'UniformOutput', false);

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
    
    % Need to initialize some fields here or else parfor will throw an error:
    Chunks(c).mfile = [];
end


%% Apply transform to each chunk:
disp('Applying transformation to each chunk...');
parfor c = 1:n_chunks
        
    disp(['apply transform to chunk ' num2str(c)]);
    
    % Create a matfile object for the current chunk:
    chunk_name = chunk_names{c}
    Chunks(c).mfile = matfile(chunk_name, 'Writable', true);
    
    % Get the frames corresponding to the current chunk:  
    start_frame = Chunks(c).start_frame;
    end_frame = Chunks(c).end_frame;
    movie_number = Chunks(c).movie_number;
    var_name = Movies(movie_number).var_name;
    n_frames_in_chunk = Chunks(c).n_frames_in_chunk;
    img_dat = Movies(movie_number).matfile.(var_name)(:,:,start_frame:end_frame); 
    
    % Retrieve the transformation for the current chunk's corresponding
    % movie:
    transform = Movies(movie_number).transform; 
    
    % Initialize matrix in matfile where transformed image data will go:
    Chunks(c).mfile.M = uint8(zeros(size(ref_img,1), size(ref_img,2), n_frames_in_chunk)); 
    
    % Apply the transformation to every frame in the chunk:
    for frame = start_frame:end_frame
        frame_dat = Movies(movie_number).matfile.(var_name)(:,:,frame) % load the frame data
        tform_frame_dat = imwarp(frame_dat, transform, 'OutputView', rfixed); % apply the transformation
        tform_frame_dat(mask == 0) = 0; % apply clipping mask
        Chunks(c).mfile.M(:,:,frame - start_frame + 1) = tform_frame_dat; 
    end 
    
    disp(['done applying transform to chunk ' num2str(c)]);
   
end
disp(' ... done applying transform to each chunk.');

end