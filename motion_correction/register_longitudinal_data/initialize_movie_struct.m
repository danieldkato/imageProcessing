function Movies = initialize_movie_struct(files, max_chunk_size)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:
% This function initialized a struct array called Movies, each element of
% which represents a single movie to be concatenated by
% register_longitudinal_data.m. Each element contains information about
% chunk boundaries within the corresponding movie.


%% II. REQUIREMENTS:
% None


%% III. INPUTS:
% 1) files - a cell array of absolute paths to .mat files containing data
%    from motion-corrected movies to be concatenated.

% 2) max_chunk_size - int specifying maximum number of frames per 'chunk',
%    where 'chunks' will later be processed in parallel. At a later,
%    intermediate state of processing, each chunk will be saved as its own
%    temporary .mat file.


%% IV. OUTPUTS:
% 1) Movies, a struct array representing the movies to be concatenated.
%    Each element has the following fields:

%   a) path - char array specifying the path to the .mat file containing
%      the data for the movie

%   b) matfile - a matfile variable linked to the source .mat file

%   c) var_name - char array specifying the name of the variable in the
%      source .mat file where the image data is stored. It is assumed that
%      the .mat file contains only one variable, and that this is where the
%      image data is stored.

%   d) n_frames - int specifying number of frames in the corresponding movie

%   e) n_chunks - int specifying number of chunks into which the corresponding movie will
%      be divided

%   f) chunk_sizes - 1 x c int vector, where c equals n_chunks, and each
%      element specifies the number of frames in the corresponding chunk

%   g) chunk_start_frames - 1 x c int vector, where c equals n_chunks, and
%      each element specifies the start frame number of the corresponding
%      chunk within the movie 

% This function also initializes a few additional fields that are left
% empty for the time being (these will be populated by functions invoked
% later by register_longitudinal_data).


%%
disp("Getting movie metadata...");

n_movies = length(files);

for m = 1:n_movies
    Movies(m).path = files{m};
    Movies(m).matfile = matfile(Movies(m).path);
    
    details = whos(Movies(m).matfile);
    var_name = details(1).name;
    
    Movies(m).var_name = var_name;
    Movies(m).n_frames = size(Movies(m).matfile.(var_name),3);
    Movies(m).n_chunks = ceil(Movies(m).n_frames/max_chunk_size);
    
    disp(max_chunk_size)
    disp(Movies(m).n_chunks)
    
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

end