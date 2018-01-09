function stitch_temp_files(Movies, Chunks, output_path)
%% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29

%% I. OVERVIEW:
% This function stitches together the transformed image data in the struct
% array Chunks into a single .mat file. 


%% II. REQUIREMENTS:
% None


%% III. INPUTS:
% 1) Movies - 1 x m struct array where each element represents a movie to be
%    concatenated. For more detail on how this should be formatted, see the
%    OUTPUTS section of the documentation for initialize_movie_struct.m.

% 2) Chunks - 1 x c array of structs where each element represents a single
%    chunk of the concatenated data. For more detail on how this should be
%    formatted, see the OUTPUTS section of the documentation for
%    apply_transforms.m. 

% 3) output_path - char array specifying the name of the .mat file to which
%    to save the registered and concatenated image data


%% IV. OUTPUTS:
% This function has no formal return, but saves to secondary storage a .mat
% file containing an m x n x t matrix of registered and concatenated image
% data, where m is the width of the reference movie in pixels, n is the
% height of the reference movie in pixels, and t is the sum of the number
% of frames across all movies. 


%%
n_frames_total = sum([Movies.n_frames]);
ref_img = Movies(1).mean_img;

output = matfile(output_path, 'Writable', true); % create master output file
output.M = uint8(NaN(size(ref_img,1), size(ref_img,2), n_frames_total)); % initialize variable within master output file

% Create a 1 x m vector where m is the number of movies, and each element
% is the corresponding movie's start frame within the composite movie
frameoffsets = cumsum([Movies.n_frames]) - [Movies.n_frames] + 1;

for c = 1:length(Chunks)
    
    movie_number = Chunks(c).movie_number; % get movie number corresponding to current chunk
    mfile = Chunks(c).mfile;
    chunk_start_frame = Chunks(c).start_frame;
    n_frames_in_chunk = Chunks(c).n_frames_in_chunk;
    
    % get absoulte start and end frames of current chunk:
    absolute_start_frame = frameoffsets(movie_number) + chunk_start_frame - 1;
    absolute_end_frame = absolute_start_frame + n_frames_in_chunk - 1;
    
    % Write data from the temporary, chunk-specific matfile to the overall
    % output file:
    output.M(:,:,absolute_start_frame:absolute_end_frame) = mfile.M(:,:,:);
end

end