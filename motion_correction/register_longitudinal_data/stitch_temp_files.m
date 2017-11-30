function stitch_temp_files(Movies, Chunks, output_path)

output = matfile(output_path, 'Writable', true); % create master output file
n_frames_total = sum([Movies.n_frames]);
output.M = uint8(NaN(size(ref_img,1), size(ref_img,2), n_frames_total)); % initialize variable within master output file

% Create a 1 x m vector where m is the number of movies, and each element
% is the corresponding movie's start frame within the composite movie
frameoffsets = cumsum([Movies.n_frames]) - [Movies.n_frames] + 1;

for c = 1:n_chunks
    
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