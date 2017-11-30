function Chunks = apply_transforms(Movies, mask)
%% Initialize struct array `Chunks` here (before passing it to parfor, or else it will raise an error):

n_chunks = sum([Movies.n_chunks]);

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