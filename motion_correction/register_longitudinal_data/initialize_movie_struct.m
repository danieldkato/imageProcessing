function Movies = initialize_movie_struct(files, max_chunk_size)

disp("Getting movie metadata...");

n_movies = length(files);

for m = 1:n_movies
    Movies(m).path = files{m};
    Movies(m).matfile = matfile(Movies(m).path);
    
    details = whos(Movies(m).matfile);
    var_name = details(1).name;
    
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