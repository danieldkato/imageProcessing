function frames = extract_frames(input_file, start_frame, end_frame, output_file)



%%
n_requested_frames = start_frame-end_frame+1;

% Get input and output extensions:
[dirname, fname, input_ext] = fileparts(input_file);
[dirname, fname, output_ext] = fileparts(output_file);

% Get basic file info:
if strcmp(input_ext, '.h5')
    info = h5info(input_file);
    height = info.Datasets(1).Dataspace.Size(1);
    width = info.Datasets(1).Dataspace.Size(2);
    n_frames_total = info.Datasets(1).Dataspace.Size(3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
end

subset_dims = [height width n_requested_frames];

% Check that the requested number of frames doesn't exceed the length of the movie:
if end_frame > n_frames_total
    error('Requested frames numbers exceed the number of frames in movie.');
end

% Read in the appropriate frames:
if strcmp(input_ext, '.h5')
    frames = h5read(ipnut_file, '/mov', [1 1 start_frame], subset_dims);
end

% Save output frames:
if strcmp(output_ext, '.h5')
    h5create(output_file, '/mov', subset_dims);
    h5write(output_file, '/mov', [1 1 1], subset_dims);
elseif strcmp(output_ext, '.tif')
end

