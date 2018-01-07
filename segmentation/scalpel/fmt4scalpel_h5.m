function fmt4scalpel_h5(varargin)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-05


%% I. OVERVIEW:
% This function reshapes a w x h x f matrix of imaging data from an input
% HDF5 file into a (w*h) x f matrix and saves it in another HDF5, where w
% is the video width, h is the video height, and f is the number of frames.
% This vectorized format is necessary for use with certain segmentation
% algorithms such as SCALPEL.


%% II. USAGE:

% In MATLAB, invoke this function with either of the following:

% fmt4scalpel_h5(input_path, output_path)
% fmt4scalpel_h5(input_path, output_path, chunk_size)
% fmt4scalpel_h5(input_path, output_path, chunk_size, multiple_outputs)

% In addition to invoking this function from another MATLAB script or from
% the MATLAB command line, it is possible to invoke this function from the
% LINUX command line with the following:

% matlab -nosplash -nodesktop -r "fmt4scalpel_h5 <input_path> <output_path> [<chunk_size>] [<multiple_outputs>]"


%% III. REQUIREMENTS:
% 1) MATLAB >= ???


%% IV. INPUTS:
% 1) input_path - path to an HDF5 file to be vectorized. The imaging data to be
% converted should be saved in a w x h x f dataset called '/mov', where w
% is the width of the movie in pixels, h is the height of the movie in
% pixels, and f is the number of frames in the movie. 

% 2) output_path - path where the output HDF5 file should be saved.  

% 3) chunk_size (optional) - integer specifying the number of frames to
% read into memory at a time. Default value is 1000.

% 4) multiple_outputs (optional) - Boolean flag specifying whether to save
% each chunk as a separate file (useful for several segmentation
% algorithms).


%% V. OUTPUT:
% This function has no formal return, but saves to disk HDF5 files
% containing vectorized imaging data. Within each HDF5 file, these data are
% saved in a dataset called '/mov'.

% If multiple_outputs is set to false (the default behavior), then all of
% the reshaped data will be saved in a single HDF5 file. The dataset '/mov'
% will have the dimensions w*h x f, where w is the width of the movie in
% pixels, h is the height of the movie in pixels, and f is the total number
% of frames in the movie. The name of the file is just the specified output
% path.

% If multiple_outputs is set to true, then the reshaped data will be saved
% in multiple HDF5 files. Within each file, the dataset '/mov' will have
% the dimensions w*h x c, where w is the width of the movie in pixels, h is
% the height of the movie in pixels, and c is the number of frames
% specified in the chunk_size argument. If output_path is specified as
% `/data/processed/Y.h5', then the output files will be named 'Y_1.h5',
% 'Y_2.h5', 'Y_3.h5', etc.


%% TODO:
% 1) Throw warning if chunk_size is greater than the length of the movie


%% Define parameters:
input_path = varargin{1};
output_basepath = varargin{2};
if nargin > 2
    chunk_size = varargin{3};
else
    chunk_size = 1000; % default chunk size
end
if nargin > 3
    multiple_outputs = varargin{4};
else
    multiple_outputs = false;
end

fs = filesep;
[dir, base, ext] = fileparts(output_basepath);
if ~exist(dir)
    mkdir(dir);
end
old = cd(dir);

%% Get the movie dimensions and number of frames:
disp('Getting input file info...');
info = h5info(input_path);
height = info.Datasets(1).Dataspace.Size(1);
width = info.Datasets(1).Dataspace.Size(2);
num_frames = info.Datasets(1).Dataspace.Size(3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
disp('... done');


%% Create output HDF5 objects and dataset:
n_chunks = ceil(num_frames/chunk_size);

disp('Creating output file...');
if ~multiple_outputs
    Outputs(1).name = h5create(output_basepath, '/mov', [height*width num_frames]);    
else
    for n = 1:n_chunks
        Outputs(n).name = [base '_' num2str(n) ext];
        disp(Outputs(n).name);
        h5create(Outputs(n).name, '/mov', [height*width chunk_size]);    
    end
end
disp('... done');


%% Reshape data and write to new HDF5 files:
for nn = 1:n_chunks
    
    if nn < n_chunks
        curr_chunk_size = chunk_size;
    else
        curr_chunk_size = mod(num_frames, chunk_size);
    end

    if multiple_outputs 
        start = [1 1];
        curr_output_file = Outputs(n).name;
    else
        start = [1 (nn-1)*chunk_size+1];
        curr_output_file = Outputs(1).name;
    end
    
    disp(['Vectorizing frames ' num2str((nn-1)*chunk_size+1) ' to ' num2str((nn-1)*chunk_size+curr_chunk_size) ' out of ' num2str(num_frames)]);
    
    disp('    ... reading from input file...');
    data = h5read(input_path, '/mov', [1 1 (nn-1)*chunk_size+1], [height width curr_chunk_size]); % Read data from source
    
    disp('    ... reshaping data...');
    data_reshaped = reshape(data, [height*width curr_chunk_size]); % Reshape
    
    disp('    ... writing to output file...');
    h5write(curr_output_file, '/mov', data_reshaped, start, [height*width curr_chunk_size]); % Write data to output
end

disp('Done.');

cd(old);