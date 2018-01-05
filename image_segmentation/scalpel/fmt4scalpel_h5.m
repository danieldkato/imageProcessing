function fmt4scalpel_h5(varargin)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-04


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

% In addition to invoking this function from another MATLAB script or from
% the MATLAB command line, it is possible to invoke this function from the
% LINUX command line with the following:

% matlab -nosplash -nodesktop -r "fmt4scalpel_h5 <input_path> <output_path> [<chunk_size>]"


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


%% V. OUTPUT:
% This function has no formal return, but saves to disk an HDF5 file
% containing vectorized imaging data. These data are saved in a w*h x f
% dataset called '/mov', where w is the width of the movie in pixels, h is
% the height of the movie in pixels, and f is the number of frames in the
% movie.


%% TODO:
% 1) Would be nice to have some way to also save output as a TIFF for
% quick visual inspection


%%
% Define parameters:
input_path = varargin{1};
output_path = varargin{2};
if nargin > 2
    chunk_size = varargin{3};
else
    chunk_size = 1000;
end

% Get the movie dimensions and number of frames:
disp('Getting input file info...');
info = h5info(input_path);
height = info.Datasets(1).Dataspace.Size(1);
width = info.Datasets(1).Dataspace.Size(2);
num_frames = info.Datasets(1).Dataspace.Size(3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
disp('... done');

% Create output HDF5 object and dataset:
disp('Creating output file...');
h5create(output_path, '/mov', [height*width num_frames]);
disp('... done');

% Reshape data and write to new HDF5:
n_chunks = ceil(num_frames/chunk_size);
for n = 1:n_chunks
    
    if n < n_chunks
        curr_chunk_size = chunk_size;
    else
        curr_chunk_size = mod(num_frames, chunk_size);
    end
        
    disp(['Writing frames ' num2str((n-1)*chunk_size+1) ' to ' num2str((n-1)*chunk_size+curr_chunk_size) ' out of ' num2str(num_frames)]);
    
    data = h5read(input_path, '/mov', [1 1 (n-1)*chunk_size+1], [height width curr_chunk_size]);
    data_reshaped = reshape(data, [height*width curr_chunk_size]);
    h5write(output_path, '/mov', data_reshaped, [1 (n-1)*chunk_size+1], [height*width curr_chunk_size]);
end