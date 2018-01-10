function fmt4scalpel_h5(varargin)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-06


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
% 1) input_path - path to a file containing data to be vectorized. The
% imaging data to be converted should be saved in a w x h x f dataset
% called '/mov', where w is the width of the movie in pixels, h is the
% height of the movie in pixels, and f is the number of frames in the
% movie.

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
% 2) Throw warning if no output_dir can be extracted from the specified
% output path


%% Define parameters:
input_path = varargin{1};
output_base = varargin{2};
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


%% Get format of input and output files:
[input_dir, input_base, input_ext] = fileparts(input_path);
[output_dir, output_name, output_ext] = fileparts(output_base);

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end


%% Get the movie dimensions and number of frames:
disp('Getting input file info...');

if strcmp(input_ext, '.h5')
    info = h5info(input_path);
    height = info.Datasets(1).Dataspace.Size(1);
    width = info.Datasets(1).Dataspace.Size(2);
    num_frames = info.Datasets(1).Dataspace.Size(3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
elseif strcmp(input_ext, '.mat')
    I = matfile(input_path);
    height = size(I.mov,1);
    width = size(I.mov,2);
    num_frames = size(I.mov,3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
end

disp('... done');



%% Create output HDF5 objects and dataset:

% cd to directory where output files will be saved:
old = cd(output_dir);

disp(size(num_frames));
disp(size(chunk_size));
n_chunks = ceil(num_frames/chunk_size);

disp('Creating output files...');
% Setup output file if the user has requested only one output file... 
if ~multiple_outputs
    Outputs(1).name = output_base;
    % ... if the user has requested an HDF5... 
    if strcmp(output_ext, '.h5')    
        h5create(Outputs(1).name, '/mov', [height*width num_frames]);    
    % ... if the user has requested a .mat... 
    elseif strcmp(output_ext, '.mat')
        Outputs(1).matfile = matfile(Outputs(1).name,'Writable',true);
        Outputs(1).matfile.mov = int16(zeros([height*width num_frames]));
    end
    
% Setup output files if the user has requested multiple output files... 
else
    for n = 1:n_chunks

        % Define the current chunk size (the last chunk may be less that chunk_size)
        if n < n_chunks
            curr_chunk_size = chunk_size;
        else
            curr_chunk_size = mod(num_frames, chunk_size);
        end        
        
        Outputs(n).name = [output_name '_' num2str(n) output_ext];
        disp(Outputs(n).name);
        % ... if the user has requested an HDF5... 
        if strcmp(output_ext, '.h5')
            h5create(Outputs(n).name, '/mov', [height*width chunk_size]);   
        % ... if the user has requested a .mat... 
        elseif strcmp(output_ext, '.mat')
            Outputs(n).matfile = matfile(Outputs(n).name,'Writable',true);
            Outputs(n).matfile.mov = int16(zeros([height*width chunk_size]));
        end
    end
end
disp('... done');


%% Reshape data and write to new HDF5 files:
for nn = 1:n_chunks
    
    start_read_frame = (nn-1)*chunk_size+1; 
    
    % Define the current chunk size (the last chunk may be less that chunk_size)
    if nn < n_chunks || mod(num_frames, chunk_size) == 0
        curr_chunk_size = chunk_size;
    else
        curr_chunk_size = mod(num_frames, chunk_size);
    end

    % Determine the current output file (depends on whether the user has requested one or multiple output files, and on requested output type):
    if multiple_outputs 
        start_write_frame = 1;
        if strcmp(output_ext, '.h5')
           curr_output_file = Outputs(nn).name;
        elseif strcmp(output_ext, '.mat')
           curr_output_file = Outputs(nn).matfile;
        end
    else
        start_write_frame = start_read_frame;
        if strcmp(output_ext, '.h5')
          curr_output_file = Outputs(1).name;
        elseif strcmp(output_ext, '.mat')
          curr_output_file = Outputs(1).matfile;
        end
    end
    
    disp(['Vectorizing frames ' num2str((nn-1)*chunk_size+1) ' to ' num2str((nn-1)*chunk_size+curr_chunk_size) ' out of ' num2str(num_frames)]);
    
    % Read data from source:
    disp('    ... reading from input file...');
    if strcmp(input_ext, '.h5')    
        data = h5read(input_path, '/mov', [1 1 start_read_frame], [height width curr_chunk_size]); 
    elseif strcmp(input_ext, '.mat')
        data = Outputs(nn).mat.mov(1:height, 1:width, start_read_frame:start_read_frame + curr_chunk_size - 1);
    end
    
    % Reshape data:
    disp('    ... reshaping data...');
    data_reshaped = reshape(data, [height*width curr_chunk_size]); 
    
    % Write data to output:
    disp(['    ... writing to  output file ' num2str(nn) '...']);
    if strcmp(output_ext, '.h5')
        h5write(curr_output_file, '/mov', data_reshaped, [1 start_write_frame], [height*width curr_chunk_size]); 
    elseif strcmp(output_ext, '.mat')
        curr_output_file.mov(:,start_write_frame:start_write_frame + curr_chunk_size -1) = data_reshaped; 
    end
end

disp('Done.');

cd(old);
