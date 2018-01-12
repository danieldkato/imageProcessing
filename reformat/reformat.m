function reformat(varargin)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-10


%% I. OVERVIEW:
% This function reformats an input movie in any of the following
% ways:

% 1) split it up into multiple smaller movies, each saved in its own output file 
% 2) vectorize each frame as a wh-element vector, where w is the movie width in pixels and h is the movie height in pixels
% 3) change the file type from .mat to .h5 or vice-versa


%% II. USAGE:

% In MATLAB, invoke this function with any of the following:

% reformat(input_path, output_path) - change file format by specifying different extensions for input_path and output_path
% reformat(input_path, output_path, chunk_size) - split input file up into multiple files
% reformat(input_path, output_path, [], vectorize) - vectorize frames
% reformat(input_path, output_path, chunk_size, vectorize) - vectorize data and split into multiple file 

% In addition to invoking this function from another MATLAB script or from
% the MATLAB command line, it is possible to invoke this function from the
% LINUX command line with the following:

% matlab -nosplash -nodesktop -r "reformat <input_path> <output_path> [<chunk_size>] [<multiple_outputs>]"


%% III. REQUIREMENTS:
% 1) MATLAB >= ???


%% IV. INPUTS:
% 1) input_path - path to an HDF5 or .mat file containing data to be
% reformatted. The imaging data to be converted should be saved in a w x h
% x f matrix, where w is the width of the movie in pixels, h is the height
% of the movie in pixels, and f is the number of frames in the movie. If
% the input file is an HDF5, this matrix should be saved in a dataset
% called '/mov'. If the input file is a .mat, this matrix should be saved
% in a variable called 'mov'.

% 2) output_path - path where the output file should be saved. This can be
% either an HDF5 or a .mat.

% 3) chunk_size (optional) - integer specifying the number of frames to
% read into memory at a time. Default value is 1000. Set to [] to save the
% output to one file.

% 4) vectorize (optional) - Boolean flag specifying whether to vectorize
% each frame. Default value is false.


%% V. OUTPUT:
% This function has no formal return, but saves to disk one or more output
% files. The file type of the outputs is determined by the extension in the
% output_path argument, which can be either .h5 or .mat. The number of the
% output files is determined by the chunk_size_argument, which specifies
% the number of frames to be written in each output file. If chunk_size is
% set to [], then all output will be written to one output file.

% If the output is saved to multiple output files, the name of each file
% will be a numbered version of the name specified in output_path. For
% example, if output_path  is specified as `/data/processed/Y.h5', then the
% output files will be named 'Y_1.h5', 'Y_2.h5', 'Y_3.h5', etc.

% If the output files are written as HDF5s, all data will be contained in a
% dataset called '/mov'. If the output files are written as .mat files, all
% data will be contained in a variable called mov.

% If vectorize is set to true, then each frame will be stored as a
% wh-element vector, where w is the movie width in pixels and h is the
% movie height in pixels.


%% TODO:
% 1) Throw warning if chunk_size is greater than the length of the movie
% 2) Throw warning if no output_dir can be extracted from the specified
% output path


%% Define parameters:

% Define some defaults:
chunk_size = 1000; % default chunk size
multiple_outputs = true;
vectorize = false;

% Get user-defined parameters:

input_path = varargin{1}; % input path
output_base = varargin{2}; % output path

% Determine whether to split output into multiple files, and if so, how many:
if nargin > 2
    chunk_size = varargin{3};
    if isempty(chunk_size)
        multiple_outputs = false;
    end
end

% Determine whether to vectorize each frame:
if nargin > 3
    vectorize = varargin{4};
end

% Get format of input and output files:
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

if ~multiple_outputs 
    chunk_size = num_frames;
end

disp('... done');


%% Process each chunk individually:

% cd to directory where output files will be saved:
old = cd(output_dir);

if multiple_outputs
    n_chunks = ceil(num_frames/chunk_size);
else
    n_chunks = 1; 
end

% For each chunk... 
for n = 1:n_chunks

    disp(['Processing chunk ' num2str(n) ' of ' num2str(n_chunks) ':']);
    
    % Define the output file name:
    if n_chunks == 1
        Outputs(n).name = output_base;        
    else
        Outputs(n).name = [output_name '_' num2str(n) output_ext];
    end
    
    % Define the chunk_size:
    if n < n_chunks || mod(num_frames, chunk_size) == 0
        curr_chunk_size = chunk_size;
    else
        curr_chunk_size = mod(num_frames, chunk_size);
    end        

    % Define the output dimensions:
    if ~vectorize
        output_dims = [height width curr_chunk_size];
    else
        output_dims = [height*width curr_chunk_size];
    end
    
    % Create the output object:
    disp('    Creating output file...');
    % ... if the user has requested an HDF5... 
    if strcmp(output_ext, '.h5')
        h5create(Outputs(n).name, '/mov', output_dims);   
    % ... if the user has requested a .mat... 
    elseif strcmp(output_ext, '.mat')
        O = matfile(Outputs(n).name,'Writable',true);
        O.mov = int16(zeros(output_dims));
    end
    
    % Define start and stop read frames:
    start_frame = (n-1)*chunk_size+1;
    end_frame = start_frame + curr_chunk_size - 1;
    
    % Read in the data:
    disp('    Reading input...');
    if strcmp(input_ext, '.h5')    
        data = h5read(input_path, '/mov', [1 1 start_frame], [height width curr_chunk_size]); 
    elseif strcmp(input_ext, '.mat')
        data = I.mov(1:height, 1:width, start_frame:end_frame);
    end    
    
    % Re-shape data if requested by user:
    if vectorize
        disp('    Vectorizing frames...');
        data = reshape(data, [height*width curr_chunk_size]);     
    end
    
    % Write data to output:
    disp(['    Writing output...']);
    if strcmp(output_ext, '.h5')
        if ~vectorize
            start_pos = [1 1 1];
        else
            start_pos = [1 1];
        end
        h5write(Outputs(n).name, '/mov', data, start_pos, output_dims); 
    elseif strcmp(output_ext, '.mat')
        if ~vectorize
            O.mov(:,:,1:curr_chunk_size) = data;
        else
            O.mov(:,1:curr_chunk_size) = data;
        end 
    end    
    
end

disp('Done.');


%% Define and write metadata:
disp('Writing metadata...')
Metadata.inputs(1).path = input_path;
for m = 1:length(Outputs)
       Metadata.outputs(m).path = fullfile(output_dir, Outputs(m).name);	
end

Metadata.parameters.chunk_size = chunk_size;
Metadata.parameters.vectorize = vectorize;

write_metadata(Metadata, 'reformat_metadata.json');

disp('... done.')

cd(old);
