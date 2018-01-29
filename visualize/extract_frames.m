function data = extract_frames(input_path, start_frame, end_frame, output_path)
%% DOCUMENTATION TABLE OF CONTENTS:      
% I. OVERVIEW
% II. USAGE
% III. REQUIREMENTS
% IV. INPUTS
% IV. OUTPUTS

% last updated DDK 2018-01-29


%% I. OVERVIEW:
% This function extracts a subset of frames from an imaging dataset and
% saves them to secondary storage in a format that can easily be
% visualized, e.g., using ImageJ. This is especially useful for visualizing
% subsets of datasets that would take to long to be conveniently saved as
% TIFFs in their entirety (like the output of normcorre) or loaded into
% memory for visualization.


%% II. USAGE:

% In MATLAB, invoke this function with any of the following:

% reformat(input_path, output_path, start_frame, end_frame)

% In addition to invoking this function from another MATLAB script or from
% the MATLAB command line, it is possible to invoke this function from the
% LINUX command line with the following:

% matlab -nosplash -nodesktop -r "extract_frames <input_path> <output_path> start_frame end_frame"


%% III. REQUIREMENTS:
% 1) MATLAB >= ???


%% IV. INPUTS:
% 1) input_path - path to an input file from which frames should be
% extracted. Can be formatted as an HDF5. The input file is an HDF5, this
% matrix should be saved in a dataset called '/mov'. If the input file is a
% .mat, this matrix should be saved in a variable called 'mov'.

% 2) start_frame - integer specifying the first frame in the range of
% frames to be extracted.

% 3) end_frame - integer specifying the last frame in the range of frames
% to be extracted.

% 4) output_path (optional) - path where the output file should be saved. This can be
% specified either an HDF5 or a TIFF using the '.h5' or '.tif' extensions,
% respectively.


%% V. OUTPUT:
% This function has no formal return, but saves to secondary storage an
% output file containing the extracted frames. This can be specified either
% as an HDF5 or a TIFF using the '.h5' or '.tif' extensions, respectively.

% If the output files are written as HDF5s, all data will be contained in a
% dataset called '/mov'.


%% TODO:
% 1) Add support for reading from .mat


%%
n_requested_frames = end_frame-start_frame+1;

% Get input extension:
[dirname_in, fname_in, input_ext] = fileparts(input_path);

% Get basic file info:
if strcmp(input_ext, '.h5') || strcmp(input_ext, '.hdf5')
    info = h5info(input_path);
    height = info.Datasets(1).Dataspace.Size(1);
    width = info.Datasets(1).Dataspace.Size(2);
    n_frames_total = info.Datasets(1).Dataspace.Size(3); % assuming there's only one dataset; TODO; include code for when there are multiple datasets
elseif strcmp(input_ext, '.tif')
    info = imfinfo(input_path);
    width = info(1).Width;
    height = info(1).Height;
    n_frames_total = numel(info);
end

subset_dims = [height width n_requested_frames]; % define size of extracted image data matrix 
data = nan(subset_dims); % initialize extracted image data matrix

% Check that the requested number of frames doesn't exceed the length of the movie:
if end_frame > n_frames_total
    error('Requested frames numbers exceed the number of frames in movie.');
end

% Read in the appropriate frames:
if strcmp(input_ext, '.h5') || strcmp(input_ext, '.hdf5')
    data = h5read(input_path, '/mov', [1 1 start_frame], subset_dims);
elseif strcmp(input_ext, '.tif')
    tiff_obj = Tiff(input_path);
    for j = start_frame:end_frame
        tiff_obj.setDirectory(j);
        data(:,:,j) = double(tiff_obj.read());
    end
end


%% If requested by user, save output frames:
if nargin > 3
    [dirname_out, fname_out, output_ext] = fileparts(output_path);
    if strcmp(output_ext, '.h5')
        h5create(output_path, '/mov', subset_dims);
        h5write(output_path, '/mov', data, [1 1 1], subset_dims);
    elseif strcmp(output_ext, '.tif')
        data = int16(data);
        saveastiff(data, output_path);
    end
end
