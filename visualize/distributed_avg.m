function avg_img = distributed_avg(path, bin_size, output_path)
% DOCUMENTATION TABLE OF CONTENTS:

% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-18


%% I. OVERVIEW: 
% This function creates an average image of a video by sampling
% evenly-spaced frames throughout the video and optionally saves it to
% secondary storage. 

% This is convenient when a video is too large to load
% into memory (e.g., a motion-corrected HDF5 generated by NoRMCorre).


%% II. REQUIREMENTS:
% 1) saveastiff, avaialble at https://github.com/flatironinstitute/NoRMCorre


%% III. INPUTS: 
% 1) path - path to a video file to be averaged. If the file is an HDF5,
%    this function will use the first dataset in the file.

% 2) bin_width - number of frames between sample frames 

% 3) output_path (optional) - path where the average image should be saved
%    to secondary storage. Currently, this can be either an HDF5 for a
%    TIFF.


%% IV. OUTPUTS:
% 1) avg_img - h x w image matrix containing the averaged image data, where
%    h is the height of the original video in pixels and w is the width of the
%    original video in pixels. 

% In addition to formally returining avg_img, this function saves avg_img
% to secondary storage if the output_path argument is set. 


%%
[indir, inname, inext] = fileparts(path);

if strcmp(inext, '.h5')
    info = h5info(path);
    
    if length(info.Datasets) > 1
        warning('More than one dataset detected in input HDF5 file; using first dataset.');
    end
    
    dims = info.Datasets(1).Dataspace.Size;
    if length(dims) ~= 3
        error('Number of dimensions in requested dataset not equal to 3. Please check that path to correct input file has been specified.');
    end
    height = dims(1);
    width = dims(2);
    num_frames = dims(3);
end

if bin_size > num_frames
    error('Requested bin size greater than number of frames in video; please specify a bin size smaller than number of frames in video.');
end

sample_frames = (1:bin_size:num_frames);
num_sample_frames = length(sample_frames);
disp(['Sampling ' num2str(num_sample_frames) ' frames spaced ' num2str(bin_size) ' frames apart...']);

M = int16(zeros(height, width, num_sample_frames));
for i = 1:num_sample_frames
    
    frame_num = sample_frames(i);
    
    % read data:
    disp(['Reading frame ' num2str(i) ' out of ' num2str(num_sample_frames) '...']);
    if strcmp(inext, '.h5')
        img_dat = h5read(path, ['/' info.Datasets(1).Name], [1 1  frame_num], [height width 1]);
        disp(class(img_dat));
        M(:,:,i) = img_dat;
    end
end

disp('Averaging sample frames...');
disp(class(M));
avg_img = int16(mean(M,3));
disp(class(avg_img));
disp(size(avg_img));

if ~isempty(output_path)
    
    disp('Saving mean image to secondary storage...');
    
    [outdir, outname, outext] = fileparts(output_path);
    
    if strcmp(outext, '.h5')
        h5create(output_path, '/avg', size(avg_img));
        h5write(output_path, '/avg', avg_img);
    elseif strcmp(outext, '.tif')
        saveastiff(avg_img, output_path);
    end
    
end

