function Z_processed = process_zstack(path)
% DOCUMENTATION TOC:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-03


%% I. OVERVIEW
% This function averages together frames from the same anatomical slice in
% a z-stack saved as a multi-page TIFF then returns the averaged slice
% images as a 3D matrix.


%% II. REQUIREMENTS
% 1) MATLAB v >= ???


%% IV. INPUTS
% 1) path - path to a z-stack saves as a multi-page TIFF with s x f frames,
% where s is the number of slices and f is the number of frames per slice.


%% V. OUTPUTS
% 1) Z_processed - h x w x s matrix of averaged z-stack image data, where h
% is the image height, w is the image width, and s is the number of
% anatomical slices spanned by the z-stack. 


%% Initialize Z_processed:
Z_processed = [];


%% Get z-stack metadata:
disp('Getting z-stack metadata...')
z_metadata = loadjson('metadata.json'); 
% TODO: throw warning if metadata doesn't exist
% TODO: validate metadata
zinfo = imfinfo(path);
num_zframes = length(zinfo);
frames_per_slice = z_metadata.framesPerSlice;  
disp('... done.')
  

%% Load the z-stack image data:
disp('Loading z-stack image data...');
z_tiff = Tiff(path);
Z = NaN(zinfo(1).Height, zinfo(2).Width, num_zframes);
for slice = 1:num_zframes
    z_tiff.setDirectory(slice);
    Z(:,:,slice) = z_tiff.read(); 
end
disp('... done');
    

%% Average together frames from each slice:
('Averaging frames from the same slice...');

if frames_per_slice > 1
    if mod(num_zframes, frames_per_slice) == 0 % confirm that the recorded number of frames per slice divides evenly into the number of frames in the z-stack

        num_slices = num_zframes/frames_per_slice;
        Z_processed = NaN(zinfo(1).Height, zinfo(2).Width, num_slices);

        ssi = frames_per_slice * ones(1, num_slices);
        slice_start_indices = cumsum(ssi) - slices_per_frame + 1;

        for i = 1:length(slice_start_indices)
            Z_processed(:,:,i) = mean(Z(:,:,slice_start_indices(i):start_slice_indices(i)+frames_per_slice-1),3);
        end

        disp('... done');

    else
        warn('Frames per slice must divide evenly into number of frames in z-stack; check that frames per slice recorded in z-stack metadata file is correct; skipping registration of average images to z-stack');
    end
else 
    Z_processed = Z; 
end    
