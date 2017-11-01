function tiff_2_mat(input_path, output_path)

% Get data dimensions:
info = imfinfo(input_path);

% Initialize output object:
disp('Initializing memory-mapped output file...');
Y = matfile(output_path, 'Writable', true);
Y.Y = uint16(NaN(info(1).Width,info(2).Height,length(info)));
disp('... done.');

% Copy data from input TIFF to memory-mapped output file:
disp('Copying image data...');
for f = 1:length(info) 
    Y.Y(:,:,f) = imread(input_path,'Index',f);
end
disp('... done.');

end