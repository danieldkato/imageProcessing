function A = extract_activity(input_directory, output_path)
% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II.REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% Last updated DDK 2018-01-11


%% I. OVERVIEW:
% This function reads in raw fluorescence data from multiple .csv files
% containing ImageJ multi measure data, consolidates it into a single raw
% activity matrix A, and, if requested by the user, saves A to secondary
% storage. The dimensions of A are n x f, where n is the number of neurons
% in the movie, and f is the number of frames in the movie.


%% II. REQUIREMENTS:
% 1) write_metadata.m, available at https://github.com/danieldkato/utilities/blob/master/metadata/MATLAB/write_metadata.m


%% III. INPUTS:
% 1) input_directory - path to a directory containing one or more .csv
% files containing output from ImageJ's multi measure tool. Each .csv
% should have the following format:

%          area(ROI_1) | mean(ROI_1) | min(ROI_1) | max(ROI_1) | area(ROI_2) | mean(ROI2) | min(ROI_2) | max(ROI_2) | ...
% frame 1
% frame 2
% frame 3
% ...

% If there are multiple .csv files, they should be numbered to reflect
% their order. E.g., 'Y_1.csv' should contain the first part of the data,
% 'Y_2.csv' should contain the second part of the data, etc.


% 2) output_path (optional) - path to file where consolidated raw
% fluorescence data should be saved. Can be specified as a .mat or .h5. If
% not specified, then A will not be saved to secondary storage.


%% IV. OUTPUTS: 
% 1) A - n x f raw fluorescence matrix, where n is the number of neurons in
% the movie and f is the number of frames in the movie.

% In addition to returning A formally, this functiona also saves A to
% secondary storage if the output_path argument is specified.


%% Harcode some things:
cols_per_frame = 4; % in the MultiMeasure output from ImageJ, each ROI is represented by 4 columns (area, mean, min, and max); unlikely to change


%% Find any .csv files in the specified input directory: 
cd(input_directory);
ls = dir();
names = arrayfun(@(x) x.name, ls, 'UniformOutput', false);
is_csv = cell2mat(cellfun(@(c) ~isempty(regexp(c, '.csv', 'ONCE')), names, 'UniformOutput', false));

% Throw a warning if there are no .csv files in the directory:
if ~find(is_csv)
    error('No .csv files detected in specified input directory. Please ensure that correct input directory has been specified and multi-meausre files have been formatted properly.');
end 

% Get a struct array of the .csv files:
csvs = ls(is_csv);


%% Make sure that .csv files are ordered appropriately:

% Try to retrieve the numbers contained in the filenames of any .csv files:
filenum_labels = arrayfun(@(a) regexp(a.name, '[0-9]*', 'match'), csvs, 'UniformOutput', false);

% If there are any un-numbered .csv files, throw an error:
not_numbered = cellfun(@(c) isempty(c), filenum_labels);
if find(not_numbered)
    error('Non-numbered .csv files detected. Please ensure that all .csv files are numbered correctly.');
end

% Sort the .csv files by the number in the title:
filenums = cellfun(@(c) str2double(c{:}), filenum_labels);
[fnums_sorted, I] = sort(filenums); 
csvs = csvs(I);

% Display the .csv file names in order to confirm that they are ordered
% correctly:
for c = 1:length(csvs)
    disp(csvs(c).name);
end


%% Load the data from each .csv into activity matrix A: 
A = [];
for c = 1:length(csvs)
    
    disp(['Reading in data from file ' num2str(c) ' out of ' num2str(length(csvs)) '...']);
    
    % Load the data from the .csv:
    num = xlsread(csvs(c).name);

    % Get the dimensions:
    ncols = size(num, 2);

    % Define the columns to read; we only want the mean columns; see INPUTS
    % for how input data is formatted:
    indices = (1:cols_per_frame:ncols-cols_per_frame); % create a vector with the appropriate number of elements
    roi_start_indices = indices + 1; % offset by 1 because first column is just the "slice" (i.e. frame) number
    read_indices = roi_start_indices + 1; % offset by 1 more because first column of each ROI is area
    num_rois = length(read_indices); 

    % Confirm that all files have the same number of ROIs; throw an error if not:
    if c == 1
        first_file_rois = num_rois;
    elseif num_rois ~= first_file_rois
        error([csvs(c).name ' and ' csvs(c-1).name ' contain different numbers of ROIs. Please make sure that these files come from the same dataset.']);
    end

    % Extract the data from the relevant columns:
    A_temp = num(read_indices)';

    % Append to overall data matrix:
    A = [A A_temp];
end
   

%% Write activity matrix A to secondary storage if requested by user:

if nargin > 1
    
    disp('Writing data to secondary storage...');
    
    % Get the specified file type:
    [output_dir, output_name, output_ext] = fileparts(output_path);

    % Create the output directory if it doesn't exist:
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    
    old = cd(output_dir);
    
    % Save A as appropriate file type:
    if strcmp(output_ext, '.mat')
        save(output_path, 'A');
    elseif strcmp(output_ext, '.h5')
        %h5_name = [output_name output_ext];
        h5create([output_name output_ext], '/A', size(A));
        h5write([output_name output_ext], '/A', A, [1 1], size(A));
    end
    
    % Write metadata:
    for m = 1:length(csvs)
        Metadata.inputs(m).path = fullfile(input_directory, csvs(m).name);
    end    
    Metadata.outputs(1).path = output_path;
    write_metadata(Metadata, 'activity_extraction_metadata.json');
    
    cd(old);
end