function dir_name = generate_mc_dir_name()
%%  
% DOCUMENTATION TABLE OF CONTENTS:
% I. OVERVIEW
% II. REQUIREMENTS
% III. INPUTS
% IV. OUTPUTS

% last updated DDK 2017-11-29


%% I. OVERVIEW:
% This function generates a directory name for the output of a motion
% correction script or function. The output directory name has the format
% 'motion_correction<num>', where <num> is some number.

% The function checks whether there already exist directories with names of
% this format in the current directory, and if so, increments the number by
% 1, pads it appropriately, and appends it to 'motion_correction.'


%% II. REQUIREMENTS:
% None.

%% III. INPUTS:
% None.

%% IV. OUTPUTS:
% 1) dir_name - char array of the form 'motion_correction<num>', where
% <num> is some number.


%%
ls = dir();
mc_basename = 'motion_correction';

% Get the names of all directories of the form 'motion_correction<num>':
mc_directory_indices = arrayfun(@(x) length(x.name) > length(mc_basename) && strcmp(x.name(1:length(mc_basename)), mc_basename), ls); 
mc_directories = ls(mc_directory_indices);

% Get numbers at end of directory names:
nums = arrayfun(@(x) x.name(length(mc_basename)+1:end), mc_directories, 'UniformOutput', false);

max = 0;

for n = 1:length(nums)
    val = str2double(nums{n});
    if val > max
        max = val;
    end
end

num_suffix = pad(num2str(max+1), 3, 'left', '0');
dir_name = [mc_basename num_suffix];

end