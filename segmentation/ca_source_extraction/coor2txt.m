%Last updated DDK 2016-09-26

% OVERVIEW:
% This script is used for exporting automatically-drawn ROIs from Eftykios
% Pnevmatikakis' ca_source_extraction software to a .txt format that can
% then be read into Fiji for manual refinement. 

% REQUIREMENTS:
% This script is meant to be run on the output of Eftykios Pnevmatikakis'
% ca_source_extraction software, available at:

% https://github.com/epnev/ca_source_extraction/blob/master/CSHL16/demo_somatic.m(commit ff3c6b6) 


% INSTRUCTIONS:
% Before running this script, run ca_source_extraction on the multi-page
% TIFF to be segmented. Save the output cell array Coor as a .mat file. 

% Run this script from the directory in which Coor.mat is saved. For each
% ROI, the script will create a two-column .txt file, where each row is a
% point defining the contours of the ROI, and the first and second columns
% represent X- and Y- coordinates of each point, respectively. The script
% will then save these .txt files into a directory called Coor2txt. If the
% directory already exists, the script will ask the user whether to
% overwrite the existing directory, create a new directory, or cancel the
% operation. 

% These .txt files can then be loaded one at a time into Fiji using File>Import>XY
% Coordinates... . Alternatively, they can be imported as a batch using the
% Batch Import ROIs macro, available at: 

% https://github.com/danieldkato/cse2Fiji/blob/master/batchImportROIs.txt



%%
struct = load('Coor.mat');
Coor = struct.Coor;
outputBaseName = 'Coor2txt';
exportComplete = 0;

while exportComplete == 0;
    A = exist(outputBaseName, 'dir'); %Check if the directory Coor2txt already exists
    if A == 7 %If the directory Coor2txt already exists, ask the user for appropriate action
        button = questdlg('The directory Coor2txt already exists. Select whether to overwrite existing directory, create a new directory, or cancel the operation.', 'Confirm folder replace', 'Overwrite', 'Create new directory', 'Cancel', 'Cancel');
        switch button
            
            case 'Overwrite'
                %If the user indicates that it's ok to overwrite the
                %existing directory, then just name the output directory
                %the base name.
                rmdir(outputBaseName, 's');
                outputDirectory = outputBaseName;
            
            case 'Create new directory'
                %If the user says to create a new directory, then add a
                %numeric suffix to the output directory base name. The
                %numeric suffix must be higher than the suffix of any
                %existing directory that has the output directory base name
                %(e.g., if Coor2txt_2 already exists, then name the output
                %directory Coor2txt_3).
                maxSuffix = 0;
                listing = dir(pwd); %Get a listing of all of the contents of the current working directory
                for i = 1:length(listing) 
                    name = listing(i).name;
                    if listing(i).isdir && length(name)>=length(outputBaseName) && strcmp(name(1:length(outputBaseName)), outputBaseName) == 1 %If the item is a directory and its name matches the output directory base name...
                        lastChar = name(length(name)); %... then get the last character of the name...
                        if ~isempty(str2num(lastChar)) && str2num(lastChar)>maxSuffix %... and if the last character is a number and it's higher than the highest suffix recorded so far...
                            maxSuffix = str2num(lastChar); %... then record it in maxSuffix
                        end
                    end
                end
                newSuffix = maxSuffix + 1;
                outputDirectory = strcat(outputBaseName, '_', num2str(newSuffix));
            
            case 'Cancel'
                %If the user says to cancel the operation, then mark the
                %export as complete and break the the while loop so that
                %the script completes execution with no further action
                %being taken.
                exportComplete = 1;
                break
        end
    end
    
    mkdir(outputDirectory);
    cd(outputDirectory);
    for i = 1:length(Coor)
        id = fopen(strcat(['ROI_', num2str(i), '_coords.txt']), 'wt');
        currentROI = Coor{i};
        for j = 1:length(currentROI)
            fprintf(id, [num2str(currentROI(1,j)), ' ']);
            fprintf(id, [num2str(currentROI(2,j)), '\n']);        
        end
        fclose(id);
    end
    cd('../')
    exportComplete = 1; 
end