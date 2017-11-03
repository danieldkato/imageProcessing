# ddk_scalpel0_wrapper.R

# DOCUMENTATION TABLE OF CONTENTS:
# I. OVERVIEW
# II. USAGE
# III. REQUIREMENTS
# IV. INPUTS
# V. OUTPUTS

# last updated DDK 2017-11-03


####################################################################################################
# I. OVERVIEW:

# This script is a wrapper for the first step ("step0") of the SCALPEL automated segmentation
# pipeline. This step is responsible for spatially and temporally smoothing the data and converting
# raw gray values to dF/F values. For more detail, see the SCALPEL documentation at https://rdrr.io/cran/scalpel.


####################################################################################################
# II. USAGE

# R < ddk_scalpel0_wrapper.R </path/to/params/file> --no-save


####################################################################################################
# III. REQUIREMENTS:

# 1) R.
# 2) The R package SCALPEL. For installation instructions, see https://rdrr.io/cran/scalpel.
# 3) The R package rjson. To install, call install.packages("rjson") from inside R. 


####################################################################################################
# IV. INPUTS:

# This command-line function takes a single argument, namely, the path to a parameters file (see USAGE
# above). This parameters file must include the following fields:

# 1) params$video_height - the height of the video, in pixels
# 2) params$output_folder - path to the directory where the output of this step should be saved
# 3) params$raw_data_folder - path to a .mat file containing the un-segmented input data. This
#    file must be named Y_1.mat and contain a single (m x n)-by-t matrix, where m is the video
#    height, n is the video width, and t is the number of frames in the video. For more detail,
#    see the documentation for scalpelStep0 at https://rdrr/io/cran/scalpel/man/scalpelStep0.html.

# For an example of how the JSON parameters file should be formatted, see
# https://github.com/danieldkato/imageProcessing/blob/master/image_segmentation/scalpel/step0/scalpel0_params.json


####################################################################################################
# V. OUTPUTS:

# This script is not (yet) a function and thus has no formal return. However, this wrapper  saves the
# following to secondary storage:

# 1) step0out.Rdata - R object containing information about the output of this processing step. This 
#    can subsequently be loaded into memory and passed as an argument to SCALPEL functions corresponding
#    to later processing steps. 
# 2) metadata.json - a JSON file containing step 0 metadata, including paths and SHA1 digests for inputs 
#    and outputs, as well as parameters.

# In addition to the files mentioned above, scalpelStep0 itself saves a number of files in a directory called
# Step0Data in the directory spcified in the `output_folder` parameter. 


####################################################################################################

# Load necessary libraries:
library("scalpel")
library("rjson")

# Get path to parameters file from command line:
args = commandArgs()

# Load parameters: 
jsondata <- fromJSON(file=args[2])
params <- jsondata$params

# Set step 0 parameters:
raw_data_folder = params$raw_data_folder
output_folder = params$output_folder
video_height = params$video_height

setwd(raw_data_folder)

# Perform scalpel step 0:
step0out = scalpelStep0(rawDataFolder=raw_data_folder,outputFolder=output_folder,videoHeight=video_height,fileType="matlab")

# Save step0out object in most recent Step0 directory:
step0_indices = unlist(grepl("Step0",dir())) # get indices of all directories containing the string "Step0"
step0_dirs = dir()[step0_indices] # get the names of all directories containing the string "Step0"
mtimes = file.info(step0_dirs)$mtime # get the time that each directory was last modified
most_recent_step0 = step0_dirs[which.max(mtimes)] # find the most recently modified directory containing string "Step0"
output_path = paste(most_recent_step0,"step0out.Rdata",sep="/")
save(step0out,file=output_path)
start = Sys.time()
date = format(start,format="%Y-%m-%d")
time = strftime(start,format="%H:%M:%S")

# Save step 0 metadata:
step0_metadata = list()

input_path = paste(raw_data_folder,"Y_1.mat",sep="/")
sys_out = system(paste("sha1sum",input_path,sep=" "),intern=TRUE)
input_sha1 = substr(sys_out,1,40)
input_1 = list(path=input_path,sha1=input_sha1)
step0_metadata$inputs <- list(input_1)

step0_metadata$params$raw_data_folder = raw_data_folder
step0_metadata$params$output_folder = output_folder
step0_metadata$params$video_height = video_height

output_full_path = paste(raw_data_folder,most_recent_step0,"step0out.Rdata",sep="/")
sys_out2 = system(paste("sha1sum",output_full_path,sep=" "),intern=TRUE)
output_sha1 = substr(sys_out2,1,40)
output_1 = list(path=output_full_path,sha1=output_sha1)
step0_metadata$outputs <-list(output_1)

step0_metadata$date = date
step0_metadata$time = time
m = toJSON(step0_metadata)
write(m,file=paste(most_recent_step0,"metadata.json",sep="/"))
