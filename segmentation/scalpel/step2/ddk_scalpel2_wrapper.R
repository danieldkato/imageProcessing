# ddk_scalpel2_wrapper.R

# DOCUMENTATION TABLE OF CONTENTS:
# I. OVERVIEW
# II. USAGE
# III. REQUIREMENTS
# IV. INPUTS
# V. OUTPUTS

# last updated DDK 2017-11-02


####################################################################################################
# I. OVERVIEW:
# This script is a wrapper for the third step ("step 2") of the SCALPEL automated segmentation
# pipeline. This step is responsible for merging the elements of the preliminary spatial components
# dictionary produced by step 1. For more detail, see the SCALPEL documentation at https://rdrr.io/cran/scalpel.


####################################################################################################
# II. USAGE:

# R < ddk_scalpel2_wrapper.R <path/to/params/file> --no-save


####################################################################################################
# III. REQUIREMENTS:

# 1) R.
# 2) The R package SCALPEL. For installation instructions, see https://rdrr.io/cran/scalpel.
# 3) The R package rjson. To install, call install.packages("rjson") from inside R. 
# 4) The R package R.matlab. To install, call install.packages("R.matlab") from inside R. 


####################################################################################################
# IV. INPUTS:
# This command-line function takes one argument, namely, a path to a parameters file (see USAGE above).
# This parameters file must include the following fields:

# 1) params$omega - value in [0,1] indicating how to weight spatial vs. temporal information in the
#    dissimilarity metric used for clustering.
# 2) params$cutoff - value in [0,1] indicating where to cut the dendrogram that results from 
#    hierarchical clustering of the preliminary dictionary elements.

# For more detail on these parameters, see https://rdrr.io/cran/man/scalpelStep2.html and
# https://arxiv.org/abs/1703.06946.

# For an example of how the JSON parameters file should be formatted, see
# https://github.com/danieldkato/imageProcessing/blob/master/image_segmentation/scalpel/step2/scalpel2_params.json


####################################################################################################
# V. OUTPUTS:
# This script is not (yet) a function and thus has no formal return. However, this wrapper  saves the
# following to secondary storage:

# 1) step2out.Rdata - R object containing information about the output of this processing step. This 
#    can subsequently be loaded into memory and passed as an argument to SCALPEL functions corresponding
#    to later processing steps.
# 2) A.mat - a .mat file containg a (w x h)-by-n binary matrix, where w is the video width, h is the video
#    height, and n is the number of refined spatial components identified by step 2. Each column represents
#    an individual spatial component, and each row within that column indicates whether the corresponding
#    pixel is part of the spatial component.
# 3) metadata.json - a JSON file containing step 2 metadata, including paths and SHA1 digests for inputs 
#    and outputs, as well as parameters.

# In addition to the files mentioned above, scalpelStep2 itself saves a number of files in a directory called
# Step2_omega_<w>_cuoff_<c>, where <w> is the omega value and <c> is the cutoff value for the current analysis. 


####################################################################################################
library("scalpel")
library("R.matlab")
library("rjson")

args = commandArgs()
print(args)

jsondata <- fromJSON(file=args[2])
inputPath = jsondata$inputs[[1]]$path
cutoff = jsondata$params$cutoff
omega = jsondata$params$omega

load(inputPath)
step1_dir = dirname(inputPath)
setwd(step1_dir)

print(cutoff)
print(class(cutoff))
print(omega)
print(class(omega))
str(step1out)
print(class(step1out))

# Perform SCALPEL step 2:
step2out = scalpelStep2(step1Output=step1out,cutoff=cutoff,omega=omega)
# this will automatically save some output into a dedicated sub-directory of the directory for the step 1 output; I also want to save the step2out object itself, along with some metadata, into this dedicated sub-directory. This should be the most recent directory whose name contains the string "Step2"

# Find the most recent "Step2" directory:
step2_indices = unlist(grepl("Step2",dir())) # get indices of all directories whose names contain the string "Step2"
step2_dirs = dir()[step2_indices] # get the names of all directories containing the string "String2"
mtimes = file.info(step2_dirs)$mtime # get the time each of these folders was last modified
now = Sys.time() # get current time
time_diffs = now - mtimes # find out how long ago each was last modified
most_recent_ind = which.min(time_diffs)
newest_step2_dir = step2_dirs[most_recent_ind]
setwd(newest_step2_dir) # cd to the msot recent step 2 directory to save some additional outputs there

# Save step2out as an .Rdata object:
s2obj_full_path = paste(step1_dir,newest_step2_dir,"step2out.Rdata",sep="/")
save(step2out,file=s2obj_full_path)

# Save step2out$A as a .mat (will be useful for converting to txt files of coordinates for importing into ImageJ):
mat_full_path = paste(step1_dir,newest_step2_dir,"A.mat",sep="/")
writeMat(mat_full_path,A=step2out$A)

# TODO: Convert A.mat to .txt?

now = Sys.time()
date = format(now,format="%Y-%m-%d") 
time = strftime(now,format="%H:%M:%S")

# Write metadata:
step2_metadata = list()

sys_out = system(paste("sha1sum",inputPath,sep=" "),intern=TRUE)
input_sha1 = substr(sys_out,1,40)
inpt1 <- list(path=inputPath,sha1=input_sha1)
step2_metadata$inputs <- list(inpt1)

step2_metadata$params <- jsondata$params

sys_out2 = system(paste("sha1sum",s2obj_full_path,sep=" "),intern=TRUE)
output1_sha1 = substr(sys_out2,1,40)
output1 <- list(path=s2obj_full_path,sha1=output1_sha1)
sys_out3 = system(paste("sha1sum",s2obj_full_path,sep=" "),intern=TRUE)
output2_sha1 = substr(sys_out3,1,40)
output2 <- list(path=mat_full_path,sha1=output2_sha1)
step2_metadata$outputs <- list(output1, output2)

step2_metadata$date <- date
step2_metadata$time <- time

m <- toJSON(step2_metadata)
write(m,file=paste(step1_dir,newest_step2_dir,"metadata.json",sep="/"))

# Return to previous working directory:
setwd(step1_dir)

