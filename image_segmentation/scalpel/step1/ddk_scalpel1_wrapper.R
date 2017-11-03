# ddk_scalpel1_wrapper.R

# DOCUMENTATION TABLE OF CONTENTS:
# I. OVERVIEWi
# II. USAGE
# III. REQUIREMENTS
# IV. INPUTS
# V. OUTPUTS

# last updated DDK 2017-11-03


####################################################################################################
# I. OVERVIEW:

# This script is a wrapper for the second step ("step 1") of the SCALPEL automated segmentation
# pipeline. This step is responsible for creating a preliminary dictionary of spatial components
# by thresholding each dF/F transformed frame of the given input movie. For more detail, see the
# SCALPEL documentation at https://rdrr.io/cran/scalpel.


####################################################################################################
# II. USAGE:

# R < ddk_scalpel1_wrapper.R </path/to/params/file> --no-save


####################################################################################################
# III. REQUIREMENTS:

# 1) R.
# 2) The R package SCALPEL. For installation instructions, see https://rdrr.io/cran/scalpel.
# 3) The R package rjson. To install, call install.packages("rjson") from inside R. 


####################################################################################################
# IV. INPUTS:

# This command line function takes a single argument, namely, the path to a JSON parameters file (see 
# USAGE above). This parameters file must minimally include the following fields:

# 1) params$min_size - the minimum allowable area, in pixels, of a preliminary dictionary component.
# 2) params$max_width - maximum allowable width, in pixels, of a preliminary dictionary component.
# 3) params$max_height - maximum allowable hieght, in pixels, of a preliminary dictionary component.
# 4) params$thresh_quantiles - vector of quantiles that will be used for thresholding in SCALPEL step 1.
#    For example, if params$thresh_quantiles is set to [0, 0.001], then the thresholds used by SCALPEL
#    step 1 will be the negative of the minimum of the dF/F-corrected data, the negative of the 0.1%
#    quantile of the dF/F-corrected data, and the average of the two.

# For more detail on these parameteres, see https://rdrr/io/cran/scalpel/man/scalpelStep1.html
# and https://arxiv.org/pdf/1703.06946.pdf.

# For an example of how the JSON parameters file should be formatted, see
# https://github.com/danieldkato/imageProcessing/blob/master/image_segmentation/scalpel/step1/scalpel1_params.json


####################################################################################################
# V. OUTPUTS:

# This script is not (yet) a function and thus has no formal return. However, this wrapper  saves the
# following to secondary storage:

# 1) step1out.Rdata - R object containing information about the output of this processing step. This 
#    can subsequently be loaded into memory and passed as an argument to SCALPEL functions corresponding
#    to later processing steps. 
# 2) metadata.json - a JSON file containing step 1 metadata, including paths and SHA1 digests for inputs 
#    and outputs, as well as parameters.

# In addition to the files mentioned above, scalpelStep1 itself saves a number of files in a directory called
# Step1_<VVVVV>, where <VVVVV> is a 5-digit version number. 


####################################################################################################
library("scalpel")
library("rjson")

# Get path to parameters file from command line input:
args = commandArgs()

jsondata <- fromJSON(file=args[2])
params <- jsondata$params

# step 1 parameters:
step0_dir = jsondata$inputs[[1]]$path
min_size = params$min_size
max_height = params$max_height
max_width = params$max_width
thresh_quantiles = params$thresh_quantiles
print(thresh_quantiles)

# Set working directory to top-level segmentation directory:
setwd(dirname(step0_dir))

# Load dFF data (necessary for computing thresholds for step 1):
dFF_data_name = paste(step0_dir,"Ydeltaf_part1.rds",sep="/")
Y = readRDS(dFF_data_name)
load(paste(step0_dir,"step0out.Rdata",sep="/"))

# Compute thresholds for step 1:
thresh_values = unname(quantile(Y,probs=thresh_quantiles))
thresh_values <- c(thresh_values, mean(thresh_values))
thresh_values <- -1*thresh_values
print(thresh_values)

# Perform scalpel step 1:
step1out = scalpelStep1(step0Output=step0out,minSize=min_size,maxHeight=max_height,maxWidth=max_width,thresholdVec=thresh_values)

# Save step1out object in most recent Step1 directory:
step1_indices = unlist(grepl("Step1",dir()))
step1_dirs = dir()[step1_indices]
mtimes = file.info(step1_dirs)$mtime
most_recent_step1 = step1_dirs[which.max(mtimes)]
save(step1out,file=paste(most_recent_step1,"step1out.Rdata",sep="/"))
start = Sys.time()
date = format(start,format="%Y-%m-%d")
time = strftime(start,format="%H:%M:%S")

# Save step1 metadata:
step1_metadata = list()

step0_path = paste(step0_dir,"step0out.Rdata",sep="/")
sys_out = system(paste("sha1sum",step0_path,sep=" "),intern=TRUE)
step0_sha1 = substr(sys_out,1,40)
inpt1 = list(path=step0_path,sha1=step0_sha1)
step1_metadata$inputs <- list(inpt1)

step1_metadata$params$max_height = max_height
step1_metadata$params$max_width = max_width
step1_metadata$params$thresh_quantiles = thresh_quantiles
step1_metadata$params$thresh_values = thresh_values

step1_full_path = paste(dirname(step0_dir),most_recent_step1,"step1out.Rdata",sep="/")
sys_out2 = system(paste("sha1sum ",step1_full_path,sep=" "),intern=TRUE)
step1_sha1 = substr(sys_out2,1,40)
outpt1 = list(path=step1_full_path,sha1=step1_sha1)
step1_metadata$outputs <- list(outpt1)

step1_metadata$date = date
step1_metadata$time = time
m <- toJSON(step1_metadata)
write(m,file=paste(most_recent_step1,"metadata.json",sep="/"))
