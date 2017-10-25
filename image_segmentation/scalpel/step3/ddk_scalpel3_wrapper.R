# ddk_scalpel3_wrapper.R

# DOCUMENTATION TABLE OF CONTENTS:
# I. OVERVIEW
# II. REQUIREMENTS
# III. INPUTS
# IV. OUTPUTS

# last updated DDK 2017-10-24


####################################################################################################
# I. OVERVIEW:

# This script is a wrapper for the fourth step ("step 3") of the SCALPEL automated fluorescence movie segmentation
# pipeline. This step is responsible for performing non-negative matrix factorization on the complete movie
# data matrix using the spatial components from the refined dictionary produced by SCALPEL step 2. For more
# detail, see the SCALPEL documentation at https://rdrr.io/cran/scalpel.


####################################################################################################
# II. REQUIREMENTS:

# 1) R.
# 2) The R package SCALPEL. For installation instructions, see https://rdrr.io/cran/scalpel.
# 3) The R package rjson. To install, call install.packages("rjson") from inside R. 


####################################################################################################
# III. INPUTS:

# This script is not (yet) a function and thus takes no formal input arguments. Instead, the user
# must specify a path to a JSON file containing the desired parameters. This parameters file must
# include the following fields:

# 1) params$lambda_method - description of how lambda should be chose. See man pages (link below) for details.
# 2) params$lambda - value of lambda to use when fitting sparse group lasso. If not NULL, lambda_method will not be used.
# 3) params$min_cluster_size - minimum number of preliminary spatial compnents (from step 1) that a cluster must 
#    contain in order to be included in the sparse group lasso. 
# 4) params$alpha - alpha value for fitting sparse group lasso. 
# 5) params$remove_border - Boolean variable specifying whether or not to exclude sptial components within 10 pixels of the frame border.
# 6) params$exclude_reps - variable specifying whether or not to include spatial components flagged for exclusion during manual curation
#    between steps 2 and 3. Set to "dicarded" to do so,"NULL" otherwise.

# For more detail on these parameters, see https://rdrr.io/cran/scalpel/man/scalpelStep3.html and
# https://arxiv.org/pdf/1703.06946.pdf.

# For an example of how the JSON parameters file should be formatted, see
# https://github.com/danieldkato/imageProcessing/blob/master/image_segmentation/scalpel/step0/scalpel0_params.json


####################################################################################################
# IV. OUTPUTS:

# This script is not (yet) a function and thus has no formal return. However, this wrapper  saves the
# following to secondary storage:

# 1) step3out.Rdata - R object containing information about the output of this processing step. This 
#    can subsequently be loaded into memory and passed as an argument to SCALPEL functions corresponding
#    to later processing steps. 
# 2) metadata.json - a JSON file containing step 3 metadata, including paths and SHA1 digests for inputs 
#    and outputs, as well as parameters. This includes a field called excluded_ROIs, which is a vector
#    of column indices specifying which columns of step2out$A (i.e.n spatial components) were excluded 
#    from the current analysis after having been flagged as false neurons during manual curation.

# In addition to the files mentioned above, scalpelStep3 itself saves a number of files in a directory called
# Step3_<params>, where <params> is a long string specifying the parameters used in the current analysis. 


####################################################################################################
#Load packages:
library("scalpel")
library("R.matlab")
library("rjson")

# Load parameters:
jsondata <- fromJSON(file="/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/image_segmentation/scalpel/step3/scalpel3_params.json")
inputPath = jsondata$inputs[[1]]$path
params = jsondata$params

# Take care of some JSOIN->R formatting things:
if (params$lambda=="NULL"){
	params$lambda=NULL
}
if (params$exclude_reps=="NULL"){
	params$exclude_reps=NULL
}

# Load step 2 output:
load(inputPath)
step2_dir = dirname(inputPath)
setwd(step2_dir)

# Perform SCALPEL step 3:
step3out = scalpelStep3(step2Output=step2out,lambdaMethod=params$lambda_method,minClusterSize=params$min_cluster_size,alpha=params$alpha,removeBorder=as.logical(params$remove_border),excludeReps=params$exclude_reps)

# Find the most recent "Step2" directory:
step3_indices = unlist(grepl("Step3",dir())) # get indices of all directories whose names contain the string "Step2"
step3_dirs = dir()[step3_indices] # get the names of all directories containing the string "String2"
mtimes = file.info(step3_dirs)$mtime # get the time each of these folders was last modified
now = Sys.time() # get current time
time_diffs = now - mtimes # find out how long ago each was last modified
most_recent_ind = which.min(time_diffs)
newest_step3_dir = step3_dirs[most_recent_ind]
setwd(newest_step3_dir) # cd to the msot recent step 2 directory to save some additional outputs there

# Save step3out as an .Rdata object:
s3obj_full_path = paste(step2_dir,newest_step3_dir,"step3out.Rdata",sep="/")
save(step3out,file=s3obj_full_path)
now = Sys.time()
date = format(now,format="%Y-%m-%d") 
time = strftime(now,format="%H:%M:%S")

# Write metadata:

# Get input metadata:
step3_metadata = list()
sys_out = system(paste("sha1sum",inputPath,sep=" "),intern=TRUE)
input_sha1 = substr(sys_out,1,40)
inpt1 <- list(path=inputPath,sha1=input_sha1)
step3_metadata$inputs <- list(inpt1)

# Write parameters metadata:
step3_metadata$params <- params

# Include which ROIs from original step2out$A ROI set were excluded:
keep = readRDS(paste(step2_dir,"Step2Data","keep.rds",sep="/"))
excluded = which(keep == "no")
step3_metadata$excluded_ROIs <- excluded

# Write outpyut metadata:
sys_out2 = system(paste("sha1sum",s3obj_full_path,sep=" "),intern=TRUE)
output1_sha1 = substr(sys_out2,1,40)
output1 <- list(path=s3obj_full_path,sha1=output1_sha1)
step3_metadata$outputs <- list(output1)

step3_metadata$date <- date
step3_metadata$time <- time

m <- toJSON(step3_metadata)
write(m,file=paste(step2_dir,newest_step3_dir,"metadata.json",sep="/"))

# Return to previous working directory:
setwd(step2_dir)

