library("scalpel")
library("R.matlab")
library("rjson")

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

# Save step2out as an .Rdata object:
s3obj_full_path = paste(step2_dir,newest_step3_dir,"step3out.Rdata",sep="/")
save(step3out,file=s3obj_full_path)

# TODO: Convert A.mat to .txt?

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

