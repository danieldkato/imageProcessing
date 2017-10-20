library("scalpel")
library("R.matlab")
library("rjson")

jsondata <- fromJSON(file="/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/image_segmentation/scalpel/step2/scalpel2_params.json")
inputPath = jsondata$inputs[[1]]$path
cutoff = jsondata$params$cutoff
omega = jsondata$params$omega

load(inputPath)
step1_dir = dirname(inputPath)
setwd(step1_dir)

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

