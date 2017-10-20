library("scalpel")
library("rjson")

jsondata <- fromJSON(file="/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/image_segmentation/scalpel/step1/scalpel1_params.json")
params <- jsondata$params

# step 1 parameters:
step0_dir = jsondata$inputs[[1]]$path
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
step1out = scalpelStep1(step0Output=step0out,maxHeight=max_height,maxWidth=max_width,thresholdVec=thresh_values)

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
