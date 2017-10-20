library("scalpel")
library("rjson")


jsondata <- fromJSON(file="/mnt/nas2/homes/dan/code_libraries/ddk_image_processing/image_segmentation/scalpel/step0/scalpel0_params.json")
params <- jsondata$params

# ste 0 parameters:
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
