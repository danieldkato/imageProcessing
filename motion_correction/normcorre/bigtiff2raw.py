import skimage.external.tifffile
#import matplotlib.pyplot as plt
import numpy as np
import os
import pyText

# where to put the output files in binary format
output_directory = '/home/georgia/ToBeMotionCorrected/PM04-2/170404'

# output filename (rawfile)
output_filename = 'rawimage.raw'

# input filename (bigtiff)
filename = '/mnt/nas2/homes/georgia/Imaging/PM04-2/170404/170404-beh4x_00002.tif'

# create an object to read the data
print "creating object"
tf = skimage.external.tifffile.TiffFile(filename) #, pages=range(5000))

# write to raw binary format
print "writing"
chunk_size_frames = 1000
n_frames = len(tf.pages)

# this needs to be less than the minimum value encountered in any frame
# an offset will be added to the data to make this minimum value zero
min_data_value = -5000

with file(os.path.join(output_directory, output_filename), 'wb') as fi:
    for chunk_start in range(0, n_frames, chunk_size_frames):
        # get chunk
        chunk = tf.asarray(slice(chunk_start, chunk_start+chunk_size_frames))
        if chunk.min() < min_data_value:
            raise ValueError("minimum value was %d" % chunk.min())
        
        # Apply offset
        chunk = chunk - min_data_value
        
        # convert to big-endian unsigned short (Eftychios format)
        chunk = chunk.astype('>H')
        
        # Matlab reads in column ("Fortran") order
        # Or maybe not
        bytestring = chunk.tobytes()
        
        # write chunk
        fi.write(bytestring)


#text when finished
"""
msg = '\nConverting BIG tiff complete!'
recipient = '@vtext.com'
pyText.textme(msg,recipient)
"""