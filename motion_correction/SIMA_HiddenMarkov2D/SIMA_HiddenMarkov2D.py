"""
Last updated DDK 2016-10-12

OVERVIEW:
This script defines a wrapper function for motion-correcting a raw multi-page TIFF using SIMA's 2D Hidden Markov approach. For full documentation on SIMA, see:

http://www.losonczylab.org/sima/1.2/index.html

The purpose of this wrapper is to  provide a layer of abstraction that is low enough to permit processing a single input movie (to be used by different batch-processing scripts, or for debugging purposes, for example), but high enough that the processing adheres to the metadata protocols entailed by my modular workflow, described here:  

10.112.43.46\Public\dank\multiSens\analysis\README.txt

REQUIREMENTS:
1) SIMA - for installation instructions, see http://www.losonczylab.org/sima/1.2/install.html
2) writeMetadata.py, 

"""

import sima
import os

def SIMA_HiddenMarkov2D (input, output, gran, md):
    
    rawTiffName = os.path.split(input)[1][0:-4]
    mcTiffName = rawTiffName + '_mc'
    
    mc_approach = sima.motion.HiddenMarkov2D(granularity=gran, max_displacement=md, verbose=True)
    sequence = [sima.Sequence.create('TIFF', input)]
    dataset = mc_approach.correct(sequence, output)
    dataset.export_frames([[[mcTiffName + '.tif']]])
    


    
    