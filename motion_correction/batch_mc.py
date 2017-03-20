"""
Last updated DDK 9/12/2016

OVERVIEW:
This script defines a function for motion-correcting batches of raw tiff stacks using the SIMA software package. 


REQUIREMENTS:
1) Python >=2.7
2) SIMA


INSTRUCTIONS: 
Launch a Python instance from which it is possible to import this script as a module, then import the function defined herein using

from batch_mc import batchMC

The function batchMC takes as its sole argument a list of pairs. The first element of each pair should be the name of a mouse, and the second element of each pair should be a list of imaging sessions for that mouse and from which the raw data should be motion-corrected. 

All mouse and session names should be formatted as strings. By convention, I have been naming sessions after the date they were performed using the format 'YYMMDD', where YY stands for the last two digits of the year, MM is the two-digit month number, and DD is the two-digit day number. Names of mice and sessions should exactly match the names of folders containing data from those mice and sessions, respectively. 


DESCRIPTION:
This script searches for raw tiffs from a list of requested imaging sessions. In order to do this, it makes a number of assumptions about the file structure used to store the data, which may be useful to describe here.

Raw data are currently (9/12/2016) stored on 10.112.43.46 (a.k.a. 'NAS2', which is mapped to as network drive 'Y:' on host 'Build3', i.e. my desktop) under the following file structure:

Public/
	dank/
		mouse1/
			2P/
				session1/
					site1/
						grab001/
							dat/
							ardulines
							grab001.tif
							meta.txt
						grab002/
						zstack/
					site2/
					site3/
					craniotomy.tif
					sitemap.jpg
					sitemap.psd
				session2/
				session3/
			intrinsic
			simtemap.jpg
			sitemap.psd
		mouse2/
		mouse3/

Public/dank contains one folder for each mouse. Each mouse folder contains one folder for each imaging session taken from that mouse (possibly along with some image files depicting a master site map that includes barrel locations and surface images from every 2-photon imaging location registered to a widefield image of the craniotomy, which can be ignored for these purposes).  

"""


import os
import sima

basePath = 'Y:\\dank' 
mice = os.listdir(basePath)
mc_approach = sima.motion.HiddenMarkov2D(granularity='row', max_displacement=[20,30], verbose=True)

def batchMC(requestedData):
	
	mcQueue = [] #BEFORE actually doing any motion correction, attempt to find all of the requested files and assemble their locations into a list. If a mouse or session is missing,  throw up an error and ask the user to either re-enter the name of the mouse or session omit it from the list. Only AFTER assembling the list, iterate through it and motion correct all of the grabs. 
	
	for pair in requestedData:
		requestedMouse = pair[0]
		requestedSessions = pair[1] 
		mouseFound = 0
		
		while not mouseFound: #For each mouse in requestedData, search the base directory for the folder corresponding to that mouse...
			if requestedMouse in mice:
				mousePath = os.path.join(basePath, requestedMouse)
				mouseFound = 1
			else:
				response = raw_input('Directory for ' + mouse + ' not found. Re-enter mouse name here or hit Enter to skip mouse.')
				if response is not '':
					requestedMouse = response
				elif response is '':
					print(mouse + ' not found, will omit from motion correction.')
					break 
		
		if mouseFound: # ... and if the folder corresponding to that mouse is found, then for each session requested for that mouse, search for the corresponding session folder...
			sessions = os.listdir(mousePath)
			for requestedSession in requestedSessions: 
				sessionFound = 0
				while not sessionFound:
					if requestedSession in sessions:
						sessionPath = os.path.join(mousePath, requestedSession)
						sessionFound = 1
					else:
						response = raw_input('Imaging session data for ' + date + ' not found. Re-enter session name or press Enter to skip session.')
						if response is not '':
							requestedSession = response
						elif response is '':
							print('Imaging session data for ' + date + ' not found, will omit from motion correction.')
							break
				
				if sessionFound: #... and if the folder corresponding to that session is found... 
					imgSites = [x for x in os.listdir(sessionPath) if 'site' in x] #... get a list of every site imaged in that session, and for each site...
					for imgSite in sites:
						sitePath = os.path.join(sessionPath, imgSite)
						grabs = [y for y in os.listdir(sitePath) if 'grab' in y] #... get a list of every grab taken at that site, and for each grab... 
						for grab in grabs: 
							grabPath = os.path.join(sitePath, grab)
							tiffs = [z for z in os.listdir(grabPath) if '.tif' in z] #... get a list of every .tif in that grab folder (there should really only be one, and it should have the same name as the grab folder, but it's good to check)...
							for tiff in tiffs:
								mcQueue.append(os.path.join(grabPath, tiff)) #... then finally, append the path to each .tif to the queue of .tif files to be motion-corrected. 
								
	for rawTiffPath in mcQueue: #Once the list of tiffs has been assembled, iterate through it and motion correct all of the tiffs
		rawTiffDirectory = os.path.dirname(rawTiffPath)
		rawTiffName = os.path.split(rawTiffPath)[1][0:-4] #omit the .tif file extension
		mcTiffName = rawTiffName + '_mc'
		sequence = [sima.Sequence.create('TIFF', rawTiffPath)]
		dataset = mc_approach.correct(sequence, os.path.join(rawTiffDirectory, mcTiffName))
		os.chdir(os.path.join(rawTiffDirectory, mcTiffName+'.sima'))
		dataset.export_frames([[[mcTiffName + '.tif']]])
								
						
				
						
			
		
		
		
		
		
		
		
		
		
		
		
		
	