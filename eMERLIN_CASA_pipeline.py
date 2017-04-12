## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
import numpy as np
from casa import table as tb
from casa import ms
from Tkinter import *
import getopt
import logging

# Find path of pipeline to find external files (like aoflagger strategies or emerlin-2.gif)
pipeline_path = os.path.dirname(sys.argv[np.where(np.asarray(sys.argv)=='-c')[0][0] + 1]) + '/'
sys.path.append(pipeline_path)
import functions.eMERLIN_CASA_functions as em
import functions.eMERLIN_CASA_GUI as emGUI

# Setup logger
logger = logging.getLogger('logger')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('eMCP.log') # create a file handler
handler.setLevel(logging.INFO)
formatter = logging.Formatter(fmt='%(asctime)s.%(msecs)1d | %(levelname)s | %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler) # add the handlers to the logger
consoleHandler = logging.StreamHandler() # Stream errors to terminal also
consoleHandler.setFormatter(formatter)
logger.addHandler(consoleHandler)

logger.info('Starting pipeline')
logger.info('Running pipeline from: {}'.format(pipeline_path))

##Inputs##
inputs = em.check_in(pipeline_path)
data_dir = em.backslash_check(inputs['data_dir'])
plots_dir = em.backslash_check(inputs['plots_dir'])
calib_dir = em.backslash_check(inputs['calib_dir'])
logger.info('Inputs used: {}'.format(inputs))

refant = inputs['refant']

## Create directory structure ##
em.makedir(plots_dir)
em.makedir(calib_dir)

if inputs['quit'] == 1: #Check from GUI if quit is needed
    logger.debug('Pipeline exit')
    sys.exit()

fitsfile = inputs['inbase']+'.fits'
vis = inputs['inbase']+'.ms'


## Check for measurement sets in current directory otherwise drag from defined data directory
if os.path.isdir(inputs['inbase']+'.ms') == False and os.path.isdir(inputs['inbase']+'.mms') == False:
	if os.path.isdir(data_dir+inputs['inbase']+'.mms') == True:
		os.system('rsync -ar --progress {0} ./'.format(data_dir+inputs['inbase']+'.mms'))
	elif os.path.isdir(data_dir+inputs['inbase']+'.ms') == True:
		os.system('rsync -ar --progress {0} ./'.format(data_dir+inputs['inbase']+'.ms'))
	else:
		logger.info('No measurement set found, assuming you need to importfits or change data dir')
else:
	logger.info('Measurement set found: {}. Continuing with your inputs'.format(vis))

## Pipeline processes, inputs are read from the inputs dictionary
if inputs['run_importfits'] == 1:
	em.run_importfitsIDI(data_dir,vis)

if inputs['hanning'] == 1:
	em.hanning(inputvis=vis,deloriginal=True)

if inputs['rfigui'] == 1:
    em.run_rfigui(vis)

### Convert MS to MMS ###
if inputs['ms2mms'] == 1:
    if em.check_aoflagger_version():
        em.ms2mms_fields(msfile=vis)
    else:
        em.ms2mms(vis=vis,mode='parallel')

## check for parallelisation
if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
    vis = inputs['inbase']+'.mms'

if inputs['autoflag'] == 1:
    em.run_aoflagger_fields(vis=vis,fields='all', pipeline_path = pipeline_path)

### Produce some initial plots ###
if inputs['do_prediag'] == 1:
	em.do_prediagnostics(vis,plots_dir)


### Delay calibration ###
if inputs['do_delay'] == 1:
    caltable_name = inputs['inbase']+'_delay.K'
    em.solve_delays(vis,caltable_name=caltable_name,calsources='',solint='600s',refant=refant,combine='spw',spw='',caldir=calib_dir,plotdir=plots_dir)




logger.info('Pipeline finished')
logger.info('#################')


