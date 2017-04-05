## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
import functions.eMERLIN_CASA_functions as em
import functions.eMERLIN_CASA_GUI as emGUI
from casa import table as tb
from casa import ms
from Tkinter import *
import getopt
from tasks import *
from casa import *


###############

##Inputs##
inputs = em.check_in()
data_dir = em.backslash_check(inputs['data_dir'])
plots_dir = em.backslash_check(inputs['plots_dir'])
print inputs
##########

if inputs['quit'] == 1: #Check from GUI if quit is needed
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
		print 'No measurement set found, assuming you need to importfits or change data dir'
else:
	print 'Measurement set found: '+vis+'. Continuing with your inputs'

## Pipeline processes, inputs are read from the inputs dictionary
if inputs['run_importfits'] == 1:
	em.run_importfitsIDI(data_dir,vis)

if inputs['hanning'] == 1:
	em.hanning(inputvis=vis,deloriginal=True)

### Flagging and parallelisation ###
### Check AOflagger version. Decide if old or new procedure is needed. ###
### if v2.7< do ms2mms fields
if em.check_aoflagger_version(): #Use J. Moldon's default flagger for mms architecture when -field parameter doesnt exist
	if inputs['ms2mms'] == 1:
		em.ms2mms_fields(msfile=vis)
		if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
			vis = inputs['inbase']+'.mms'
	if inputs['autoflag'] == 1:
	    em.run_aoflagger_fields(vis=vis,fields='all')

### if v2.9+ use -field parameter and can generate source specific rfi strategies
else: ##run aoflagger on .ms file first so that gui works properly if aoflagger >2.9
	if inputs['autoflag'] == 1:
		if inputs['rfigui'] == 1:
			os.system('rfigui '+vis)
			em.run_aoflagger(vis=vis,mode='user')
		if inputs['rfigui']== 0:
			em.run_aoflagger(vis=vis,mode='default')
	if inputs['ms2mms'] == 1:
		em.ms2mms(vis=vis,mode='parallel')

## check for parallelisation
if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
			vis = inputs['inbase']+'.mms'


### Produce some initial plots ###
if inputs['do_prediag'] == 1:
	em.do_prediagnostics(vis,plots_dir)
