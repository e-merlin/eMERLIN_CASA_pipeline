## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
sys.path.insert(0,'./CASA_eMERLIN_pipeline') # add github path at runtime
import os,sys,math
import eMERLIN_CASA_functions as em
import eMERLIN_CASA_GUI as emGUI
from casa import *
from casa import table as tb
from casa import ms
from Tkinter import *
import getopt
################

##Inputs##
inputs = em.check_in()
data_dir = em.backslash_check(inputs['data_dir'])
plots_dir = em.backslash_check(inputs['plots_dir'])
print inputs
##########

if inputs['quit'] == 1: #Check from GUI if quit is needed
	sys.exit()

fitsfile = inputs['inbase']+'.fits'
vis = data_dir+inputs['inbase']+'.ms'

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
	print data_dir
	em.run_importfitsIDI(data_dir,vis)

if inputs['hanning'] == 1:
	em.hanning(inputvis=vis,deloriginal=True)

if inputs['autoflag'] == 1:
    if inputs['rfigui'] == 1:
	os.system('rfigui '+vis)
        em.run_aoflagger(vis=vis,mode='user')
    if inputs['rfigui']== 0:
	em.run_aoflagger(vis=vis,mode='default')


### Convert to mms for parallelisation ###
if inputs['ms2mms'] == 1:
	em.ms2mms(vis=vis,mode='parallel')

if os.path.isdir('./'+inputs['inbase']+'.mms') == True: #takes into account parallel or not
	vis = inputs['inbase']+'.mms'

if inputs['do_prediag'] == 1:
	em.do_prediagnostics(vis,plots_dir)

'''

	os.system('rm -rf '+inbase+'.mms.K0')
	gaincal(vis=inbase+'.mms', gaintype='K',field=','.join(list(set(phsrefs+fluxcals+bpasscals+pointcals))), caltable=inbase+'.mms.K0', refant=refant, solint='inf', minblperant=3, minsnr=3)
'''
