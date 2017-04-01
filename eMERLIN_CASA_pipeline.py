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
data_dir = inputs['data_dir']
plots_dir = inputs['plots_dir']
#inputs = emGUI.GUI_pipeline().confirm_parameters()
print inputs
##########

if inputs['quit'] == 1:
	sys.exit()

fitsfile = inputs['inbase']+'.fits'
vis = inputs['inbase']+'.ms'

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

vis = inputs['inbase']+'.mms'

if inputs['do_prediag'] == 1:
	em.do_prediagnostics(vis)

'''

	os.system('rm -rf '+inbase+'.mms.K0')
	gaincal(vis=inbase+'.mms', gaintype='K',field=','.join(list(set(phsrefs+fluxcals+bpasscals+pointcals))), caltable=inbase+'.mms.K0', refant=refant, solint='inf', minblperant=3, minsnr=3)
'''
