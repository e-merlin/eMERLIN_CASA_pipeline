## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
import eMERLIN_CASA_functions as em
import eMERLIN_CASA_GUI as emGUI
from casa import *
import os
from casa import table as tb
from casa import ms
from Tkinter import *
################

##Inputs##
inputs, processes = emGUI.GUI_pipeline().confirm_parameters()
print inputs
print processes
##########

if inputs['quit'] == 1:
	sys.exit()

fitsfile = inputs['inbase']+'.fits'
vis = inputs['inbase']+'.ms'

if processes['run_importuvfits'] == 1:
	em.run_importuvfits(fitsfile,vis)


if processes['hanningflag'] == 1:
	em.hanningflag(inputvis=vis,deloriginal=True)

if processes['autoflag'] == 1:
    if processes['rfigui'] == 1:
	os.system('rfigui '+vis)
        em.run_aoflagger(vis=vis,mode='user')
    if processes['rfigui']== 0:
	em.run_aoflagger(vis=vis,mode='default')


### Convert to mms for parallelisation ###
if processes['ms2mms'] == 1:
	em.ms2mms(vis=vis,mode='parallel')

vis = inputs['inbase']+'.mms'

if processes['do_prediag'] == 1:
	em.do_prediagnostics(vis)

'''

	os.system('rm -rf '+inbase+'.mms.K0')
	gaincal(vis=inbase+'.mms', gaintype='K',field=','.join(list(set(phsrefs+fluxcals+bpasscals+pointcals))), caltable=inbase+'.mms.K0', refant=refant, solint='inf', minblperant=3, minsnr=3)
'''
