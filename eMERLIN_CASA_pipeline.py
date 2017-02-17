## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math 
import eMERLIN_CASA_functions as em
from casa import *
import Tkinter,tkFileDialog
################

##Inputs##
inbase = '20nov2016_CASA_calib'          # Prefix of measurement set
aoflagger_execute = 'aoflagger'         # Command line executable of aoflagger
wsclean_execute   = 'wsclean-2.2'       # Command line executale of wsclean
refant = 'Mk2'                          # Reference antenna
targets = ['0716+4708']	                # List of targets (comma-separated if more than one)
phsrefs = ['0720+0720'] 	        # List pf phase cals (comma-separated if more than one)
fluxcals = ['1331+305','0319+415']	# List of flux cals (comma-separated... you get the idea)
bpasscals = ['1407+284']	        # List of bandpass cals (as above)
pointcals = ['1407+284']                # List of point cals (as above, although unlikely to be >1)
##########


thesteps = []
step_title = {1: 'Convert into measurement set',
              2: 'Hanning smoothing',
              3: 'Rfigui strategies',
	      4: 'AOflag with defined strategies',
              5: 'Convert into mms',
              6: 'What\'s in your data?',
	      7: 'Delay correction',
}


thesteps = []
for i in range(len(step_title)):
	print 'Step ', i+1, step_title[i+1]
print '\n'
try:
	print 'List of steps to be executed ...', mysteps
	thesteps = mysteps
except:
	print 'global variable mysteps not set'

while True:
	for i in range(len(mysteps)):
		print 'Step ', mysteps[i], step_title[mysteps[i]]
	print '\n'
	s = raw_input('Are these the steps you want to conduct (yes or no): ')
	if s == 'yes' or s == 'y':
		break
	if s == 'no' or s == 'n':
		sys.exit('Please restate the mysteps parameter')

if (thesteps==[]):
	thesteps = range(0,len(step_title))
	print 'Executing all steps: ', thesteps

### Import the data into a measurement set ###

fitsfile = inbase+'.fits'
vis = inbase+'.ms'

mystep = 1
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	em.run_importuvfits(fitsfile,vis)


mystep = 2
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
	em.hanningflag(inputvis=vis,deloriginal=True)
	


mystep = 3
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
	
	os.system('rfigui '+inbase+'_han.ms')

mystep = 4
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
	
	em.run_aoflagger(vis=vis,mode='user')
	
### Convert to mms for parallelisation ###
mystep = 5
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	em.ms2mms(vis=vis,mode='parallel')


mystep = 6
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	os.system('rm -rf '+inbase+'.mms.listobs')
	listobs(vis=inbase+'.mms',
            listfile=inbase+'.mms.listobs')

mystep = 7
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
	
	os.system('rm -rf '+inbase+'.mms.K0')
	gaincal(vis=inbase+'.mms', gaintype='K',field=','.join(list(set(phsrefs+fluxcals+bpasscals+pointcals))), caltable=inbase+'.mms.K0', refant=refant, solint='inf', minblperant=3, minsnr=3) 
	
