## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math 
import eMERLIN_CASA_functions as em
from casa import *

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

	x = au.timeOnSource(inbase+'_han.ms')
	y = []

	sourcenames = []
	for i in range(len(x.keys())-6):
		sourcenames = sourcenames + [x[i]['source_name']+'.ms']
		if os.path.isfile(x[i]['source_name']+'.rfis')==False:
			y=y+[x[i]['source_name']]
			print y
	if len(y) != 0:
		for i in range(len(y)):
			print 'Missing rfistrategy for: '+y[i]
		print 'Please run step 3 again!'
	else:
		while True:
			s = raw_input('All rfistrategys are there: Proceed?:\n')
			if s == 'yes' or s == 'y':
				for i in range(len(x.keys())-6):
					print 'Flagging: '+x[i]['source_name']+'.ms'+' with strategy: '+x[i]['source_name']+'.rfis'
					os.system('aoflagger -fields '+str(i)+' -strategy '+x[i]['source_name']+'.rfis '+inbase+'_han.ms')
				break
			if s == 'no' or s == 'n':
				sys.exit('Please restart when you are happy')
	
### Convert to mms for parallelisation ###
mystep = 5
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	os.system('rm -r '+inbase+'.mms')
	partition(vis=inbase+'_han.ms',outputvis=inbase+'.mms',createmms=True,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)


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
	
