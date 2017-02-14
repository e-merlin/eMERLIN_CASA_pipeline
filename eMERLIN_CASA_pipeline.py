## v0.00001 of an eMERLIN CASA pipeline ##
import os,sys
import dfluxpy-perleybutler, delay_checker
##Inputs##
inbase ='20nov2016_concat_UV20nov2016_concat_UV'
targets = '0716+4708, 1557+3721'	#1236+621 	# List of targets (comma-separated if more than one)
widetargets = 'B1938+6648'	#1236+621 # List of targets for wide-field imaging (comma-separated if more than one)
phsrefs = '0720+0720' 	# List pf phase cals (comma-separated if more than one)
fluxcals = '1331+305','0319+415'	# List of flux cals (comma-separated... you get the idea)
bpasscals = '1407+284'	# List of bandpass cals (as above)
pointcals = '1407+284'	# List of point cals (as above, although unlikely to be >1)
##########
thesteps = []
step_title = {1: 'What\'s in the data',
              2: 'Convert to mms for parallelisation',
	          3: 'Setjy fluxes',
	          4: 'Dropouts: Look for bad data amp vs. uvdist',
	          5: 'Rfigui',
	          6: 'AOFlag the data',
	          7: 'Run VLA pipeline',
	          8: 'RFI post cal: amp vs. uvdist',
	          9: 'RFI post cal: amp vs. frequency',
	         10:'RFI post cal: amp vs. uvdist',
	         11:'Re-calibrate?',
	         12:'Split target',
	         13:'Image target',
             14:'Self Calibration: image',
	         15:'Self Calibration: gaincal',
	         16:'Self Calibration: apply',
             17:'Closure phase',
             18:'Direction-dependent?'
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
	if s == 'yes':
		break
	if s == 'no':
		sys.exit('Please resabort tate the mysteps parameter')

if (thesteps==[]):
	thesteps = range(0,len(step_title))
	print 'Executing all steps: ', thesteps

mystep = 1
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	os.system('rm -rf '+inbase+'.ms.listobs')
    listobs(vis=inbase+'.ms',
            listfile=inbase+'.ms.listobs')

mystep = 2
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

    os.system('rm -r '+inbase+'.mms')
    partition(vis=inbase+'ms',outputvis=inbase+'.mms',createmms=True,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)

mystep = 3
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
