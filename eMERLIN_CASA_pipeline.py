## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
################

##Inputs##
inbase ='20nov2016_CASA_calib'
refant = 'Mk2'
targets = ['0716+4708']	#1236+621 	# List of targets (comma-separated if more than one)
widetargets = 'B1938+6648'	#1236+621 # List of targets for wide-field imaging (comma-separated if more than one)
phsrefs = ['0720+0720'] 	# List pf phase cals (comma-separated if more than one)
fluxcals = ['1331+305','0319+415']	# List of flux cals (comma-separated... you get the idea)
bpasscals = ['1407+284']	# List of bandpass cals (as above)
pointcals = ['1407+284']# List of point cals (as above, although unlikely to be >1)
##########

##Functions##
def dfluxpy(freq,baseline):


    lowest_freq = 300.0;
    highest_freq = 50000.0;
    if (freq < lowest_freq or freq > highest_freq):
        print "Frequency must be between $lowest_freq and $highest_freq MHz. \n"

    # Old values for 3C286
    # A = 1.23734
    # B = -0.43276
    # C = -0.14223
    # D = 0.00345

    # New values taken from AIPS SETJY 31DEC11
    # Values as of 2010

    # A = 1.2361
    # B = -0.4127
    # C = -0.1864
    # D = 0.0294

    # Perley & Butler 2012 values
    A = 1.2515
    B = -0.4605
    C = -0.1715
    D = 0.0336

    log10f = (math.log(freq)/2.3025851) - 3.0; # Why the -3? Because freq has to be GHz for the formula to work.
    log_flux = A + B*log10f + C*log10f*log10f + D*log10f*log10f*log10f
    vlaflux = math.pow(10.0,log_flux)




    # The VLA flux must now be corrected to account for the higher resolving power of merlin. The formula used was obtained with the help of Peter Thomasson. If we assume that 3C286 is represented by a gaussian of angular size theta_s, and represent the resolving power of the given baseline as a function of frequency f and antenna separation L by theta_b(f,L), then the reduction in central flux density A(0) due to the finite theta_s is given by
    #
    #	                           1
    #	            -----------------------------------
    #	 A'(0)       2 pi (theta_b(f,L)^2 + theta_s^2)
    #	------- = --------------------------------------- ,
    #	 A(0)                      1
    #	                   ---------------------
    #	                    2 pi theta_b(f,L)^2
    #
    #	               1
    #	        = -------------- ,
    #	           1 + rho(f,L)
    #
    # where the resolved fraction rho(f,L) is given by
    #
    #	              theta_s^2
    #	rho(f,L) = ---------------- .
    #	            theta_b(f,L)^2
    #
    # Use of theta_b(f,L) = k/(fL) allows this to be written
    #
    #	           (   f*L     )^2
    #	rho(f,L) = (-----------)   * rho_ref .
    #	           (f_ref*L_ref)
    #
    # The reference value of rho is fixed at 0.04 for the MK-TA baseline at 5 GHz (Peter Thomasson).

    ref_bl_length = 11236.79 # MK-TA separation in metres.
    ref_freq = 5000.0
    ref_rho = 0.04
    thisbl = "this baseline (Mk-Ta)"

    bl_length = baseline
    # bl_str = sprintf "%8.2f", $ref_bl_length;


    frac = (freq / ref_freq) * (bl_length / ref_bl_length)
    rho = frac * frac * ref_rho
    merlinflux = vlaflux / (1.0 + rho)

    # Another useful quantity is the resolved percentage:
    #
    resolved_percent = 100.0 * rho / (1.0 + rho)
    caution_res_pc = 10.0

    return vlaflux, merlinflux, resolved_percent, caution_res_pc, thisbl





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
mystep = 1
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]

	os.system('rm -r '+inbase+'.ms')
	importuvfits(fitsfile=inbase+'.fits',vis=inbase+'.ms')


mystep = 2
if(mystep in thesteps):
	print 'Step ', mystep, step_title[mystep]
	
	os.system('rm -r '+inbase+'_han.ms')
	hanningsmooth(vis=inbase+'.ms',outputvis=inbase+'_han.ms',datacolumn='data')
	flagdata(vis=inbase+'_han.ms',mode='manual',autocorr=True)


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
	gaincal(vis=inbase+'.mms', gaintype='K',field=','.join(list(set(phsrefs+fluxcals+bpasscals+pointcals))), caltable=inbase+'mms.K0', refant=refant, solint='inf', minblperant=3, minsnr=3) 
	
