#!/usr/local/python
import os
from casa import table as tb
from casa import ms
import numpy as np
from Tkinter import *
import tkMessageBox
import sys, shutil
import getopt
#from task_importfitsidi import *
from eMERLIN_CASA_GUI import GUI_pipeline
# Need to be in this order:
from tasks import *
from casa import *

import logging
logger = logging.getLogger('logger')

def check_in(pipeline_path):
	try:
		opts, arg = getopt.getopt(sys.argv[1:],'i:c:hg',['help','input=','gui'])
		logger.debug(sys.argv[1:])
	except getopt.GetoptError as err:
		logger.error(err)
		sys.exit(2)
	for o,a in opts:
		logger.debug('{0} {1}'.format(o,a))
		if o in ('-i','--input'):
			inputs = headless(a) ## read input file
			inputs['quit'] = 0 ##needed to add to be compatible with GUI
			logger.info('inputs from file: {}'.format(inputs))
		elif o in ('-g','--gui'):
			inputs = GUI_pipeline(pipeline_path).confirm_parameters() ## read input file
			logger.info('inputs from GUI: {}'.format(inputs))
		elif o in ('-h','--help'):
			logger.debug('help will be written soon')
			sys.exit()
		elif o == '-c':
			logger.debug('Executing!')
		else:
			assert False, "rerun with either headless -i or gui" #if none are specifed run GUI
	return inputs

def backslash_check(directory):
	if directory[-1] != '/':
		return directory+'/'
	else:
		return directory


def headless(inputfile):
	''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
	INPUTFILE = open(inputfile, "r")
	control = {}
	# a few useful regular expressions
	newline = re.compile(r'\n')
	space = re.compile(r'\s')
	char = re.compile(r'\w')
	comment = re.compile(r'#.*')
	# parse the input file assuming '=' is used to separate names from values
	for line in INPUTFILE:
		if char.match(line):
			line = comment.sub(r'', line)
			line = line.replace("'", '')
			(param, value) = line.split('=')
			param = newline.sub(r'', param)
			param = param.strip()
			param = space.sub(r'', param)
			value = newline.sub(r'', value)
			value = value.strip()
			valuelist = value.split(', ')
			if len(valuelist) == 1:
				if valuelist[0] == '0' or valuelist[0]=='1':
					control[param] = int(valuelist[0])
				else:
					control[param] = str(valuelist[0])
			else:
				control[param] = str(valuelist)
	return control

def Tkinter_select():
	root = Tkinter.Tk()
	root.withdraw()
	file = tkFileDialog.askdirectory(parent=root,mode='rb',title='Choose a file')
	if file != None:
    		print file
	return file

def check_history(vis):
	tb.open(vis+'/HISTORY')
	x = tb.getcol('MESSAGE')
	y = [i for i, item in enumerate(x) if 'eMER_CASA_Pipeline:' in item]
	if len(y) == 0:
		print 'Measurement set has not been processed \n'
	else:
		print 'WARNING: Some pipeline processes have already been run'
		for i in range(len(y)):
			print x[y[i]]

def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass

def rmdir(pathdir,message='Deleted:'):
    try:
        shutil.rmtree(pathdir)
        logger.info('{0} {1}'.format(message, pathdir))
    except:
        logger.debug('Could not delete: {0} {1}'.format(message, pathdir))
        pass

def rmfile(pathdir,message='Deleted:'):
    try:
        os.remove(pathdir)
        logger.info('{0} {1}'.format(message, pathdir))
    except:
        logger.debug('Could not delete: {0} {1}'.format(message, pathdir))
        pass

def run_importfitsIDI(data_dir,vis):
	logger.info('Starting importfitsIDI procedure')
	os.system('rm -r '+vis)
	fitsfiles =[]
	for file in os.listdir(data_dir):
		if file.endswith('fits') or file.endswith('FITS'):
			fitsfiles = fitsfiles + [data_dir+file]
			logger.info('FITS file found to be imported: {0}'.format(file))
	logger.info('Start importfitsIDI')
	importfitsidi(fitsidifile=fitsfiles, vis=vis, constobsid=True, scanreindexgap_s=15.0)
	ms.writehistory(message='eMER_CASA_Pipeline: Import fitsidi to ms, complete',msname=vis)
	logger.info('End importfitsIDI')
	logger.info('Start UVFIX')
	fixvis(vis=vis,outputvis=vis+'.uvfix')
	logger.info('End UVFIX')
	os.system('rm -r {0}'.format(vis))
	os.system('mv {0} {1}'.format(vis+'.uvfix', vis))
	logger.info('Start flagdata_autocorr')
	flagdata(vis=vis,mode='manual',autocorr=True)
	ms.writehistory(message='eMER_CASA_Pipeline: Fixed uv coordinates & remove autocorr',msname=vis)
	logger.info('End flagdata_autocorr')
	logger.debug('You have been transformed from an ugly UVFITS to beautiful MS')
	return

##Hanning smoothing and flag of autocorrelations, will delete original and rename
def hanning(inputvis,deloriginal):
    logger.info('Start hanning')
    if inputvis[-3:].lower() == '.ms':
        outputvis = inputvis[:-3]+'_hanning'+inputvis[-3:]
        os.system('rm -r '+outputvis)
        hanningsmooth(vis=inputvis,outputvis=outputvis,datacolumn='data')
    elif inputvis[-3:].lower() == 'mms':
        outputvis = inputvis[:-4]+'_hanning'+inputvis[-4:]
        os.system('rm -r '+outputvis)
        mstransform(vis=inputvis,outputvis=outputvis,hanning=True,datacolumn='data')
    if deloriginal==True:
        os.system('rm -r {0}'.format(inputvis))
        os.system('rm -r {0}.flagversions'.format(inputvis))
    else:
        os.system('mv {0} {1}'.format(inputvis, inputvis+'_prehanning'))
        os.system('mv {0}.flagversions {1}.flagversions'.format(inputvis, inputvis+'_prehanning'))
    os.system('mv {0} {1}'.format(outputvis, inputvis))
    os.system('mv {0}.flagversions {1}.flagversions'.format(outputvis, inputvis))
    ms.writehistory(message='eMER_CASA_Pipeline: Hanning smoothed data, complete',msname=inputvis)
    logger.info('End hanning')
    return

def run_rfigui(vis):
    logger.info('Start run_rfigui')
    """This function should output the new strategies to /aoflagger_strategies/user/<field>.rfis in 
    either the local folder or the pipeline folder."""
    os.system('rfigui '+vis)
    logger.info('End run_rfigui')

##Run aoflagger. Mode = auto uses best fit strategy for e-MERLIN (credit J. Moldon), Mode=user uses custon straegy for each field
#def run_aoflagger(vis,mode):
#	if mode == 'user':
#		x = vishead(vis,mode='list',listitems='field')['field'][0]
#		os.system('touch pre-cal_flag_stats.txt')
#		y = []
#		for i in range(len(x)):
#			if os.path.isfile(x[i]+'.rfis')==False:
#				y=y+[x[i]]
#		if len(y) != 0:
#			for i in range(len(y)):
#				print 'Missing rfistrategy for: '+y[i]
#			print 'Please run step 3 again!'
#		else:
#			while True:
#				s = raw_input('All rfistrategys are there: Proceed?:\n')
#				if s == 'yes' or s == 'y':
#					ms.writehistory(message='eMER_CASA_Pipeline: AOFlag with user specified strategies:',msname=vis)
#					for i in range(len(x)):
#						ms.writehistory(message='eMER_CASA_Pipeline: Flagging field, '+x[i]+' with strategy: '+x[i]+'.rfis',msname=vis)
#						print 'Flagging field, '+x[i]+' with strategy: '+x[i]+'.rfis'
#						os.system('aoflagger -fields '+str(i)+' -strategy '+x[i]+'.rfis  '+vis+ '| tee -a pre-cal_flag_stats.txt')
#					break
#				if s == 'no' or s == 'n':
#					sys.exit('Please restart when you are happy')
#	elif mode == 'default':
#		print '---- Running AOflagger with eMERLIN default strategy ----\n'
#		os.system('aoflagger -strategy eMERLIN_default_ao_strategy_v1.rfis '+vis)
#		ms.writehistory(message='eMER_CASA_Pipeline: AOFlag with default strategy, complete',msname=vis)
#	else:
#		print 'Error: Please use either mode=user or mode=default'
#		sys.exit()

def run_aoflagger_fields(vis,fields='all', pipeline_path='./'):
    """This version of the autoflagger iterates through the ms within the mms structure selecting individual fields. It uses pre-defined strategies. The np.unique in the loop below is needed for single source files. The mms thinks there are many filds (one per mms). I think it is a bug from virtualconcat."""
    logger.info('Start run_aoflagger_fields')
    vis_fields = vishead(vis,mode='list',listitems='field')['field'][0]
    fields_num = {f:i for i,f in enumerate(vis_fields)}
    if fields == 'all':
        fields = vis_fields
    else:
        fields = np.atleast_1d(fields)
    old_aoflagger = check_aoflagger_version()
    for field in np.unique(fields):
        # First, check if user has a new strategy for this field in the local folder.
        # If not, check if user has produced a new strategy for this field in the pipeline folder (for typical sources, etc).
        # If not, check for default strategies for this field
        # If nothing is found, just use the default strategy
        if os.path.isfile('./aoflagger_strategies/user/{0}.rfis'.format(field))==True:
            aostrategy = './aoflagger_strategies/user/{0}.rfis'.format(field)
        elif os.path.isfile(pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field))==True:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field)
        elif os.path.isfile(pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field))==True:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field)
        else:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format('default_faint')
        logger.info('Running AOFLagger for field {0} ({1}) using strategy {2}'.format(field,fields_num[field], aostrategy))
        if old_aoflagger: # < 2.9
            flagcommand = 'time aoflagger -strategy {0} {1}'.format(aostrategy, vis+'/SUBMSS/*{0}.mms.*.ms'.format(field))
        else: # >= 2.9
            flagcommand = 'time aoflagger -fields {2} -strategy {0} {1}'.format(aostrategy, vis, fields_num[field])
        os.system(flagcommand+' | tee -a pre-cal_flag_stats.txt')
        ms.writehistory(message='eMER_CASA_Pipeline: AOFlag field {0} with strategy {1}:'.format(field, aostrategy),msname=vis)
    logger.info('End run_aoflagger_fields')

def check_aoflagger_version():
	logger.info('Checking AOflagger version')
	from subprocess import Popen, PIPE
	process = Popen(['aoflagger'], stdout=PIPE)
	(output, err) = process.communicate()
	exit_code = process.wait()
	version = output.split()[1]
	version_list = version.split('.')
	if (version_list[0] == '2') and (int(version_list[1]) < 9):
		old_aoflagger = True
	else:
		old_aoflagger = False
	logger.info('AOflagger version is {0}'.format(version))
	return old_aoflagger

def ms2mms(vis,mode):
	logger.info('Start ms2mms')
	if mode == 'parallel':
		partition(vis=vis,outputvis=vis[:-3]+'.mms',createmms=True,separationaxis="baseline",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)
		if os.path.isdir(vis[:-3]+'.mms') == True:
			os.system('rm -r '+vis)
			os.system('rm -r '+vis+'.flagversions')
		ms.writehistory(message='eMER_CASA_Pipeline: Converted MS to MMS for parallelisation',msname=vis[:-3]+'.mms')

	## Need to use single if you need to aoflag the data later
	if mode == 'single':
		partition(vis=vis,outputvis=vis[:-3]+'.ms',createmms=False,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)
		if os.path.isdir(vis[:-3]+'.ms') == True:
			os.system('rm -r '+vis)
			os.system('rm -r '+vis+'.flagversions')
	logger.info('End ms2mms')

def ms2mms_fields(msfile):
    logger.info('Start ms2mms_fields')
    output_mmsfile = msfile[:-3]+'.mms'
    fields = vishead(msfile, mode = 'list', listitems = 'field')['field'][0]
    mmsfiles = []
    for field in fields:
        logger.info('Running partition on field found: {}'.format(field))
        mmsfile = msfile[:-3]+'_'+field+'.mms'
        mmsfiles.append(mmsfile)
        partition(vis=msfile, outputvis=mmsfile, createmms=True, separationaxis="baseline", numsubms="auto", flagbackup=False, datacolumn="all", field= field, spw="", scan="", antenna="", correlation="", timerange="", intent="", array="", uvrange="", observation="", feed="", disableparallel=None, ddistart=None, taql=None)
    # Virtual concatenation. No data copied, just moved to SUBMMS directory
    if len(mmsfiles) == 1: # No need to concatenate because there is only one file
        os.system('mv {0} {1}'.format(mmsfiles[0], output_mmsfile))
    elif len(mmsfiles) > 1:
        virtualconcat(vis = mmsfiles, concatvis = output_mmsfile, copypointing=True)
    else:
        logger.critical('No MMS files found.')
    if os.path.isdir(msfile) == True:
        os.system('rm -r '+msfile)
        os.system('rm -r '+msfile+'.flagversions')
    logger.info('End ms2mms_fields')

def do_prediagnostics(vis,plot_dir):
	##Pre diagnostics for measurement sets##
	## Includes:
	## - Antenna positions
	## - Amplitude vs. Time
	## - Amplitude vs. Frequency
	## - Phase vs. Time
	## - Phase vs. Frequency
	## - Closures (if task is available)
	## - Listobs summary
	logger.info('Start prediagnostics')
	if os.path.isdir(plot_dir) == False:
		os.system('mkdir '+plot_dir)
	if os.path.isdir('./'+plot_dir+'pre-calibration') == False:
		os.system('mkdir ./'+plot_dir+'pre-calibration')
	directory = plot_dir
	## Get information from ms
	x = vishead(vis,mode='list',listitems='field')['field'][0]
	tb.open(vis+'/SPECTRAL_WINDOW')
	nChan = str(tb.getcol('NUM_CHAN')[0])
	time = str(10E6)
	for i in range(len(x)):
		## - Amplitude vs. Time
		plotms(vis=vis,xaxis='time',yaxis='amplitude',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel=str(nChan), iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time.pdf',highres=True ,dpi=1200,expformat='pdf',exprange='all', showgui=False)
		os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time_* '+directory+'Pre-cal_amp_vs_time_'+vis+'_'+x[i]+'.pdf')
		os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time_*')
		## - Amplitude vs frequency
		plotms(vis=vis,xaxis='frequency',yaxis='amplitude',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel='1', avgtime=time, iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency.pdf',highres=True ,dpi=1200,expformat='pdf',exprange='all', showgui=False)
		os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency_* '+directory+'Pre-cal_amp_vs_frequency_'+vis+'_'+x[i]+'.pdf')
		os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency_*')

		## - Phase vs time
		plotms(vis=vis,xaxis='time',yaxis='phase',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel=str(nChan), iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time.pdf', expformat='pdf',highres=True,dpi=1200,exprange='all', showgui=False)
		os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time_* '+directory+'Pre-cal_phase_vs_time_'+vis+'_'+x[i]+'.pdf')
		os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time_*')

		## - Phase vs frequency
		plotms(vis=vis,xaxis='frequency',yaxis='phase',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgtime=time, iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency.pdf', expformat='pdf',highres=True,dpi=1200,exprange='all', showgui=False)
		os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency_* '+directory+'Pre-cal_phase_vs_frequency_'+vis+'_'+x[i]+'.pdf')
		os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency_*')
	#vishead(vis=vis,listfile=directory+vis+'.listobs')
	#plotants(vis=vis,figfile=directory+vis+'.plotants.png')
	## Amplitude vs Time:
	logger.info('End prediagnostics')

def solve_delays(vis,caltable_name,calsources,solint,refant,caldir,plotdir,combine='',spw='',timerange='',minblperant=2,minsnr=2):
    logger.info('Start solve_delays')
    caltable = caldir+caltable_name
    caltableplot = plotdir+caltable_name+'.png'
    rmdir(caltable)
    rmfile(caltableplot)
    logger.info('Delay calibration in: {0}'.format(caltable))
    logger.info('Delay calibration plot in: {0}'.format(caltableplot))
    gaincal(vis=vis, gaintype='K', caltable=caltable, field=calsources, solint=solint, combine=combine, refant=refant, spw=spw, timerange=timerange, minblperant=minblperant, minsnr=minsnr)
    logger.info('caltable: {0}, figfile: {1}'.format(caltable, caltableplot))
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',showgui=False,figfile=caltableplot, fontsize = 8)
    logger.info('End solve_delays')
    return


def dfluxpy(freq,baseline):
	#######
	# Python version of 3C286 flux calculation program (original author unknown)
	# ..............................
	# Author DMF       20/10/2011
	# ..............................
	#
	# Update to use Perley & Butler 2012 coefficients
	# 10/04/2013
	# DMF
	########

	# Reworked to use the 1999 VLA flux formula, and a 2nd formula to give a continuous estimate of the resolved fraction, by Ian Stewart, JBO, 8 Aug 2007.
	# Minor changes by amsr, 8 Aug 2007

	# my $program_name = 'dflux'; # $0 returns the './' prefix if this is used.

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
