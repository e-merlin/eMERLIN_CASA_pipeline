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

def check_mixed_mode(vis,mode):
	logger.info('Check for mixed mode')
	tb.open(vis + '/SPECTRAL_WINDOW')
	bw_spw = np.array(tb.getcol('TOTAL_BANDWIDTH'))
	tb.close()
	if len(np.unique(bw_spw)) != 1:
		if mode == 'split':
			logger.info('Splitting continuum from spectral line')
			cont_spw = np.where(bw_spw==np.max(np.unique(bw_spw)))[0]
			print np.array2string(cont_spw, separator=',')[1:-1]
			split(vis=vis, outputvis=vis+'.continuum', spw=np.array2string(cont_spw, separator=',')[1:-1], datacolumn='data')
			spec_line = np.delete(bw_spw, cont_spw)
			logger.info('Splitting spectral line')
			for i in range(len(np.unique(spec_line))):
				spec_line_spw = np.where(bw_spw==np.unique(spec_line)[i])[0]
				split(vis=vis, outputvis=vis+'.sp{0}'.format(i), spw=np.array2string(spec_line_spw, separator=',')[1:-1],datacolumn='data')
				ms.writehistory(message='eMER_CASA_Pipeline: Spectral line split from {0}'.format(vis),msname=vis+'.sp{0}'.format(i))
			ms.writehistory(message='eMER_CASA_Pipeline: Spectral lines split from this ms',msname=vis)
			os.system('mv {0} {1}'.format(vis, vis+'.original'))
			os.system('mv {0} {1}'.format(vis+'.continuum', vis))
			logger.info('Will continue with continuum, original data is {0}'.format(vis+'.original'))
			return_variable = ''
		if mode == 'check':
			logger.info('MS is mixed mode. Please split')
			return_variable = True
	else:
		if mode == 'split':
			logger.info('Not mixed mode, continuing')
			return_variable = ''
		if mode == 'check':
			return_variable = False
	return return_variable

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

def solve_delays(msfile, inbase, calsources, solint, refant, caldir, plotdir, caltables, combine='', spw='', timerange='', minblperant=2, minsnr=2):
    """
    Input tables: None
    Output tables: caldir+inputs['inbase']+'_delay.K0'
    Output plots: plotdir+inputs['inbase']+'_delay.K0.png'
    """
    logger.info('Start solve_delays')
    delay_caltable = inbase +'_delay.K0'
    # This should be implemented in the main script
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    spwmap_out = [0]*num_spw
    caltable = caldir+delay_caltable
    caltableplot = plotdir+delay_caltable+'.png'
    gaintype = 'K'
#    rmdir(caltable)
#    rmfile(caltableplot)
    logger.info('Running gaincal on field {0}, gaintype = {1}, solint = {2}'.format(calsources, gaintype, solint))
    logger.info('Delay calibration in: {0}'.format(caltable))
    logger.info('Delay calibration plot in: {0}'.format(caltableplot))
#    gaincal(vis=msfile, gaintype=gaintype, caltable=caltable, field=calsources, solint=solint, combine=combine, refant=refant, spw=spw, timerange=timerange, minblperant=minblperant, minsnr=minsnr)
    logger.info('caltable: {0}, figfile: {1}'.format(caltable, caltableplot))
#    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',showgui=False,figfile=caltableplot, fontsize = 8)
    caltables['delay0'] = {}
    caltables['delay0']['tables']  = [caltable]
    caltables['delay0']['spwmaps'] = [spwmap_out]
    logger.info('End solve_delays')
    return caltables


def run_gaincal(msfile, caltable, calmode, solint, field, combine, refant, spw, previous_tables, previous_spwmap, caldir, plotdir, subplot, iteration, plotrange_phs = [-1,-1,-180,180], plotrange_amp = [-1,-1,-1,-1], timerange='', minblperant=2, minsnr=2):
    # This should be implemented in the main script
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    if 'spw' in combine.split():
        spwmap_out = [0]*num_spw
    else:
        spwmap_out = range(num_spw)
    basename = os.path.basename(caltable)
    caltableplot_phs = plotdir + basename +'_phs.png'
    caltableplot_amp = plotdir + basename +'_amp.png' 
    rmdir(caltable)
    rmfile(caltableplot_phs)
    rmfile(caltableplot_amp)
    logger.info('Running gaincal on field {0}, calmode = {1}, solint = {2}'.format(field, calmode, solint))
    logger.info('Previous calibration applied: {0}'.format(', '.join(previous_tables)))
    logger.info('Previous calibration spwmap: {0}'.format(previous_spwmap))
    logger.info('Generating calibration table: {0}'.format(caltable))
    gaincal(vis=msfile, calmode=calmode, field=field, caltable=caltable, solint=solint, combine=combine, refant=refant, spw=spw, timerange=timerange, gaintable=previous_tables, spwmap = previous_spwmap, minblperant=minblperant, minsnr=minsnr)
    logger.info('caltable: {0}, figfiles: {1}, {2}'.format(caltable, caltableplot_phs, caltableplot_amp))
    plotcal(caltable=caltable, xaxis='time', yaxis='phase', subplot=subplot, iteration=iteration, showgui=False, figfile=caltableplot_phs, fontsize = 8, plotrange = plotrange_phs)
    plotcal(caltable=caltable, xaxis='time', yaxis='amp', subplot=subplot, iteration=iteration, showgui=False, figfile=caltableplot_amp, fontsize = 8, plotrange = plotrange_amp)
    return caltable, spwmap_out

def run_bandpass(msfile, bptable, bpcal, refant, previous_tables, previous_spwmap, caldir, plotdir, spw='', solint='inf', combine='scan'):
    # This should be implemented in the main script
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    if 'spw' in combine.split():
        spwmap_out = [0]*num_spw
    else:
        spwmap_out = range(num_spw)
    basename = os.path.basename(bptable)
    bptableplot_phs = plotdir+basename+'_phs'+'.png'
    bptableplot_amp = plotdir+basename+'_amp'+'.png'
    rmdir(bptable)
    rmfile(bptableplot_phs)
    rmfile(bptableplot_amp)
    logger.info('Running bandpass on field {0}, solint = {1}, combine = {2}'.format(bpcal, solint, combine))
    logger.info('Previous calibration applied: {0}'.format(', '.join(previous_tables)))
    logger.info('Previous calibration spwmap: {0}'.format(previous_spwmap))
    logger.info('Generating bandpass table: {0}'.format(bptable))
    bandpass(vis=msfile, caltable=bptable, field=bpcal, fillgaps=16, solint=solint, combine=combine, solnorm=True, refant=refant, minblperant=2, gaintable=previous_tables, spwmap = previous_spwmap, minsnr=3)
    logger.info('bptable: {0}, figfile: {1}'.format(bptable, ', '.join([bptableplot_phs, bptableplot_amp])))
    plotcal(caltable=bptable, xaxis='freq', yaxis='phase', subplot=321,iteration='antenna', showgui=False, figfile=bptableplot_phs, fontsize = 8, plotrange = [-1,-1,-180,180])
    plotcal(caltable=bptable, xaxis='freq', yaxis='amp',  subplot=321, iteration='antenna', showgui=False, figfile=bptableplot_amp, fontsize = 8, plotrange = [-1,-1,-1,-1])
    return bptable, spwmap_out

def smooth_caltable(msfile, tablein, plotdir, caltable='', field='', smoothtype='median', smoothtime=120.):
    logger.info('Smoothing table: {0}, field {1}, smoothtype {2}, smoothtime {3}'.format(tablein, field, smoothtype, smoothtime))
    basename = os.path.basename(tablein)
    caltableplot_phs = plotdir + basename +'_phs.png'
    caltableplot_amp = plotdir + basename +'_amp.png'
    if caltable=='':
        os.system('mv {0} {1}'.format(caltableplot_phs, plotdir + basename +'_phs_pre_smooth.png'))
        os.system('mv {0} {1}'.format(caltableplot_amp, plotdir + basename +'_amp_pre_smooth.png'))
    smoothcal(vis=msfile, tablein=tablein, caltable=tablein+'smooth', field='', smoothtype='median', smoothtime=60*20.)
    logger.info('Pre-smoothing table saved to: {0}'.format(tablein+'_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein, tablein+'_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein+'smooth', tablein))
    plotcal(caltable=tablein, xaxis='time', yaxis='phase', subplot=321, iteration='antenna', showgui=False, figfile=caltableplot_phs, fontsize = 8, plotrange=[-1,-1,-180,180])
    plotcal(caltable=tablein, xaxis='time', yaxis='amp', subplot=321, iteration='antenna', showgui=False, figfile=caltableplot_amp, fontsize = 8, plotrange=[-1,-1,-1,-1])
    return

def initial_bp_cal(msfile, inbase, bpcal, refant, caldir, plotdir, caltables, minblperant=2, minsnr=2):
    """
    Input tables (caldir+inputs['inbase']+): '_delay.K0'
    Output tables (caldir+inputs['inbase']+): ['bpcal0.G0', 'bpcal0.G1', 'bpcal0.B0']
    Output plots (plotdir+inputs['inbase']+): ['bpcal0.G0.png', 'bpcal0.G1.png', 'bpcal0.B0.png']
    """
    logger.info('Start initial_bp_cal')
    # This should be implemented in the main script
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    previous_tables = caltables['delay0']['tables']
    previous_spwmap = caltables['delay0']['spwmaps']

    # 1 Phase calibration
    calmode1 = 'p'
    solint1 = '10s'
    caltable1 = caldir+inbase+'_bpcal0.G0'
    spwmap1 = range(num_spw)
#    run_gaincal(msfile=msfile, caltable=caltable1, calmode=calmode1, solint=solint1, field=bpcal, combine='', refant=refant, spw='', previous_tables=previous_tables, previous_spwmap=previous_spwmap, caldir=caldir, plotdir=plotdir, subplot=321, iteration='antenna', plotrange_phs = [-1,-1,-180,180])
    previous_tables.append(caltable1)
    previous_spwmap.append(spwmap1)

    # 2 A&P calibration
    calmode2 = 'ap'
    solint2 = '120s'
    caltable2 = caldir+inbase+'_bpcal0.G1'
    spwmap2 = range(num_spw)
#    run_gaincal(msfile=msfile, caltable=caltable2, calmode=calmode2, solint=solint2, field=bpcal, combine='', refant=refant, spw='', previous_tables=previous_tables, previous_spwmap=previous_spwmap, caldir=caldir, plotdir=plotdir, subplot=321, iteration='antenna', plotrange_phs = [-1,-1,-180,180])
    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    previous_tables.append(caltable2)
    previous_spwmap.append(spwmap2)

    # 3 Bandpass calibration
    bptable0 = caldir+inbase+'_bpcal0.B0'
    bptableplot0_phs = plotdir+'bpcal.'+bpcal+'_precal.B0_phs'+'.png'
    bptableplot0_amp = plotdir+'bpcal.'+bpcal+'_precal.B0.amp'+'.png'
#    run_bandpass(msfile=msfile, bptable=bptable0, bpcal=bpcal, refant=refant, previous_tables=previous_tables, previous_spwmap=previous_spwmap, caldir=caldir, plotdir=plotdir, spw='', solint='inf', combine='scan')
    logger.info('End initial_bp_cal')

    caltables['bpcal0'] = {}
    caltables['bpcal0']['tables']  = [caltable1, caltable2]
    caltables['bpcal0']['spwmaps'] = [spwmap1, spwmap2]
    logger.info('End solve_delays')

    return caltables



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
