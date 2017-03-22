#!/usr/local/python
import os
from casa import *
from casa import table as tb
from casa import ms
from Tkinter import *
import tkMessageBox

def check_in():
	inputs, processes = GUI_pipeline().confirm_parameters()


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


def run_importuvfits(fitsfile,vis):
	os.system('rm -r '+vis)
	importuvfits(fitsfile=fitsfile,vis=vis)
	ms.writehistory(message='eMER_CASA_Pipeline: Import uvfits to ms, complete',msname=vis)
    flagdata(vis=vis,mode='manual',autocorr=True)
	print 'You have been transformed from an ugly UVFITS to beautiful MS'
	return

##Hanning smoothing and flag of autocorrelations, will delete original and rename
def hanning(inputvis,deloriginal):
    if inputvis[-3:].upper() == '.MS':
        outputvis = inputvis[:-3]+'_hanning'+inputvis[-3:]
        os.system('rm -r '+outputvis)
        hanningsmooth(vis=inputvis,outputvis=outputvis,datacolumn='data')
    elif inputvis[-3:].upper() == 'MMS':
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
    return


##Run aoflagger. Mode = auto uses best fit strategy for e-MERLIN (credit J. Moldon), Mode=user uses custon straegy for each field
def run_aoflagger(vis,mode):

	if mode == 'user':
		x = vishead(vis,mode='list',listitems='field')['field'][0]
		os.system('touch pre-cal_flag_stats.txt')
		y = []
		for i in range(len(x)):
			if os.path.isfile(x[i]+'.rfis')==False:
				y=y+[x[i]]
		if len(y) != 0:
			for i in range(len(y)):
				print 'Missing rfistrategy for: '+y[i]
			print 'Please run step 3 again!'
		else:
			while True:
				s = raw_input('All rfistrategys are there: Proceed?:\n')
				if s == 'yes' or s == 'y':
					ms.writehistory(message='eMER_CASA_Pipeline: AOFlag with user specified strategies:',msname=vis)
					for i in range(len(x)):
						ms.writehistory(message='eMER_CASA_Pipeline: Flagging field, '+x[i]+' with strategy: '+x[i]+'.rfis',msname=vis)
						print 'Flagging field, '+x[i]+' with strategy: '+x[i]+'.rfis'
						os.system('aoflagger -fields '+str(i)+' -strategy '+x[i]+'.rfis  '+vis+ '| tee -a pre-cal_flag_stats.txt')
					break
				if s == 'no' or s == 'n':
					sys.exit('Please restart when you are happy')
	elif mode == 'default':
		print '---- Running AOflagger with eMERLIN default strategy ----\n'
		os.system('aoflagger -strategy eMERLIN_default_ao_strategy_v1.rfis '+vis)
		ms.writehistory(message='eMER_CASA_Pipeline: AOFlag with default strategy, complete',msname=vis)
	else:
		print 'Error: Please use either mode=user or mode=default'
		sys.exit()


def ms2mms(vis,mode):
	if mode == 'parallel':
		partition(vis=vis,outputvis=vis[:-3]+'.mms',createmms=True,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
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

def do_prediagnostics(vis):
	##Pre diagnostics for measurement sets##
	## Includes:
	## - Antenna positions
	## - Amplitude vs. Time 
	## - Amplitude vs. Frequency
	## - Phase vs. Time
	## - Phase vs. Frequency
	## - Closures (if task is available)
	## - Listobs summary

	if os.path.isdir('./pipeline_plots') == False:
		os.system('mkdir ./pipeline_plots')
	if os.path.isdir('./pipeline_plots/pre-calibration') == False:
		os.system('mkdir ./pipeline_plots/pre-calibration')
	directory = './pipeline_plots/pre-calibration/'
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
