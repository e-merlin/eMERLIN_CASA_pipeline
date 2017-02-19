#!/usr/local/python
import os
from casa import *
from casa import table as tb
from casa import ms
import Tkinter,tkFileDialog
##imports fitsfile to ms

def GUI_pipeline(vis):
	def get_var(w):
		print Tkinter.Entry.get(w)
	root = Tkinter.Toplevel()
	logo = Tkinter.PhotoImage(file='emerlin-2.gif')
	w1 = Tkinter.Label(root,image=logo)
	explanation = """This is the eMERLIN pipeline"""
	w2 = Tkinter.Label(root,justify='left',padx=10,text=explanation)
	w3 = Tkinter.Button(root,text='Quit',command=root.quit)
	w4 = Tkinter.Label(root,text='targets',padx=5,justify='right')
	w5 = Tkinter.Entry(root)
	w6 = Tkinter.Button(root,text='Confirm?',command=get_var(w5))
	w1.grid(row=0,column=1,columnspan=2,rowspan=2,sticky='w,e,n,s')
	w2.grid(row=0,column=0,columnspan=1,rowspan=2)
	w3.grid(row=2,column=2)
	w4.grid(row=3,column=0)
	w5.grid(row=3,column=1,columnspan=2)
	w6.grid(row=3,column=3,columnspan=1)
	root.mainloop()
	root.destroy()
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
	print 'You have been transformed from an ugly UVFITS to beautiful MS'
	return

##Hanning smoothing and flag of autocorrelations, will delete original and rename
def hanningflag(inputvis,deloriginal):
	os.system('rm -r '+inputvis+'_hanning.ms')
	hanningsmooth(vis=inputvis,outputvis=inputvis+'_hanning.ms',datacolumn='data')
	flagdata(vis=inputvis+'_hanning.ms',mode='manual',autocorr=True)
	if deloriginal==True:
		os.system('rm -r '+inputvis)
		os.system('mv '+inputvis+'_hanning.ms '+inputvis)
		os.system('mv '+inputvis+'_hanning.ms.flagversions '+inputvis+'.flagversions')
		ms.writehistory(message='eMER_CASA_Pipeline: Hanning smoothed data, complete',msname=inputvis)
	else:
		print 'Original not deleted, '+inputvis+'_hanning.ms is the new measurement set'
		ms.writehistory(message='eMER_CASA_Pipeline: Hanning smoothed data, complete',msname=inputvis+'_hanning.ms')
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
					for i in range(len(x)):
						print 'Flagging field, '+x[i]+' with strategy: '+x[i]+'.rfis'
						os.system('aoflagger -fields '+str(i)+' -strategy '+x[i]+'.rfis  '+vis+ '| tee -a pre-cal_flag_stats.txt')
					break
				if s == 'no' or s == 'n':
					sys.exit('Please restart when you are happy')
	elif mode == 'default':
		print '---- Running AOflagger with eMERLIN default strategy ----\n'
		os.system('aoflagger -strategy eMERLIN_default_ao_strategy_v1.rfis '+vis)
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

	## Need to use single if you need to aoflag the data later	
	if mode == 'single':
		partition(vis=vis,outputvis=vis[:-3]+'.ms',createmms=False,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)
		if os.path.isdir(vis[:-3]+'.ms') == True:
			os.system('rm -r '+vis)
			os.system('rm -r '+vis+'.flagversions')

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



