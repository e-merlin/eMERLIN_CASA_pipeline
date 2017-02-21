#!/usr/local/python
import os
from casa import *
from casa import table as tb
from casa import ms
from Tkinter import *
##imports fitsfile to ms

class GUI_pipeline:
	def __init__(self):
		self.root = Toplevel()
		self.root.title('eMERLIN CASA Pipeline')
		logo = PhotoImage(file='emerlin-2.gif')
		self.w1 = Label(self.root,image=logo)
		self.w1.grid(row=0,column=2,columnspan=3,rowspan=2,sticky='w,e,n,s')
		self.w2 = Label(self.root,justify='left',padx=10,text="This is the eMERLIN pipeline")
		self.w2.grid(row=0,column=0,columnspan=1,rowspan=2)
		self.subtitle = Label(self.root,justify='center',padx=10,pady=10,text="------ Inputs ------")
		self.subtitle.grid(row=3,column=0,columnspan=3,rowspan=1)


		##file names ###
		self.inbase = StringVar()
		################################################
		##set default to be the last instance with fits#
		################################################
		x=[]
		for file in os.listdir('./'):
				if file.endswith('.ms') or file.endswith('.mms') or file.endswith('.fits'):
					x=x+[file]
		if len(x) > 3:
			self.inbase.set('Multiple ms, click check history for cwd')
		else:
			if len(x)!=0:
				x=x[0]
				if x.endswith('.ms'):
					x=x[:-3]
				if x.endswith('.mms'):
					x=x[:-4]
				if x.endswith('.fits'):
					x=x[:-5]
				self.inbase.set(x)
		self.history = StringVar()
		self.inbase_label = Label(self.root,text='UV file (without .fits/.ms/.mms):',justify='right')
		self.inbase_label.grid(row=4,column=0,columnspan=1,sticky='e')
		self.inbase_entry = Entry(self.root,textvariable=self.inbase)
		self.inbase_entry.grid(row=4,column=1,columnspan=2)
		self.inbase_button_his = Button(self.root,text='Check history',command=self.printhistory)
		self.inbase_button_his.grid(row=4,column=3)


		##################################################
		##### Set defaults to be set from targets ########
		##################################################
		fields = []
		if len(x) != 0:
			if os.path.isdir('./'+x+'.mms') == True:
				fields = vishead(x+'.mms',mode='list',listitems='field')['field'][0]
			elif os.path.isdir('./'+x+'.ms') == True:
				fields = vishead(x+'.ms',mode='list',listitems='field')['field'][0]
		## Targets ###
		self.targets = StringVar()
		self.targets_label = Label(self.root,text='Targets:',justify='right')
		self.targets_label.grid(row=5,column=0,columnspan=1,sticky='e')
		self.targets_entry = Entry(self.root,textvariable=self.targets)
		self.targets_entry.grid(row=5,column=1,columnspan=2)

		## Phase calibrators ###
		self.phscals = StringVar()
		self.phscals_label = Label(self.root,text='Phase calibrators:',justify='right')
		self.phscals_label.grid(row=6,column=0,columnspan=1,sticky='e')
		self.phscals_entry = Entry(self.root,textvariable=self.phscals)
		self.phscals_entry.grid(row=6,column=1,columnspan=2)

		## Flux calibrators ###
		x = ''
		if len(fields) != 0:
			if '1407+284' in fields:
				x = x+'1407+284'
			if '1331+305' in fields:
				if len(x)!=0:
					x=x+','
				x=x+'1331+305'
		self.fluxcal = StringVar()
		self.fluxcal.set(x)
		self.fluxcal_label = Label(self.root,text='Flux calibrators:',justify='right')
		self.fluxcal_label.grid(row=7,column=0,columnspan=1,sticky='e')
		self.fluxcal_entry = Entry(self.root,textvariable=self.fluxcal)
		self.fluxcal_entry.grid(row=7,column=1,columnspan=2)

		## Bandpass calibrators ###
		x=''		
		if len(fields) != 0:
			if '1407+284' in fields:
				x = x+'1407+284'
		self.bpcal = StringVar()
		self.bpcal.set(x)
		self.bpcal_label = Label(self.root,text='Bandpass calibrators:',justify='right')
		self.bpcal_label.grid(row=8,column=0,columnspan=1,sticky='e')
		self.bpcal_entry = Entry(self.root,textvariable=self.bpcal)
		self.bpcal_entry.grid(row=8,column=1,columnspan=2)

		## Point calibrators ###
		x=''
		if len(fields) != 0:
			if '1407+284' in fields:
				x = x+'1407+284'
		self.ptcal = StringVar()
		self.ptcal.set(x)
		self.ptcal_label = Label(self.root,text='Point calibrators:',justify='right')
		self.ptcal_label.grid(row=9,column=0,columnspan=1,sticky='e')
		self.ptcal_entry = Entry(self.root,textvariable=self.ptcal)
		self.ptcal_entry.grid(row=9,column=1,columnspan=2)

		## Refant ###
		self.refant = StringVar()
		self.refant.set(x)
		self.refant_label = Label(self.root,text='Reference antennas:',justify='right')
		self.refant_label.grid(row=10,column=0,columnspan=1,sticky='e')
		self.refant_entry = Entry(self.root,textvariable=self.refant)
		self.refant_entry.grid(row=10,column=1,columnspan=2)

		#-----------------------------------------------------------------#
		##############Processes############################################
		#-----------------------------------------------------------------#
		self.subtitle2 = Label(self.root,justify='center',padx=10,pady=10,text="----- Processes ------")
		self.subtitle2.grid(row=11,column=0,columnspan=3,rowspan=1)

		## Convert fits to ms ##
		self.run_importuvfits = IntVar()
		self.run_importuvfits_check = Checkbutton(self.root,text='Convert fits to ms',variable =self.run_importuvfits,onvalue=1,offvalue=0,justify='left')
		self.run_importuvfits_check.grid(row=12,column=0,sticky='w')
		#---------------------##
		
		## Hanning smoothing and flag autocorrelations ##
		self.hanningflag = IntVar()
		self.hanningflag_check = Checkbutton(self.root,text='Hanning smoothing',variable =self.hanningflag,onvalue=1,offvalue=0)
		self.hanningflag_check.grid(row=13,column=0,sticky='w')
		# ---------------------------------------------##		
		
		# Auto flagging #
		self.autoflag = IntVar()
		self.autoflag_check = Checkbutton(self.root,text='Autoflagging?',variable =self.autoflag,onvalue=1,offvalue=0,pady=5)
		self.autoflag_check.grid(row=14,column=0,sticky='w')

		## Rfigui to set strategies ##
		self.rfigui = IntVar()
		self.rfigui_check = Checkbutton(self.root,text='Set strategies per source?',variable =self.rfigui,onvalue=1,offvalue=0)
		self.rfigui_check.grid(row=14,column=1,sticky='w')
		# -------------------------------------------##

		## Convert to mms ##		
		self.ms2mms = IntVar()
		self.ms2mms_check = Checkbutton(self.root,text='Convert to MMS',variable =self.ms2mms,onvalue=1,offvalue=0,pady=5)
		self.ms2mms_check.grid(row=15,column=0,sticky='w')
		## ---------------##

		self.do_prediag = IntVar()
		self.do_prediag_check = Checkbutton(self.root,text='Pre-diagnostics',variable =self.do_prediag,onvalue=1,offvalue=0,pady=5)
		self.do_prediag_check.grid(row=16,column=0,sticky='w')

		## Set parameters ##
		self.w6 = Button(self.root,text='Confirm?',command=self.confirm_parameters)
		self.w6.grid(row=100,column=2,sticky='e')
		
		### Run button ###
		self.run = Button(self.root,text='Run',command=self.root.quit)
		self.run.grid(row=100,column=1,sticky='e')
		self.quit_var = IntVar()
		self.quit = Button(self.root,text='Quit',command=self.quit)
		self.quit.grid(row=100,column=3,sticky='e')
		###################
		### Destroy GUI ###
		self.root.update_idletasks()
		self.root.mainloop()
		self.root.destroy()
		###################

	def quit(self):
		self.quit_var.set(1)
		self.root.quit()

	def confirm_parameters(self):
		self.inputs = {'quit':self.quit_var.get(),'inbase':self.inbase.get(),'targets':self.targets.get(),'phscals':self.phscals.get(),'fluxcal':self.fluxcal.get(),'bpcal':self.bpcal.get(),'ptcal':self.ptcal.get(),'refant':self.refant.get()}
		self.processes = {'run_importuvfits':self.run_importuvfits.get(),'hanningflag':self.hanningflag.get(),'autoflag':self.autoflag.get(),'rfigui':self.rfigui.get(),'ms2mms':self.ms2mms.get(),'do_prediag':self.do_prediag.get()}
		print self.inputs
		return self.inputs, self.processes

	def printhistory(self):
		def check_his(msname):
			tb.open(msname+'/HISTORY')
			x = tb.getcol('MESSAGE')
			y = [i for i, item in enumerate(x) if 'eMER_CASA_Pipeline:' in item]
			if len(y) == 0:
				print 'Measurement set: '+msname+' has not been processed'
			else:
				print 'Measurement set: '+msname+' has these steps conducted:'
			for i in range(len(y)):
					print x[y[i]]
		if os.path.isdir('./'+self.inbase_entry.get()+'.ms') == True:
			msname = self.inbase_entry.get()+'.ms'
			check_his(msname)
		elif os.path.isdir('./'+self.inbase_entry.get()+'.mms') == True:
			msname = self.inbase_entry.get()+'.mms'
			check_his(msname)
		else:
			print 'Data set: '+self.inbase_entry.get()+'.fits/.ms/.mms does not exist'
			print 'Current working directory:'
			for file in os.listdir('./'):
				if file.endswith('.ms') or file.endswith('.mms') or file.endswith('.fits'):
					print file

	def default_inbase(self):
		x = []
		for file in os.listdir('./'):
				if file.endswith('.ms') or file.endswith('.mms') or file.endswith('.fits'):
					x=x+[file]
		if len(x) > 3:
			return 'Multiple ms, click check history for cwd'
		else:
			x=x[0]
			if x.endswith('.ms'):
				return x[:-3]
			if x.endswith('.mms'):
				return x[:-4]
			if x.endswith('.fits'):
				return x[:-5]

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
