from Tkinter import *
from casa import *
import pickle
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
			return 'Multiple ms, click check history for cwd'
		else:
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
		if os.path.isdir('./'+x+'.mms') == True:
			fields = vishead(x+'.mms',mode='list',listitems='field')['field'][0]
			telescopes = vishead(x+'.mms',mode='list',listitems='telescope')['telescope'][0]
		elif os.path.isdir('./'+x+'.ms') == True:
			fields = vishead(x+'.ms',mode='list',listitems='field')['field'][0]
			telescopes = vishead(x+'.ms',mode='list',listitems='telescope')['telescope'][0]
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
		if '1407+284' in fields:
			x = x+'1407+284'
		self.ptcal = StringVar()
		self.ptcal.set(x)
		self.ptcal_label = Label(self.root,text='Point calibrators:',justify='right')
		self.ptcal_label.grid(row=9,column=0,columnspan=1,sticky='e')
		self.ptcal_entry = Entry(self.root,textvariable=self.ptcal)
		self.ptcal_entry.grid(row=9,column=1,columnspan=2)

		## Refant ###
		x=''
		if 'Mk2' in telescopes:
			x=x+'Mk2'
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
		## Convert fits to ms		
		self.doconvert = IntVar()
		self.doconvert_check = Checkbutton(self.root,text='Convert fits to ms',variable =self.doconvert,onvalue=1,offvalue=0,justify='left')
		self.doconvert_check.grid(row=12,column=0,sticky='w')

		self.dohanning = IntVar()
		self.dohanning_check = Checkbutton(self.root,text='Hanning smoothing',variable =self.dohanning,onvalue=1,offvalue=0)
		self.dohanning_check.grid(row=13,column=0,sticky='w')

		self.autoflag = IntVar()
		self.autoflag_check = Checkbutton(self.root,text='Autoflagging?',variable =self.autoflag,onvalue=1,offvalue=0,pady=10)
		self.autoflag_check.grid(row=14,column=0,sticky='w')

		self.rfigui = IntVar()
		self.rfigui_check = Checkbutton(self.root,text='Set strategies pre source?',variable =self.rfigui,onvalue=1,offvalue=0)
		self.rfigui_check.grid(row=14,column=1,sticky='w')

		self.w6 = Button(self.root,text='Confirm?',command=self.confirm_parameters)
		self.w6.grid(row=100,column=3)
		### Quit button ###
		
		self.w3 = Button(self.root,text='Quit',command=self.root.quit)
		self.w3.grid(row=100,column=2)
		###################
		### Destroy GUI ###
		self.root.update_idletasks()
		self.root.mainloop()
		self.root.destroy()
		###################

	def confirm_parameters(self):
		self.inputs = [self.inbase.get(),self.targets.get(),self.phscals.get(),self.fluxcal.get(),self.fluxcal.get(),self.bpcal.get(),self.ptcal.get(),self.refant.get()]
		print self.inputs
		return self.inputs
	def printhistory(self):
		if os.path.isdir('./'+self.inbase_entry.get()+'.ms') == True:
			msname = self.inbase_entry.get()+'.ms'
			if os.path.isdir('./'+self.inbase_entry.get()+'.mms') == True:
				msname = self.inbase_entry.get()+'.mms'
			tb.open(msname+'/HISTORY')
			x = tb.getcol('MESSAGE')
			y = [i for i, item in enumerate(x) if 'eMER_CASA_Pipeline:' in item]
			if len(y) == 0:
				print 'Measurement set: '+msname+' has not been processed'
			else:
				print 'Measurement set: '+msname+' has these steps conducted:'
			for i in range(len(y)):
				print x[y[i]]
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


GUI_pipeline().confirm_parameters()
