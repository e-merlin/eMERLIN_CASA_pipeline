from Tkinter import *
import pickle
class GUI_pipeline:
	def __init__(self):
		self.root = Toplevel()
		self.root.title('eMERLIN CASA Pipeline') 
		logo = PhotoImage(file='emerlin-2.gif')
		self.w1 = Label(self.root,image=logo)
		explanation = """This is the eMERLIN pipeline"""
		self.w2 = Label(self.root,justify='left',padx=10,text=explanation)

		self.w1.grid(row=0,column=2,columnspan=2,rowspan=2,sticky='w,e,n,s')
		self.w2.grid(row=0,column=0,columnspan=1,rowspan=2)
		self.w4.grid(row=3,column=0)

		self.inbase = StringVar()
		self.inbase_label = Label(self.root,text='file name',padx=5,justify='right')
		self.inbase_entry = Entry(self.root,textvariable=self.inbase)
		self.inbase_entry.grid(row=3,column=1,columnspan=2)

		self.phscals = StringVar()
		self.w5 = Entry(self.root,textvariable=self.phscals)
		self.w5.grid(row=4,column=1,columnspan=2)

		self.w6 = Button(self.root,text='Confirm?',command=self.printmsg)
		### Quit button ###
		self.w3 = Button(self.root,text='Quit',command=self.root.quit)
		self.w3.grid(row=5,column=2)
		###################
		### Destroy GUI ###
		self.root.mainloop()
		self.root.destroy()
		###################

	def printmsg(self):
		x = self.phscals.get()
		print x
		pickle.dump(x,open("save.p",'wb'))


GUI_pipeline().printmsg()
