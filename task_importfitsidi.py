import os
import shutil
from taskinit import casalog, mstool, tbtool, write_history
import numpy as np

def importfitsidi(fitsidifile,vis,constobsid=None,scanreindexgap_s=None):
	"""Convert FITS-IDI visibility file into a CASA visibility file (MS).

	Keyword arguments:
	fitsidifile -- Name(s) of input FITS IDI file(s)
		default: None; example='3C273XC1.IDI' or ['3C273XC1.IDI1', '3C273XC1.IDI2']
	vis -- Name of output visibility file (MS)
		default: None; example: vis='3C273XC1.ms'
		
	constobsid -- If True a constant obs id == 0  of is given to all input files 
	        default = False (new obs id for each input file)

	scanreindexgap_s --  if > 0., a new scan is started whenever the gap between two
                integrations is > the given value (seconds) or when a new field starts
                default = 0. (no reindexing)
	"""

	#Python script
        retval = True
	try:
		casalog.origin('importfitsidi')
		casalog.post("")
                myms = mstool()
		mytb = tbtool()
		if(type(fitsidifile)==str):
			casalog.post('### Reading file '+fitsidifile, 'INFO')
			myms.fromfitsidi(vis,fitsidifile)
			myms.close()
		elif(type(fitsidifile)==list):
			clist = fitsidifile
			casalog.post('### Reading file '+clist[0], 'INFO')
			myms.fromfitsidi(vis,clist[0])
			myms.close()
			clist.pop(0)
			tname = '_importfitsidi_tmp_'+vis
			shutil.rmtree(tname, ignore_errors=True)
			for fidifile in clist:
				casalog.post('### Reading file '+fidifile, 'INFO')
				myms.fromfitsidi(tname,fidifile)
				myms.close()
				myms.open(vis, nomodify=False)
				myms.concatenate(msfile=tname, freqtol='', dirtol='')
				myms.close()
				shutil.rmtree(tname, ignore_errors=True)
		else:
                        raise Exception, 'Parameter fitsidifile should be of type str or list'			

		if (constobsid):
			mytb.open(vis+'/OBSERVATION', nomodify=False)
			nobs = mytb.nrows()
			cando = True
			if nobs>1:
				casalog.post('Trying to keep obsid constant == 0 for all input files', 'INFO')
				# check if all observations are from the same telescope; if not warn and leave as is
				tels = mytb.getcol('TELESCOPE_NAME')
				for i in range(1,nobs):
					if tels[i]!=tels[0]:
						cando = False

				if cando:
					# get min and max time and write them into the first row;
					casalog.post('Adjusting OBSERVATION table', 'INFO')
					timeranges = mytb.getcol('TIME_RANGE')
					ttr = timeranges.transpose()
					newmin = min(ttr[0])
					newmax = max(ttr[1])
					mytb.putcell('TIME_RANGE', 0, [newmin,newmax])
					# delete the other rows
					mytb.removerows(range(1,nobs))
				else:
					casalog.post('The input files stem from different telescopes. Need to give different obs id.', 'WARN')
			mytb.close()
			
			if cando:
				# give the same obs id == 0 to the entire output MS
				casalog.post('Setting observation ID of all integrations to 0', 'INFO')
				mytb.open(vis, nomodify=False)
				for i in xrange(0, mytb.nrows()):
					mytb.putcell('OBSERVATION_ID', i, 0)
				mytb.close()


		else: # don't want constant obs id
			if(type(fitsidifile)==list and len(fitsidifile)>1):
				casalog.post('Incrementing observation ID for each input file ...', 'INFO')
			
		if (scanreindexgap_s > 0.):
			# reindex the scan column
			mytb.open(vis, nomodify=False)
			times = mytb.getcol('TIME')
			fields = mytb.getcol('FIELD_ID')
			arrayids = mytb.getcol('ARRAY_ID')
			scannumbers = mytb.getcol('SCAN_NUMBER')

			timesorted = np.argsort(np.array(times)) 

			scannumber = 1
			prevtime = times[timesorted[0]]
			prevfield = fields[timesorted[0]]
			prevarrayid = arrayids[timesorted[0]]
			scannumbers[timesorted[0]] = scannumber

			for i in xrange(1,mytb.nrows()):
				ii = timesorted[i]
				timenow = times[ii]
				fieldnow = fields[ii]
				arrayidnow = arrayids[ii]
				if (timenow-prevtime > scanreindexgap_s) \
					    or (fieldnow != prevfield) \
					    or (arrayidnow != prevarrayid):
					scannumber += 1
					casalog.post("Starting new scan "+str(scannumber)+" at "+str(timenow)\
							     +", field "+str(fieldnow)+", array_id "+str(arrayidnow), 'INFO')
				scannumbers[ii] = scannumber
				prevtime = timenow
				prevfield = fieldnow
				prevarrayid = arrayidnow

			mytb.putcol('SCAN_NUMBER', scannumbers)	
			mytb.close()
		
	        # write history
                try:
                        param_names = importfitsidi.func_code.co_varnames[:importfitsidi.func_code.co_argcount]
                        param_vals = [eval(p) for p in param_names]   
                        retval &= write_history(myms, vis, 'importfitsidi', param_names,
                                                param_vals, casalog)

                except Exception, instance:
                        casalog.post("*** Error \'%s\' updating HISTORY" % (instance),
                                     'WARN')

	except Exception, instance: 
		print '*** Error ***',instance
		shutil.rmtree('_importfitsidi_tmp_'+vis, ignore_errors=True)
		raise Exception, instance


