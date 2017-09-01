## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
import numpy as np
import pickle
from casa import table as tb
from casa import ms
from Tkinter import *
import getopt
import logging

# Find path of pipeline to find external files (like aoflagger strategies or emerlin-2.gif)
pipeline_path = os.path.dirname(sys.argv[np.where(np.asarray(sys.argv)=='-c')[0][0] + 1]) + '/'
sys.path.append(pipeline_path)
import functions.eMERLIN_CASA_functions as em
import functions.eMERLIN_CASA_GUI as emGUI

# Setup logger
logger = logging.getLogger('logger')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('eMCP.log', mode = 'a') # create a file handler
handler.setLevel(logging.INFO)
formatter = logging.Formatter(fmt='%(asctime)s | %(levelname)s | %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler) # add the handlers to the logger
consoleHandler = logging.StreamHandler() # Stream errors to terminal also
consoleHandler.setFormatter(formatter)
logger.addHandler(consoleHandler)

logger.info('Starting pipeline')
logger.info('Running pipeline from: {}'.format(pipeline_path))

##Inputs##
inputs = em.check_in(pipeline_path)
data_dir = em.backslash_check(inputs['data_dir'])
plots_dir = em.backslash_check(inputs['plots_dir'])
calib_dir = em.backslash_check(inputs['calib_dir'])
logger.info('Inputs used: {}'.format(inputs))

# Functions to save and load dictionaries
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


# Sources to process
sources = {}
sources['targets'] = inputs['targets'].replace(' ','')
sources['phscals'] = inputs['phscals'].replace(' ','')
sources['fluxcal'] = inputs['fluxcal']
sources['bpcal']   = inputs['bpcal']
sources['ptcal']   = inputs['ptcal']
sources['calsources'] = ','.join(set([sources['phscals'], sources['fluxcal'],
                                  sources['bpcal'], sources['ptcal']]))
sources['maincal'] = ','.join([sources['fluxcal'],sources['bpcal'], sources['ptcal']])
sources['allsources'] = sources['calsources'] +','+ sources['targets']
sources['no_fluxcal'] = sources['allsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
sources['cals_no_fluxcal'] = sources['calsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
sources['targets_phscals'] = sources['targets']+','+sources['phscals']

logger.info('Targets:   {0}'.format(sources['targets']))
logger.info('Phasecals: {0}'.format(sources['phscals']))
logger.info('Fluxcal:   {0}'.format(sources['fluxcal']))
logger.info('Bandpass:  {0}'.format(sources['bpcal']))
logger.info('Pointcal:  {0}'.format(sources['ptcal']))

# Flags applied to the data by the pipeline
try:
    flags = load_obj('./flags')
    logger.info('Loaded previous flags list from: {0}'.format('./flags.pkl'))
except:
    flags = []
    logger.info('Generating empty flags list')


## Create directory structure ##
em.makedir(plots_dir)
em.makedir(calib_dir)


if inputs['quit'] == 1: #Check from GUI if quit is needed
    logger.debug('Pipeline exit')
    sys.exit()

fitsfile = inputs['inbase']+'.fits'
msfile = inputs['inbase']+'.ms'

save_obj(sources, 'sources')

## Check for measurement sets in current directory otherwise drag from defined data directory
if os.path.isdir(inputs['inbase']+'.ms') == False and os.path.isdir(inputs['inbase']+'.mms') == False:
	if os.path.isdir(data_dir+inputs['inbase']+'.mms') == True:
		os.system('rsync -ar --progress {0} ./'.format(data_dir+inputs['inbase']+'.mms'))
	elif os.path.isdir(data_dir+inputs['inbase']+'.ms') == True:
		os.system('rsync -ar --progress {0} ./'.format(data_dir+inputs['inbase']+'.ms'))
	else:
		logger.info('No measurement set found, assuming you need to importfits or change data dir')
else:
	logger.info('Measurement set found: {}. Continuing with your inputs'.format(msfile))

## Pipeline processes, inputs are read from the inputs dictionary
if inputs['run_importfits'] == 1:
    em.run_importfitsIDI(data_dir,msfile)
    em.check_mixed_mode(msfile,mode='split')


### Retrieve MS information
# Antenna list and reference antenna
ms.open(msfile)
d = ms.getdata(['axis_info'],ifraxis=True)
ms.close()
antennas = np.unique('-'.join(d['axis_info']['ifr_axis']['ifr_name']).split('-'))

logger.info('Antennas in MS: {0}'.format(antennas))
refant = inputs['refant']
logger.info('Refant: {}'.format(refant))
if refant not in antennas:
    logger.warning('Selected reference antenna {0} not in MS!'.format(refant))


if inputs['hanning'] == 1:
	em.hanning(inputvis=msfile,deloriginal=True)

if inputs['rfigui'] == 1:
    em.run_rfigui(msfile)

### Convert MS to MMS ###
if inputs['ms2mms'] == 1:
    if em.check_aoflagger_version():
        em.ms2mms_fields(msfile=msfile)
    else:
        em.ms2mms(vis=msfile,mode='parallel')

## check for parallelisation
if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
    msfile = inputs['inbase']+'.mms'

if inputs['flagdata0_aoflagger'] == 1:
    flags = em.run_aoflagger_fields(vis=msfile, flags=flags, fields='all', pipeline_path = pipeline_path)

### Produce some initial plots ###
if inputs['do_prediag'] == 1:
	em.do_prediagnostics(msfile,plots_dir)


### A-priori flagdata: Lo&Mk2, edge channels, standard quack
if inputs['flagdata1_apriori'] == 1:
    flags = em.flagdata1_apriori(msfile=msfile, sources=sources, flags=flags,
                                 antennas=antennas, do_quack=True)


### Average data ###
if inputs['do_average1'] == 1:
    em.run_split(msfile, fields=sources['allsources'], width=4, timebin='2s')

# Check if averaged data already generated
if os.path.isdir('./'+inputs['inbase']+'_avg.mms') == True:
    msfile = './'+inputs['inbase']+'_avg.mms'
elif os.path.isdir('./'+inputs['inbase']+'_avg.ms') == True:
    msfile = './'+inputs['inbase']+'_avg.ms'
else:
    pass

logger.info('Using MS file: {0}'.format(msfile))



###################
### CALIBRATION ###
###################

# All the calibration steps will be saved in the dictionary caltables.pkl
# located in the calib directory. If it does not exist a new one is created.
try:
    caltables = load_obj(calib_dir+'caltables')
    logger.info('Loaded previous calibration tables from: {0}'.format(calib_dir+'caltables.pkl'))
except:
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    caltables = {}
    caltables['inbase'] = inputs['inbase']
    caltables['plots_dir'] = plots_dir
    caltables['calib_dir'] = calib_dir
    caltables['num_spw'] = num_spw
    caltables['refant'] = refant
    logger.info('New caltables dictionary created. Saved to: {0}'.format(calib_dir+'caltables.pkl'))

### Initialize models ###
if inputs['do_initialize_models'] == 1:  # Need to add parameter to GUI
    models_path = pipeline_path+'calibrator_models/'
    em.run_initialize_models(msfile=msfile, fluxcal=sources['fluxcal'],
                             models_path=models_path,
                             delmod_sources=sources['no_fluxcal'])


### Initial BandPass calibration ###
if inputs['do_initial_bandpass'] > 0:
    caltables = em.initial_bp_cal(msfile=msfile, caltables=caltables,
                                  previous_cal=[], bpcal=sources['bpcal'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_bandpass0')
    if inputs['do_initial_bandpass'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0'],
           previous_cal_targets=['bpcal.B0'])

### Flagdata using TFCROP and bandpass shape B0
if inputs['flagdata2_tfcropBP'] == 1:
    # If B0 has not been applied before, do it now
    if inputs['do_initial_bandpass'] != 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0'],
           previous_cal_targets=['bpcal.B0'])
    flags = em.flagdata2_tfcropBP(msfile=msfile, sources=sources, flags=flags)


### Delay calibration ###
if inputs['do_delay'] > 0:
    caltables = em.solve_delays(msfile=msfile, caltables=caltables,
                                previous_cal=['bpcal.B0'], calsources=sources['calsources'])
    # Should the previous_cal be bpcal.B0? Probably better delay fit, but later
    # delay.K1 is applied without bpcal.B0, when bpcal_sp.B1 is computed
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_delay')
    if inputs['do_delay'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0','delay.K1'],
           previous_cal_targets=['bpcal.B0','delay.K1'])


### Gain calibration ###
if inputs['do_gain_calibration'] > 0:
    caltables = em.initial_gaincal(msfile=msfile, caltables=caltables,
                                  previous_cal=['delay.K1', 'bpcal.B0'],
                                  calsources=sources['calsources'], phscals=sources['phscals'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_gaincal')
    if inputs['do_gain_calibration'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1','bpcal.B0'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1','bpcal.B0'])

### Flux scale ###
if inputs['do_fluxscale'] > 0:
    caltables = em.eM_fluxscale(msfile=msfile, caltables=caltables,
                                sources=sources,
                                ampcal_table='allcal_ap.G1', antennas=antennas)
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_fluxscale')
    if inputs['do_fluxscale'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal.B0'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal.B0'])

### BandPass calibration ###
if inputs['do_bandpass_sp'] > 0:
    caltables = em.bandpass_sp(msfile=msfile, caltables=caltables,
                               previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled'],
                               bpcal=sources['bpcal'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_bandpass_sp')
    if inputs['do_bandpass_sp'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal_sp.B1'])

### Apply calibration  ###
if inputs['do_applycal_all'] > 0:
    em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
       previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1'],
       previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal_sp.B1'])



logger.info('Pipeline finished')
logger.info('#################')
