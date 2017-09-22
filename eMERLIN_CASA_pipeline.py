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

casalog.setlogfile('casa_eMCP.log')

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


#################################
### LOAD AND PREPROCESS DATA  ###
#################################

## Pipeline processes, inputs are read from the inputs dictionary
if inputs['run_importfits'] == 1:
    em.run_importfitsIDI(data_dir,msfile)
    em.check_mixed_mode(msfile,mode='split')

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

### check for parallelisation
if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
    msfile = inputs['inbase']+'.mms'

### Create dictionary with ms information.
# It will try to remake the msinfo dictionary from the msfile. If the msfile is
# not there, it will try to load the pkl file. If nothing there (means that we
# don't care about the unaveraged data), nothing is done.
if os.path.isdir(msfile):
    msinfo = em.get_msinfo(msfile, inputs)
else:
    try:
        msinfo = load_ob(msfile)
    except:
        pass

### Run AOflagger
if inputs['flag_0_aoflagger'] == 1:
    flags = em.run_aoflagger_fields(vis=msfile, flags=flags, fields='all', pipeline_path = pipeline_path)

### Produce some initial plots ###
if inputs['prediag'] == 1:
	em.do_prediagnostics(msfile,plots_dir)

### A-priori flagdata: Lo&Mk2, edge channels, standard quack
if inputs['flag_1_apriori'] == 1:
    flags = em.flagdata1_apriori(msfile=msfile, sources=msinfo['sources'], flags=flags,
                                 antennas=msinfo['antennas'], do_quack=True)

### Load manual flagging file
if inputs['flag_2a_manual'] == 1:
    flags = em.flagdata2_manual(msfile=msfile, inpfile=inputs['manual_flags_a'], flags=flags)


### Average data ###
if inputs['average_1'] == 1:
    em.run_split(msfile, msinfo, width=4, timebin='2s')

# Check if averaged data already generated
if os.path.isdir('./'+inputs['inbase']+'_avg.mms') == True:
    msfile = './'+inputs['inbase']+'_avg.mms'
elif os.path.isdir('./'+inputs['inbase']+'_avg.ms') == True:
    msfile = './'+inputs['inbase']+'_avg.ms'
else:
    pass

logger.info('Using MS file: {0}'.format(msfile))


### Load manual flagging file
if inputs['flag_2b_manual'] == 1:
    flags = em.flagdata2_manual(msfile=msfile, inpfile=inputs['manual_flags_b'], flags=flags)


### Retrieve MS information
# Sources in the MS
sources['msfile_fields'] = ','.join(vishead(msfile,mode='list',listitems='field')['field'][0])
logger.info('Sources in MS {0}: {1}'.format(msfile, sources['msfile_fields']))

# Antenna list and reference antenna
ms.open(msfile)
d = ms.getdata(['axis_info'],ifraxis=True)
ms.close()
antennas = np.unique('-'.join(d['axis_info']['ifr_axis']['ifr_name']).split('-'))

logger.info('Antennas in MS: {0}'.format(antennas))

# Defining reference antenna
refant = inputs['refant']
refant_user = refant.replace(' ', '').split(',')
refant_in_ms = (np.array([ri in antennas for ri in refant_user])).all()

if not refant_in_ms:
    if refant != '':
        logger.warning('Selected reference antenna(s) {0} not in MS! User selection will be ignored'.format(refant))
    # Finding best antennas for refant
    refant, refant_pref = em.find_refant(msfile, field=sources['bpcal'],
                                         antennas='Mk2,Pi,Da,Kn', spws='2,3', scan='')
else:
    refant = ','.join(refant_user)

logger.info('Refant: {}'.format(refant))


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
if inputs['init_models'] == 1:  # Need to add parameter to GUI
    models_path = pipeline_path+'calibrator_models/'
    em.run_initialize_models(msfile=msfile, fluxcal=sources['fluxcal'],
                             models_path=models_path,
                             delmod_sources=sources['no_fluxcal'])


### Initial BandPass calibration ###
if inputs['bandpass_0'] > 0:
    caltables = em.initial_bp_cal(msfile=msfile, caltables=caltables,
                                  previous_cal=[], bpcal=sources['bpcal'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_bandpass0')
    if inputs['bandpass_0'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0'],
           previous_cal_targets=['bpcal.B0'])

### Flagdata using TFCROP and bandpass shape B0
if inputs['flag_3_tfcropBP'] == 1:
    # If B0 has not been applied before, do it now
    if inputs['bandpass_0'] != 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0'],
           previous_cal_targets=['bpcal.B0'])
    flags = em.flagdata3_tfcropBP(msfile=msfile, sources=sources, flags=flags)


### Delay calibration ###
if inputs['delay'] > 0:
    caltables = em.solve_delays(msfile=msfile, caltables=caltables,
                                previous_cal=['bpcal.B0'], calsources=sources['calsources'])
    # Should the previous_cal be bpcal.B0? Probably better delay fit, but later
    # delay.K1 is applied without bpcal.B0, when bpcal_sp.B1 is computed
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_delay')
    if inputs['delay'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['bpcal.B0','delay.K1'],
           previous_cal_targets=['bpcal.B0','delay.K1'])


### Initial gain calibration ###
if inputs['gain_0_p_ap'] > 0:
    caltables = em.initial_gaincal(msfile=msfile, caltables=caltables,
                                  previous_cal=['delay.K1', 'bpcal.B0'],
                                  calsources=sources['calsources'], phscals=sources['phscals'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_gaincal')
    if inputs['gain_0_p_ap'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1','bpcal.B0'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1','bpcal.B0'])

### Flux scale ###
if inputs['fluxscale'] > 0:
    caltables = em.eM_fluxscale(msfile=msfile, caltables=caltables,
                                sources=sources,
                                ampcal_table='allcal_ap.G1', antennas=antennas)
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_fluxscale')
    if inputs['fluxscale'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal.B0'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal.B0'])

### BandPass calibration with spectral index information ###
if inputs['bandpass_1_sp'] > 0:
    caltables = em.bandpass_sp(msfile=msfile, caltables=caltables,
                               previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled'],
                               bpcal=sources['bpcal'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_bandpass_sp')
    if inputs['bandpass_1_sp'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1'],
           previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal_sp.B1'])

### Amplitude calibration including spectral information ###
if inputs['gain_1_amp_sp'] > 0:
    caltables = em.sp_amp_gaincal(msfile=msfile, caltables=caltables,
                                  previous_cal=['delay.K1','allcal_p.G0','bpcal_sp.B1'],
                                  calsources=sources['calsources'])
    save_obj(caltables, calib_dir+'caltables')
    save_obj(caltables, calib_dir+'caltables_gaincal')
    if inputs['gain_1_amp_sp'] == 2:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
           previous_cal=['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_ap.G3'],
           previous_cal_targets=['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap.G3'])


### Apply calibration  ###
if inputs['applycal_all'] > 0:
    em.run_applycal(msfile=msfile, caltables=caltables, sources=sources,
       previous_cal=['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_ap.G3'],
       previous_cal_targets=['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap.G3'])


### Run monitoring for bright sources:
try:
    if inputs['monitoring'] == 1:
        flags, caltables = em.monitoring(msfile=msfile, data_dir=data_dir, sources=sources,
                      flags=flags, caltables=caltables, previous_cal=[''],
                      calsources=sources['calsources'], antennas=antennas)
except:
    pass


logger.info('Pipeline finished')
logger.info('#################')
