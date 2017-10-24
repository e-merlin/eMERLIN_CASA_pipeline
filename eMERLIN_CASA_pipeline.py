## v0.00001 of an eMERLIN CASA pipeline ##
##Dependencies##
import os,sys,math
import numpy as np
import pickle
from Tkinter import *
import getopt
import logging

# CASA imports
from taskinit import *
from tasks import *

pipeline_version = 'v0.5.1'

# Find path of pipeline to find external files (like aoflagger strategies or emerlin-2.gif)
try:
    pipeline_path = os.path.dirname(sys.argv[np.where(np.asarray(sys.argv)=='-c')[0][0] + 1]) + '/'
except:
    pass

sys.path.append(pipeline_path)
import functions.eMERLIN_CASA_functions as em
import functions.eMERLIN_CASA_GUI as emGUI

casalog.setlogfile('casa_eMCP.log')

# Functions
# Save and load dictionaries
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_pipeline_version(pipeline_path):
    headfile = pipeline_path + '.git/HEAD'
    branch = open(headfile, 'rb').readlines()[0].strip().split('/')[-1]
    commit = open(pipeline_path + '.git/refs/heads/'+branch, 'rb').readlines()[0].strip()
    short_commit = commit[:7]
    return branch, short_commit

def run_pipeline(inputs=None, inputs_path=''):
    # Setup logger
    logger = logging.getLogger('logger')
    logging.Formatter.converter = time.gmtime
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler('eMCP.log', mode = 'a') # create a file handler
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(fmt='%(asctime)s | %(levelname)s | %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler) # add the handlers to the logger
    consoleHandler = logging.StreamHandler() # Stream errors to terminal also
    consoleHandler.setFormatter(formatter)
    logger.addHandler(consoleHandler)

    try:
        branch, short_commit = get_pipeline_version(pipeline_path)
    except:
        branch, short_commit = 'unknown', 'unknown'
    logger.info('Starting pipeline')
    logger.info('Running pipeline from: {}'.format(pipeline_path))
    logger.info('Pipeline version: {}'.format(pipeline_version))
    logger.info('Using github branch: {}'.format(branch))
    logger.info('github last commit: {}'.format(short_commit))
    logger.info('This log uses UTC times')

    # Inputs
    if inputs_path == '': # Running pipeline
        inputs = em.check_in(pipeline_path)
    else: # Running pipeline from within CASA
        inputs = em.headless(inputs_path)

    # Paths to use
    data_dir = em.backslash_check(inputs['data_dir'])
    plots_dir = em.backslash_check(inputs['plots_dir'])
    calib_dir = em.backslash_check(inputs['calib_dir'])

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
    # It will try to load the msfile if it is there. If not, it will try to
    # remake it. If the MS file is missing, nothing is done (we ignore
    # unaveraged data set).
    if os.path.isfile(msfile+'.msinfo'):
        msinfo = load_obj(msfile+'.msinfo')
    elif os.path.isdir(msfile):
        msinfo = em.get_msinfo(msfile, inputs)
        msinfo['Lo_dropout_scans'] = inputs['Lo_dropout_scans']
        save_obj(msfile, msfile+'.msinfo')
    else:
        logger.info('No unaveraged data or msinfo found. Pre-processing will not work.')
        pass

    ### Run AOflagger
    if inputs['flag_0_aoflagger'] == 1:
        flags = em.run_aoflagger_fields(vis=msfile, flags=flags, fields='all', pipeline_path = pipeline_path)

    ### Produce some initial plots ###
    if inputs['prediag'] == 1:
        em.do_prediagnostics(msfile,plots_dir)

    ### A-priori flagdata: Lo&Mk2, edge channels, standard quack
    if inputs['flag_1_apriori'] == 1:
        flags = em.flagdata1_apriori(msfile=msfile, msinfo=msinfo, flags=flags, do_quack=True)

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

    ### Load or create dictionary with ms information.
    if os.path.isfile(msfile+'.msinfo'):
        msinfo = load_obj(msfile+'.msinfo')
    elif os.path.isdir(msfile):
        msinfo = em.get_msinfo(msfile, inputs)
        msinfo['Lo_dropout_scans'] = inputs['Lo_dropout_scans']
        save_obj(msfile, msfile+'.msinfo')
    else:
        logger.info('No unaveraged data or msinfo found. Pre-processing will not work.')
        pass


    ### Defining reference antenna
    msinfo['refant'] = em.define_refant(msfile, msinfo, inputs)
    save_obj(msfile, msfile+'.msinfo')


    ### Load manual flagging file
    if inputs['flag_2b_manual'] == 1:
        flags = em.flagdata2_manual(msfile=msfile, inpfile=inputs['manual_flags_b'], flags=flags)



    ###################
    ### CALIBRATION ###
    ###################

    # All the calibration steps will be saved in the dictionary caltables.pkl
    # located in the calib directory. If it does not exist a new one is created.
    try:
        caltables = load_obj(calib_dir+'caltables')
        logger.info('Loaded previous calibration tables from: {0}'.format(calib_dir+'caltables.pkl'))
    except:
        caltables = {}
        caltables['inbase'] = inputs['inbase']
        caltables['plots_dir'] = plots_dir
        caltables['calib_dir'] = calib_dir
        caltables['num_spw'] = msinfo['num_spw']
        logger.info('New caltables dictionary created. Saved to: {0}'.format(calib_dir+'caltables.pkl'))

    caltables['refant'] = msinfo['refant']
    save_obj(caltables, calib_dir+'caltables')

    ### Initialize models ###
    if inputs['init_models'] == 1:  # Need to add parameter to GUI
        models_path = pipeline_path+'calibrator_models/'
        em.run_initialize_models(msfile=msfile, fluxcal=msinfo['sources']['fluxcal'],
                                 models_path=models_path,
                                 delmod_sources=msinfo['sources']['no_fluxcal'])


    ### Initial BandPass calibration ###
    if inputs['bandpass_0'] > 0:
        caltables = em.initial_bp_cal(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                      previous_cal=[])
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_bandpass0')
        if inputs['bandpass_0'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables,
                            sources=msinfo['sources'], previous_cal=['bpcal.B0'],
                            previous_cal_targets=['bpcal.B0'])

    ### Flagdata using TFCROP and bandpass shape B0
    if inputs['flag_3_tfcropBP'] == 1:
        # If B0 has not been applied before, do it now
        if inputs['bandpass_0'] != 2:
            em.run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
               previous_cal=['bpcal.B0'], previous_cal_targets=['bpcal.B0'])
        flags = em.flagdata3_tfcropBP(msfile=msfile, msinfo=msinfo, flags=flags)


    ### Delay calibration ###
    use_fringefit = False
    if inputs['delay'] > 0:
        if not use_fringefit:
            caltables = em.solve_delays(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                        previous_cal=['bpcal.B0'])
        else:
            logger.info('Full fringe fit selected.')
            caltables = em.delay_fringefit(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                           previous_cal=['bpcal.B0'])
        # Should the previous_cal be bpcal.B0? Probably better delay fit, but later
        # delay.K1 is applied without bpcal.B0, when bpcal_sp.B1 is computed
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_delay')
        if inputs['delay'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
               previous_cal=['bpcal.B0','delay.K1'],
               previous_cal_targets=['bpcal.B0','delay.K1'])


    ### Initial gain calibration ###
    if inputs['gain_0_p_ap'] > 0:
        caltables = em.initial_gaincal(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                       previous_cal=['delay.K1', 'bpcal.B0'])
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_gaincal')
        if inputs['gain_0_p_ap'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
               previous_cal=['delay.K1','allcal_p.G0', 'allcal_p_jitter.G0', 'allcal_ap.G1','bpcal.B0'],
               previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1','bpcal.B0'])

    ### Flux scale ###
    if inputs['fluxscale'] > 0:
        caltables = em.eM_fluxscale(msfile=msfile, caltables=caltables,
                                    sources=msinfo['sources'],
                                    ampcal_table='allcal_ap.G1',
                                    antennas=msinfo['antennas'])
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_fluxscale')
        if inputs['fluxscale'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables,
                            sources=msinfo['sources'],
                            previous_cal=['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled','bpcal.B0'],
                            previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal.B0'])

    ### BandPass calibration with spectral index information ###
    if inputs['bandpass_1_sp'] > 0:
        caltables = em.bandpass_sp(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                   previous_cal=['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled'])
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_bandpass_sp')
        if inputs['bandpass_1_sp'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
               previous_cal=['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1'],
               previous_cal_targets=['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal_sp.B1'])

    ### Amplitude calibration including spectral information ###
    if inputs['gain_1_amp_sp'] > 0:
        caltables = em.sp_amp_gaincal(msfile=msfile, msinfo=msinfo, caltables=caltables,
                                      previous_cal=['delay.K1','allcal_p.G0','allcal_p_jitter.G0','bpcal_sp.B1'])
        save_obj(caltables, calib_dir+'caltables')
        save_obj(caltables, calib_dir+'caltables_gaincal')
        if inputs['gain_1_amp_sp'] == 2:
            em.run_applycal(msfile=msfile, caltables=caltables,
                            sources=msinfo['sources'],
                            previous_cal=['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G3'],
                            previous_cal_targets=['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap_scan.G3'])


    ### Apply calibration  ###
    if inputs['applycal_all'] > 0:
        em.run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
           previous_cal=['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G3'],
           previous_cal_targets=['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap_scan.G3'])
        msinfo['applycal_all'] = True
        save_obj(msfile, msfile+'.msinfo')


    ### RFLAG automatic flagging ###
    if inputs['flag_4_rflag'] == 1:
        try:
            msinfo['applycal_all'] == True
            flags = em.flagdata4_rflag(msfile=msfile, msinfo=msinfo, flags=flags)
        except:
            logger.warning('flag_4_rflag selected but applycal_all has not been run. RFLAG only works on calibrated data!')
            logger.warning('flag_4_rflag will not be executed.')

    ### Run monitoring for bright sources:
    try:
        if inputs['monitoring'] == 1:
            flags, caltables = em.monitoring(msfile=msfile, msinfo=msinfo,
                                             flags=flags, caltables=caltables,
                                             previous_cal=[''])
    except:
        pass


    logger.info('Pipeline finished')
    logger.info('#################')

    return inputs, caltables, msinfo



# The  __name__ == "__main__" does not work for CASA.
try:
    if run_in_casa == True:
        # Running the pipeline from inside CASA
        print('Pipeline initialized. To run the pipeline within CASA use:')
        print('inputs, caltables, msinfo = run_pipeline(inputs_path=<input file>)')
except:
    inputs = em.check_in(pipeline_path)
    run_pipeline(inputs=inputs)

