# Dependencies
import os,sys,math
import numpy as np
import pickle
from Tkinter import *
import getopt
import logging
import collections

# CASA imports
from taskinit import *
from tasks import *
import casadef


current_version = 'v0.9.02'

# Find path of pipeline to find external files (like aoflagger strategies or emerlin-2.gif)
try:
    pipeline_filename = sys.argv[sys.argv.index('-c') + 1]
    pipeline_path = os.path.abspath(os.path.dirname(pipeline_filename))
except:
    pass

if pipeline_path[-1] != '/':
    pipeline_path = pipeline_path + '/'
sys.path.append(pipeline_path)

import functions.eMERLIN_CASA_functions as em
import functions.weblog as emwlog
import functions.eMERLIN_CASA_plots as emplt
from default_params import defaults

casalog.setlogfile('casa_eMCP.log')

# Functions
# Save and load dictionaries
def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def get_pipeline_version(pipeline_path):
    headfile = pipeline_path + '.git/HEAD'
    branch = open(headfile, 'rb').readlines()[0].strip().split('/')[-1]
    commit = open(pipeline_path + '.git/refs/heads/'+branch, 'rb').readlines()[0].strip()
    short_commit = commit[:7]
    return branch, short_commit

def run_pipeline(inputs=None, inputs_path=''):
    # Paths to use
    weblog_dir = './weblog/'
    info_dir   = './weblog/info/'
    calib_dir  = './weblog/calib/'
    plots_dir  = './weblog/plots/'
    logs_dir   = './logs/'
    images_dir = './weblog/images/'

    ## Create directory structure ##
    em.makedir(weblog_dir)
    em.makedir(info_dir)
    em.makedir(plots_dir)
    em.makedir(calib_dir)
    em.makedir(images_dir)
    em.makedir(logs_dir)
    em.makedir(plots_dir+'caltables')

    pipeline_version = current_version

    # Continue with previous pipeline configuration if possible:
    try:
        eMCP = load_obj(info_dir + 'eMCP_info.pkl')
    except:
        eMCP = collections.OrderedDict()
        eMCP['steps'] = em.eMCP_info_start_steps()
        eMCP['is_mixed_mode'] = 'unknown'
        eMCP['img_stats'] = collections.OrderedDict()

    eMCP['inputs'] = inputs
    eMCP['defaults'] = defaults

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
    logger.info('CASA version: {}'.format(casadef.casa_version))
    logger.info('Pipeline version: {}'.format(pipeline_version))
    logger.info('Using github branch: {}'.format(branch))
    logger.info('github last commit: {}'.format(short_commit))
    logger.info('This log uses UTC times')
    eMCP['pipeline_path'] = pipeline_path
    eMCP['casa_version'] = casadef.casa_version
    em.check_pipeline_conflict(eMCP, pipeline_version)
    eMCP['pipeline_version'] = pipeline_version
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    # Inputs
    if inputs_path == '': # Running pipeline
        inputs = em.check_in(pipeline_path)
    else: # Running pipeline from within CASA
        inputs = em.headless(inputs_path)


    #################################
    ### LOAD AND PREPROCESS DATA  ###
    #################################

    ## Pipeline processes, inputs are read from the inputs dictionary
    if inputs['run_importfits'] > 0:
        eMCP = em.import_eMERLIN_fitsIDI(eMCP)

    if os.path.isdir('./'+inputs['inbase']+'.ms') == True:
        msfile = inputs['inbase']+'.ms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### check for parallelisation
    if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
        msfile = inputs['inbase']+'.mms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### Run AOflagger
    if inputs['flag_aoflagger'] > 0:
        eMCP = em.run_aoflagger_fields(eMCP)

    ### A-priori flagdata: Lo&Mk2, edge channels, standard quack
    if inputs['flag_apriori'] > 0:
        eMCP = em.flagdata1_apriori(eMCP)

    ### Load manual flagging file
    if inputs['flag_manual'] > 0:
        eMCP = em.flagdata_manual(eMCP, run_name='flag_manual')

    ### Multi phase center ###
    if inputs['shift_field_pos'] > 0:
        eMCP = em.shift_all_positions(eMCP)

    ### Average data ###
    if inputs['average'] > 0:
        eMCP = em.run_split(eMCP)

    # Check if averaged data already generated
    if os.path.isdir('./'+inputs['inbase']+'_avg.mms') == True:
        msfile = './'+inputs['inbase']+'_avg.mms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)
    elif os.path.isdir('./'+inputs['inbase']+'_avg.ms') == True:
        msfile = './'+inputs['inbase']+'_avg.ms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### Produce some plots ###
    if inputs['plot_data'] == 1:
        eMCP = emplt.make_4plots(eMCP, datacolumn='data')

    ### Save flag status up to this point
    if inputs['save_flags'] == 1:
        eMCP = em.saveflagstatus(eMCP)

    ###################
    ### CALIBRATION ###
    ###################

    # All the calibration steps will be saved in the dictionary caltables.pkl
    # located in the calib directory. If it does not exist a new one is created.
    any_calsteps = ['bandpass', 'delay', 'flag_tfcropBP','gain_p_ap','fluxscale','bandpass_sp','gain_amp_sp','applycal_all']
    if np.array([inputs[cal]>0 for cal in any_calsteps]).any():
        try:
            caltables = load_obj(calib_dir+'caltables.pkl')
            logger.info('Loaded previous calibration tables from: {0}'.format(calib_dir+'caltables.pkl'))
        except:
            caltables = {}
            caltables['inbase'] = inputs['inbase']
            caltables['plots_dir'] = plots_dir
            caltables['calib_dir'] = calib_dir
            caltables['num_spw'] = msinfo['num_spw']
        all_calsteps = [
               'bpcal_d.K0',
                   'bpcal_p.G0',
                   'bpcal_ap.G1',
                   'bpcal.B0',
                   'delay.K1',
                   'allcal_p.G0',
                   'allcal_p_jitter.G0',
                   'allcal_ap.G1',
                   'bpcal_sp.B1',
                   'allcal_ap.G3',
                   'phscal_p_scan.G3',
                   'phscal_ap_scan.G3']    # This is just used to know which tables to search in weblog
        caltables['all_calsteps'] = all_calsteps
        logger.info('New caltables dictionary created. Saved to: {0}'.format(calib_dir+'caltables.pkl'))
        caltables['Lo_dropout_scans'] = eMCP['msinfo']['Lo_dropout_scans']
        caltables['refant'] = msinfo['refant']
        save_obj(caltables, calib_dir+'caltables.pkl')

    ### Restore flag status at to this point
    if inputs['restore_flags'] == 1:
        em.restoreflagstatus(msinfo)

    ### Load manual flagging file
    if inputs['flag_manual_avg'] == 1:
        eMCP = em.flagdata_manual(eMCP, run_name='flag_manual_avg')

    ### Initialize models ###
    if inputs['init_models'] > 0:  # Need to add parameter to GUI
        eMCP = em.run_initialize_models(eMCP)

    ### Initial BandPass calibration ###
    if inputs['bandpass'] > 0:
        eMCP, caltables = em.initial_bp_cal(eMCP, caltables)

    ### Flagdata using TFCROP and bandpass shape B0
    if inputs['flag_tfcropBP'] > 0:
        eMCP = em.flagdata_tfcropBP(eMCP, caltables)

    ### Delay calibration ###
    if inputs['delay'] > 0:
        if not eMCP['defaults']['delay']['use_fringefit']:
            eMCP, caltables = em.solve_delays(eMCP, caltables)
        else:
            logger.info('Full fringe fit selected.')
            eMCP, caltables = em.delay_fringefit(eMCP, caltables)

    ### Initial gain calibration ###
    if inputs['gain_p_ap'] > 0:
        eMCP, caltables = em.initial_gaincal(eMCP, caltables)

        logger.info('FLAGDATA with mode RFLAG')
        flagdata(vis=msfile,
                 mode='rflag',
                 field=eMCP['msinfo']['sources']['calsources'],
                 correlation='ABS_LL,RR',
                 combinescans = False,
                 datacolumn='corrected',
                 winsize=3,
                 timedevscale=4.0,
                 freqdevscale=4.0,
                 extendflags=False,
                 action='apply',
                 flagbackup = True)
    
        ######## REPEAT   #########
    
        logger.warning('REPEATING delay and gain_p_ap calibration!!')
        ### Delay calibration ###
        if inputs['delay'] > 0:
            if not eMCP['defaults']['delay']['use_fringefit']:
                eMCP, caltables = em.solve_delays(eMCP, caltables)
            else:
                logger.info('Full fringe fit selected.')
                eMCP, caltables = em.delay_fringefit(eMCP, caltables)
    
        ### Initial gain calibration ###
        if inputs['gain_p_ap'] > 0:
            eMCP, caltables = em.initial_gaincal(eMCP, caltables)

    ### Flux scale ###
    if inputs['fluxscale'] > 0:
        eMCP, caltables = em.eM_fluxscale(eMCP, caltables)

    ### BandPass calibration with spectral index information ###
    if inputs['bandpass_sp'] > 0:
        eMCP, caltables = em.bandpass_sp(eMCP, caltables)

    ### Amplitude calibration including spectral information ###
    if inputs['gain_amp_sp'] > 0:
        eMCP, caltables = em.gain_amp_sp(eMCP, caltables)

    ### Apply calibration  ###
    if inputs['applycal_all'] > 0:
        eMCP = em.applycal_all(eMCP, caltables)

    ### RFLAG automatic flagging ###
    if inputs['flag_rflag'] > 0:
        eMCP = em.flagdata_rflag(eMCP)

    ### Produce some visibility plots ###
    if inputs['plot_corrected'] > 0:
        eMCP = emplt.make_4plots(eMCP, datacolumn='corrected')

    ### First images ###
    if inputs['first_images'] > 0:
        logger.info('Running statwt')
        statwt(vis=msfile)
        eMCP = em.run_first_images(eMCP)

    ### Plot flagstatistics ###
    run_flag_statistics = 1
    if run_flag_statistics > 0:
        try:
            emplt.flag_statistics(msinfo)
        except:
            pass

    # Keep important files
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    os.system('cp casa_eMCP.log {}casa_eMCP.log.txt'.format(info_dir))


    ### Run monitoring for bright sources:
    try:
        if inputs['monitoring'] > 0:
            caltables = em.monitoring(msfile=msfile, msinfo=msinfo,
                                             caltables=caltables,
                                             previous_cal=[''])
    except:
        pass

    try:
        os.system('mv ipython-*.log casa-*.log *.last ./logs')
        logger.info('Moved ipython-*.log casa-*.log *.last to ./logs')
    except:
        pass
    logger.info('Pipeline finished')
    logger.info('#################')

    return



# The  __name__ == "__main__" approach does not work for CASA.
try:
    if run_in_casa == True:
        # Running the pipeline from inside CASA
        print('Pipeline initialized. To run the pipeline within CASA use:')
        print('run_pipeline(inputs_path=<input file>)')
except:
    inputs = em.check_in(pipeline_path)
    run_pipeline(inputs=inputs)

