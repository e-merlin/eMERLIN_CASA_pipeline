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

pipeline_version = 'v0.7.11'

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
#import functions.eMERLIN_CASA_GUI as emGUI

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
    fits_path = em.backslash_check(inputs['fits_path'])
    calib_dir = './calib/'
    plots_dir = './plots/'
    logs_dir  = './logs/'
    images_dir = './images/'

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
    em.makedir(logs_dir)
    em.makedir(images_dir)
    em.makedir(plots_dir+'caltables')

#    if inputs['quit'] == 1: #Check from GUI if quit is needed
#        logger.debug('Pipeline exit')
#        sys.exit()

    msfile = inputs['inbase']+'.ms'
    logger.info('Using MS file: {0}'.format(msfile))

    #################################
    ### LOAD AND PREPROCESS DATA  ###
    #################################

    ## Pipeline processes, inputs are read from the inputs dictionary
    if inputs['run_importfits'] > 0:
        em.run_importfitsIDI(fits_path, msfile, doaverage=inputs['run_importfits'])
        em.check_mixed_mode(msfile,mode='split')

    ### Write summary weblog ###
    if inputs['summary_weblog'] > 0:
        msinfo = em.get_msinfo(msfile, inputs)
        save_obj(msinfo, msfile+'.msinfo')
        logger.info('Saving information of MS {0} in: {1}'.format(msfile, msfile+'.pkl'))
        logger.info('Starting summary weblog')
        if not os.path.isdir(msfile):
            logger.info('Error finding original data: {0}'.format(msfile))
            logger.info('summary_weblog cannot be run. Exiting pipeline.')
            sys.exit()
        emplt.make_elevation(msfile, msinfo)
        emplt.make_uvcov(msfile, msinfo)
        emwlog.start_weblog(msinfo)

    if inputs['hanning'] > 0:
        run_hanning = inputs['hanning']
        em.hanning(inputvis=msfile, run_hanning=run_hanning, deloriginal=True)

    ### Convert MS to MMS ###
    if inputs['ms2mms'] > 0:
        em.ms2mms(vis=msfile,mode='parallel')

    ### check for parallelisation
    if os.path.isdir('./'+inputs['inbase']+'.mms') == True:
        msfile = inputs['inbase']+'.mms'
        msinfo = em.get_msinfo(msfile, inputs)
        save_obj(msinfo, msfile+'.msinfo')
        logger.info('Saving information of MS {0} in: {1}'.format(msfile, msfile+'.pkl'))
        logger.info('Using MS file: {0}'.format(msfile))

    ### Run AOflagger
    if inputs['flag_0_aoflagger'] > 0:
        if inputs['flag_0_aoflagger'] == 1:
            separate_bands = True
        elif inputs['flag_0_aoflagger'] == 2:
            separate_bands = False
        else:
            logger.warning('flag_0_aoflagger can only be 1 or 2')
        flags = em.run_aoflagger_fields(msfile=msfile, separate_bands=separate_bands,
                                        flags=flags, fields='all', pipeline_path = pipeline_path)

    ### A-priori flagdata: Lo&Mk2, edge channels, standard quack
    if inputs['flag_1_apriori'] > 0:
        sources = em.user_sources(inputs)
        Lo_dropout_scans =inputs['Lo_dropout_scans']
        try:
            msinfo
        except:
            msinfo = em.get_msinfo(msfile, inputs)
            save_obj(msinfo, msfile+'.msinfo')
        flags = em.flagdata1_apriori(msfile=msfile, msinfo=msinfo,
                                     Lo_dropout_scans=Lo_dropout_scans, flags=flags, do_quack=True)

    ### Load manual flagging file
    if inputs['flag_2a_manual'] > 0:
        flags = em.flagdata2_manual(msfile=msfile, inpfile=inputs['manual_flags_a'], flags=flags)

    ### Multi phase center ###
    if inputs['shift_field_pos'] > 0:
        em.shift_all_positions(msfile)

    ### Average data ###
    if inputs['average_1'] > 0:
        sources = em.user_sources(inputs)
        timebin = '{}s'.format(inputs['average_1'])
        em.run_split(msfile, sources=sources, width=4, timebin=timebin)

    # Check if averaged data already generated
    if os.path.isdir('./'+inputs['inbase']+'_avg.mms') == True:
        msfile = './'+inputs['inbase']+'_avg.mms'
        avg_file = True
    elif os.path.isdir('./'+inputs['inbase']+'_avg.ms') == True:
        msfile = './'+inputs['inbase']+'_avg.ms'
        avg_file = True
    else:
        avg_file = False

    ### Load or create dictionary with ms information.
    if avg_file == True:
        logger.info('Using MS file: {0}'.format(msfile))
        msinfo = em.get_msinfo(msfile, inputs)
        save_obj(msinfo, msfile+'.msinfo')
    else:
        try:
            msinfo   # Don't produce it again if summary_weblog was just executed.
        except:
            msinfo = em.get_msinfo(msfile, inputs)
            save_obj(msinfo, msfile+'.msinfo')

    ### Produce some plots ###
    if inputs['plot_data'] == 1:
        emplt.make_4plots(msfile, msinfo, datacolumn='data')

    ### Save flag status up to this point
    if inputs['save_flags'] == 1:
        em.saveflagstatus(msinfo)


    ###################
    ### CALIBRATION ###
    ###################

    # All the calibration steps will be saved in the dictionary caltables.pkl
    # located in the calib directory. If it does not exist a new one is created.
    any_calsteps = ['bandpass_0', 'delay', 'flag_3_tfcropBP','gain_0_p_ap','fluxscale','bandpass_1_sp','gain_1_amp_sp','applycal_all', 'weblog']
    if np.array([inputs[cal]>0 for cal in any_calsteps]).any():
        try:
            caltables = load_obj(calib_dir+'caltables')
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
                    'phscal_p_scan.G2',
                    'bpcal_sp.B1',
                    'allcal_ap.G3',
                    'allcal_ap_scan.G3']    # This is just used to know which tables to search in weblog
            caltables['all_calsteps'] = all_calsteps
            logger.info('New caltables dictionary created. Saved to: {0}'.format(calib_dir+'caltables.pkl'))
        caltables['Lo_dropout_scans'] = inputs['Lo_dropout_scans']
        caltables['refant'] = msinfo['refant']
        save_obj(caltables, calib_dir+'caltables')

    ### Restore flag status at to this point
    if inputs['restore_flags'] == 1:
        em.restoreflagstatus(msinfo)

    ### Load manual flagging file
    if inputs['flag_2b_manual'] == 1:
        flags = em.flagdata2_manual(msfile=msfile, inpfile=inputs['manual_flags_b'], flags=flags)

    ### Initialize models ###
    if inputs['init_models'] > 0:  # Need to add parameter to GUI
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
    if inputs['flag_3_tfcropBP'] > 0:
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
        save_obj(msinfo, msfile+'.msinfo')


    ### RFLAG automatic flagging ###
    if inputs['flag_4_rflag'] > 0:
        try:
            msinfo['applycal_all'] == True
            flags = em.flagdata4_rflag(msfile=msfile, msinfo=msinfo, flags=flags)
        except:
            logger.warning('flag_4_rflag selected but applycal_all has not been run. RFLAG only works on calibrated data!')
            logger.warning('flag_4_rflag will not be executed.')

    ### Produce some visibility plots ###
    if inputs['plot_corrected'] > 0:
        emplt.make_4plots(msfile, msinfo, datacolumn='corrected')
        emplt.make_uvplt(msinfo)

    ### First images ###
    if inputs['first_images'] > 0:
        em.run_first_images(msinfo)

    ### Plot flagstatistics ###
    if inputs['flag_statistics'] > 0:
        emplt.flag_statistics(msinfo)

    ### Write weblog ###
    if inputs['weblog'] > 0:
        elevplot = './plots/plots_observation/{0}_elevation.png'.format(msinfo['msfilename'])
        if os.path.isfile(elevplot):
            logger.info('Elevation plot found.')
            logger.info('To regenerate elev and uvcov plots remove {}.'.format(elevplot))
        else:
            emplt.make_elevation(msfile, msinfo)
            emplt.make_uvcov(msfile, msinfo)
        emwlog.start_weblog(msinfo)

    ### Run monitoring for bright sources:
    try:
        if inputs['monitoring'] > 0:
            flags, caltables = em.monitoring(msfile=msfile, msinfo=msinfo,
                                             flags=flags, caltables=caltables,
                                             previous_cal=[''])
    except:
        pass

    try:
        os.system('mv casa-*.log *.last ./logs')
    except:
        pass
    logger.info('Pipeline finished')
    logger.info('#################')

    return inputs, msinfo



# The  __name__ == "__main__" approach does not work for CASA.
try:
    if run_in_casa == True:
        # Running the pipeline from inside CASA
        print('Pipeline initialized. To run the pipeline within CASA use:')
        print('inputs, caltables, msinfo = run_pipeline(inputs_path=<input file>)')
except:
    inputs = em.check_in(pipeline_path)
    run_pipeline(inputs=inputs)

