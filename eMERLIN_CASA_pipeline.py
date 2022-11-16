# Dependencies
import os
import sys
import json
import argparse

from eMCP.functions import eMCP_functions as em
from eMCP.utils import eMCP_utils as emutils
from eMCP.plots import eMCP_plots as emplt

current_version = 'v2.0.20'

pipeline_filename = sys.argv[0]

casalog_name = 'casa_eMCP.log'


def run_pipeline(inputs_file='./inputs.ini', run_steps=[], skip_steps=[]):
    #Create directory structure
    pipeline_path = os.path.dirname(os.path.realpath(__file__))
    logger.info(f'Executing pipeline in: {pipeline_path}')
    calib_dir, info_dir = emutils.create_dir_structure(pipeline_path)

    # Initialize eMCP dictionary, or continue with previous pipeline configuration if possible:
    eMCP = emutils.start_eMCP_dict(info_dir)

    # Get git info about pipeline version
    installed_version = emutils.get_pipeline_version()
    #    try:
    #        branch, short_commit = emutils.get_pipeline_version(pipeline_path)
    #    except:
    #        branch, short_commit = 'unknown', 'unknown'
    pipeline_version = current_version

    logger.info('Starting pipeline')
    logger.info('Running pipeline from:')
    logger.info('{}'.format(pipeline_path))
    #'#    logger.info('CASA version: {}'.format(casalith.version_string()))
    logger.info('Pipeline version: {}'.format(pipeline_version))
    logger.info('Pipeline installed version: {}'.format(installed_version))
    logger.info('This log uses UTC times')
    eMCP['pipeline_path'] = pipeline_path
    #'#    eMCP['casa_version'] = casalith.version_string()
    emutils.check_pipeline_conflict(eMCP, pipeline_version)
    eMCP['pipeline_version'] = pipeline_version
    emutils.save_obj(eMCP, info_dir + 'eMCP_info.pkl')

    # Load default parameters
    if os.path.isfile('./default_params.json'):
        defaults_file = './default_params.json'
    else:
        defaults_file = os.path.join(pipeline_path, 'default_params.json')
    logger.info('Loading default parameters from {0}:'.format(defaults_file))
    eMCP['defaults'] = json.loads(open(defaults_file).read())  #,
    #                                  object_pairs_hook=deunicodify_hook)

    # Inputs
    if os.path.exists(inputs_file):
        inputs = emutils.read_inputs(inputs_file)
        eMCP['inputs'] = inputs
    else:
        logger.critical('No inputs file found: {}'.format(inputs_file))
        emutils.exit_pipeline(eMCP='')

    # Steps to run:
    eMCP['input_steps'] = emutils.find_run_steps(eMCP, run_steps, skip_steps)

    #    # Update casa-data if requested
    #    if eMCP['defaults']['global']['update_casa-data']:
    #        logger.info('Updating casa-data')
    #        os.system('update-data')

    ##################################
    ###  LOAD AND PREPROCESS DATA  ###
    ##################################

    ## Pipeline processes, inputs are read from the inputs dictionary
    if eMCP['input_steps']['run_importfits'] > 0:
        eMCP = em.import_eMERLIN_fitsIDI(eMCP)

    if os.path.isdir('./' + inputs['inbase'] + '.ms') == True:
        msfile = inputs['inbase'] + '.ms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### check for parallelisation
    if os.path.isdir('./' + inputs['inbase'] + '.mms') == True:
        msfile = inputs['inbase'] + '.mms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### Run AOflagger
    if eMCP['input_steps']['flag_aoflagger'] > 0:
        eMCP = em.run_aoflagger_fields(eMCP)

    ### A-priori flagdata: Lo&Mk2, edge channels, standard quack
    if eMCP['input_steps']['flag_apriori'] > 0:
        eMCP = em.flagdata1_apriori(eMCP)

    ### Load manual flagging file
    if eMCP['input_steps']['flag_manual'] > 0:
        eMCP = em.flagdata_manual(eMCP, run_name='flag_manual')

    ### Average data ###
    if eMCP['input_steps']['average'] > 0:
        eMCP = em.run_average(eMCP)

    # Check if averaged data already generated
    if os.path.isdir('./' + inputs['inbase'] + '_avg.mms') == True:
        msfile = './' + inputs['inbase'] + '_avg.mms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)
    elif os.path.isdir('./' + inputs['inbase'] + '_avg.ms') == True:
        msfile = './' + inputs['inbase'] + '_avg.ms'
        eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)
        em.plot_elev_uvcov(eMCP)

    ### Produce some plots ###
    if eMCP['input_steps']['plot_data'] == 1:
        eMCP = emplt.make_4plots(eMCP, datacolumn='data')

    ### Save flag status up to this point
    if eMCP['input_steps']['save_flags'] == 1:
        eMCP = em.saveflagstatus(eMCP)

    ###################
    ### CALIBRATION ###
    ###################

    ### Initialize caltable dictionary
    caltables = em.initialize_cal_dict(eMCP)

    ### Restore flag status at to this point
    if eMCP['input_steps']['restore_flags'] == 1:
        eMCP = em.restoreflagstatus(eMCP)

    ### Load manual flagging file
    if eMCP['input_steps']['flag_manual_avg'] == 1:
        eMCP = em.flagdata_manual(eMCP, run_name='flag_manual_avg')
        caltables['Lo_dropout_scans'] = eMCP['msinfo']['Lo_dropout_scans']
        emutils.save_obj(caltables, os.path.join(calib_dir, 'caltables.pkl'))

    ### Initialize models ###
    if eMCP['input_steps']['init_models'] > 0:  # Need to add parameter to GUI
        eMCP = em.run_initialize_models(eMCP)

    ### Initial BandPass calibration ###
    if eMCP['input_steps']['bandpass'] > 0:
        eMCP, caltables = em.initial_bp_cal(eMCP, caltables)

    ### Initial gaincal = delay, p, ap ###
    if eMCP['input_steps']['initial_gaincal'] > 0:
        eMCP, caltables = em.initial_gaincal(eMCP, caltables)

    ### Flux scale ###
    if eMCP['input_steps']['fluxscale'] > 0:
        eMCP, caltables = em.eM_fluxscale(eMCP, caltables)

    ### BandPass calibration with spectral index information ###
    if eMCP['input_steps']['bandpass_final'] > 0:
        eMCP, caltables = em.bandpass_final(eMCP, caltables)

    ### Amplitude calibration including spectral information ###
    if eMCP['input_steps']['gaincal_final'] > 0:
        eMCP, caltables = em.gaincal_final(eMCP, caltables)

    ### Apply calibration  ###
    if eMCP['input_steps']['applycal_all'] > 0:
        eMCP = em.applycal_all(eMCP, caltables)

    ### RFLAG automatic flagging ###
    if eMCP['input_steps']['flag_target'] > 0:
        em.run_flag_target(eMCP)

    ### Produce some visibility plots ###
    if eMCP['input_steps']['plot_corrected'] > 0:
        eMCP = emplt.make_4plots(eMCP, datacolumn='corrected')

    ### First images ###
    if eMCP['input_steps']['first_images'] > 0:
        eMCP = em.run_first_images(eMCP)


#
#    ### Split fields ###
#    if eMCP['input_steps']['split_fields'] > 0:
#        eMCP = em.run_split_fields(eMCP)
#
#
#    # Keep important files
#    emutils.save_obj(eMCP, info_dir + 'eMCP_info.pkl')
#    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
#    os.system('cp casa_eMCP.log {}casa_eMCP.log.txt'.format(info_dir))
#
#    emwlog.start_weblog(eMCP)

    try:
        os.system('mv casa-*.log *.last ./logs')
        logger.info('Moved casa-*.log *.last to ./logs')
    except:
        pass
    logger.info('Pipeline finished')
    logger.info('#################')

    return eMCP


def get_args():
    '''This function parses and returns arguments passed in'''
    # Assign description to the help doc
    description = 'e-MERLIN CASA pipeline. Visit: https://github.com/e-merlin/eMERLIN_CASA_pipeline'
    usage = 'casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r [steps]'
    epilog = 'Select "-r calibration" to run only the calibration part (recommended)'
    parser = argparse.ArgumentParser(description=description,
                                     usage=usage,
                                     epilog=epilog)
    parser.add_argument('-i',
                        '--inputs',
                        dest='inputs_file',
                        help='Inputs file to use. Default is inputs.ini',
                        default='./inputs.ini')
    parser.add_argument('-r', '--run-steps', dest='run_steps',
                        type=str, nargs='+',
                        help='Whitespace separated list of steps to run. '\
                        'Apart from individual steps, it also accepts "all", '\
                        '"pre_processing" and "calibration"',
                        default=[])
    parser.add_argument('-s',
                        '--skip-steps',
                        dest='skip_steps',
                        type=str,
                        nargs='+',
                        help='Whispace separated list of steps to skip',
                        default=[])
    parser.add_argument('-l',
                        '--list-steps',
                        dest='list_steps',
                        action='store_true',
                        help='Show list of available steps and exit')
    parser.add_argument('-c', help='Ignore, needed for casa')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    if args.list_steps:
        emutils.list_steps()
    else:
        # Setup logger
        logger = emutils.get_logger()
        run_pipeline(inputs_file=args.inputs_file,
                     run_steps=args.run_steps,
                     skip_steps=args.skip_steps)
