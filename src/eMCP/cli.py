#!/usr/bin/env python
"""
Command-line interface for eMERLIN CASA Pipeline
"""
import os
import sys
import argparse
import traceback
import shutil
import time
import importlib.resources
from datetime import datetime
import numpy as np
import yaml

from eMCP.functions import eMCP_functions as em
from eMCP.utils import eMCP_utils as emutils
from eMCP.plots import eMCP_plots as emplt
from eMCP.utils.weblog_config import get_weblog_function
from ._version import __version__

# Initialize logger
logger = emutils.get_logger()

def run_pipeline(inputs_file='./inputs.ini', run_steps=[], skip_steps=[]):
    """Main function to run the pipeline with specified steps"""
    # Create directory structure
    pipeline_path = os.path.dirname(os.path.realpath(__file__))
    logger.info(f'Executing pipeline in: {pipeline_path}')
    calib_dir, info_dir = emutils.create_dir_structure()

    # Initialize eMCP dictionary, or continue with previous pipeline configuration if possible:
    eMCP = emutils.start_eMCP_dict(info_dir)

    # Get git info about pipeline version
    installed_version = emutils.get_pipeline_version()
    pipeline_version = __version__

    logger.info('Starting pipeline')
    logger.info('Running pipeline from:')
    logger.info(f'{pipeline_path}')
    logger.info(f'Pipeline version: {pipeline_version}')
    logger.info(f'Pipeline installed version: {installed_version}')
    logger.info('This log uses UTC times')
    eMCP['pipeline_path'] = pipeline_path
    emutils.check_pipeline_conflict(eMCP, pipeline_version)
    eMCP['pipeline_version'] = pipeline_version
    emutils.save_obj(eMCP, info_dir + 'eMCP_info.yaml')

    # Load default parameters from YAML
    yaml_file = './default_params.yaml'
    package_yaml = os.path.join(pipeline_path, 'default_params.yaml')
    
    # Check for local file first, then fall back to package default
    if os.path.isfile(yaml_file):
        defaults_file = yaml_file
    else:
        defaults_file = package_yaml
    
    # Load the YAML file
    print('Loading default parameters from {0}:'.format(defaults_file))
    with open(defaults_file, 'r') as f:
        eMCP['defaults'] = yaml.safe_load(f)

    # Inputs
    if os.path.exists(inputs_file):
        inputs = emutils.read_inputs(inputs_file)
        eMCP['inputs'] = inputs
    else:
        print('No inputs file found: {}'.format(inputs_file))
        emutils.exit_pipeline(eMCP='')

    # Steps to run:
    eMCP['input_steps'] = emutils.find_run_steps(eMCP, run_steps, skip_steps)

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
        emutils.save_obj(caltables, os.path.join(calib_dir, 'caltables.yaml'))

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

    try:
        os.system('mv casa-*.log *.last ./logs')
        print('Moved casa-*.log *.last to ./logs')
    except:
        pass

def init_project(force=False):
    """Copy config files to current directory"""
    # Initialize success flag
    success = True
    
    # Set up destination paths
    dest_inputs = './inputs.ini'
    dest_params_yaml = './default_params.yaml'
    
    # Handle inputs.ini
    if os.path.exists(dest_inputs) and not force:
        logger.warning(f"{dest_inputs} already exists. Use --force to overwrite.")
        success = False
    else:
        try:
            # Try to get the resource from the package
            inputs_content = importlib.resources.read_text('eMCP', 'inputs.ini')
            with open(dest_inputs, 'w') as f:
                f.write(inputs_content)
            logger.info(f"Created {dest_inputs}")
        except Exception as e:
            logger.error(f"Error creating {dest_inputs}: {e}")
            success = False
    
    # Handle default_params.yaml
    if os.path.exists(dest_params_yaml) and not force:
        logger.warning(f"{dest_params_yaml} already exists. Use --force to overwrite.")
        success = False
    else:
        try:
            # Try to get the resource from the package
            params_content = importlib.resources.read_text('eMCP', 'default_params.yaml')
            with open(dest_params_yaml, 'w') as f:
                f.write(params_content)
            logger.info(f"Created {dest_params_yaml}")
        except Exception as e:
            logger.error(f"Error creating {dest_params_yaml}: {e}")
            success = False
    
    return success

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='eMERLIN CASA Pipeline',
        epilog='''
Examples:
  emcp --init               # Create default config files in current directory
  emcp --init --force       # Overwrite existing config files
  emcp -r "flag_apriori flag_manual average"  # Run specific steps (space or comma-separated)
  emcp -l                   # List available steps
''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Main options
    parser.add_argument('--init', action='store_true',
                       help='Initialize a new project directory with default config files')
    parser.add_argument('--force', action='store_true', 
                       help='Overwrite existing files when using --init')
    parser.add_argument('-i', '--inputs', dest='inputs_file', default='./inputs.ini',
                        help='Inputs file [./inputs.ini]')
    parser.add_argument('-r', '--run-steps', dest='run_steps', default='',
                        help='List of steps to run (space or comma-separated: "flag_apriori flag_manual average")')
    parser.add_argument('-s', '--skip-steps', dest='skip_steps', default='',
                        help='List of steps to skip (space or comma-separated: "plot_data save_flags")')
    parser.add_argument('-l', '--list-steps', dest='list_steps', action='store_true',
                        help='List all available steps')
    parser.add_argument('-v', '--version', action='store_true',
                        help='Show version information')
    
    return parser.parse_args()

def main():
    """Entry point for the application"""
    args = get_args()
    
    # Handle the init flag
    if args.init:
        init_project(force=args.force)
        return
    
    # Handle the version flag
    if args.version:
        logger.info(f"eMERLIN CASA Pipeline (eMCP) version {__version__}")
        return
        
    # Handle the list steps flag
    if args.list_steps:
        # Organize steps into categories for better readability
        logger.info('Available steps:')
        logger.info('pre_processing:')
        logger.info('    run_importfits, flag_aoflagger, flag_apriori, flag_manual, average, plot_data, save_flags')
        logger.info('calibration:')
        logger.info('    restore_flags, flag_manual_avg, init_models, bandpass, initial_gaincal, fluxscale, bandpass_final, gaincal_final', 'applycal_all', 'flag_target', 'plot_corrected', 'first_images')
        return

    run_steps = []
    if args.run_steps:
        # Support both comma and space as separators
        steps_str = args.run_steps.replace(',', ' ')
        run_steps = [step.strip() for step in steps_str.split() if step.strip()]
    
    skip_steps = []
    if args.skip_steps:
        # Support both comma and space as separators
        steps_str = args.skip_steps.replace(',', ' ')
        skip_steps = [step.strip() for step in steps_str.split() if step.strip()]
    
    run_pipeline(args.inputs_file, run_steps, skip_steps)

if __name__ == '__main__':
    main()
