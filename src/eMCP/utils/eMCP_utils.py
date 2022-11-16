import os
import sys
import pickle
import configparser
import time
import shutil
import casacore.tables
import eMCP

import logging

logger = logging.getLogger('logger')


# Utility functions
def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass


def rmdir(pathdir, message='Deleted:'):
    if os.path.exists(pathdir):
        try:
            shutil.rmtree(pathdir)
            logger.info('{0} {1}'.format(message, pathdir))
        except:
            logger.debug('Could not delete: {0} {1}'.format(message, pathdir))
            pass


def rmfile(pathdir, message='Deleted:'):
    if os.path.exists(pathdir):
        try:
            os.remove(pathdir)
            logger.info('{0} {1}'.format(message, pathdir))
        except:
            logger.debug('Could not delete: {0} {1}'.format(message, pathdir))
            pass


def mvdir(pathdir, outpudir):
    if os.path.exists(pathdir):
        try:
            shutil.move(pathdir, outpudir)
            logger.info('Moved: {0} {1}'.format(pathdir, outpudir))
        except:
            logger.debug('Could not move: {0} {1}'.format(pathdir, outpudir))
            pass


# Save and load dictionaries
def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def get_logger(LOG_FORMAT='%(asctime)s | %(levelname)s | %(message)s',
               DATE_FORMAT='%Y-%m-%d %H:%M:%S',
               LOG_NAME='logger',
               LOG_LEVEL=logging.INFO,
               LOG_FILE_INFO='eMCP.log'):

    log = logging.getLogger(LOG_NAME)
    log_formatter = logging.Formatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT)
    logging.Formatter.converter = time.gmtime

    #'#    if not run_in_casa:
    # comment this to suppress console output
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)

    # eMCP.log with pipeline log
    file_handler_info = logging.FileHandler(LOG_FILE_INFO, mode='a')
    file_handler_info.setFormatter(log_formatter)
    file_handler_info.setLevel(logging.INFO)
    log.addHandler(file_handler_info)

    log.setLevel(LOG_LEVEL)
    return log


def create_dir_structure(pipeline_path):
    # Paths to use
    weblog_dir = './weblog/'
    info_dir = './weblog/info/'
    calib_dir = './weblog/calib/'
    plots_dir = './weblog/plots/'
    logs_dir = './logs/'
    images_dir = './weblog/images/'

    ## Create directory structure ##
    makedir(weblog_dir)
    makedir(info_dir)
    makedir(plots_dir)
    makedir(calib_dir)
    makedir(images_dir)
    makedir(logs_dir)
    makedir(plots_dir + 'caltables')
    os.system('cp -p {0}/utils/emerlin-2.gif {1}'.format(
        pipeline_path, weblog_dir))
    os.system('cp -p {0}/utils/eMCP.css {1}'.format(pipeline_path, weblog_dir))
    os.system('cp -p {0}/utils/eMCP_logo.png {1}'.format(
        pipeline_path, weblog_dir))
    return calib_dir, info_dir


def prt_dict(d, pre=''):
    subdict = []
    for key in d.keys():
        if type(d[key]) == dict:
            subdict.append(key)
        else:
            print('{0:20s}: {1}'.format(pre + key, d[key]))
    if subdict != []:
        for key_inner in subdict:
            print(pre + key_inner)
            prt_dict(d[key_inner], pre=pre + '   ')


def prt_dict_tofile(d, tofilename=None, addfile='', pre=' '):
    if tofilename != None:
        f = open(tofilename, 'w')
    else:
        f = addfile
    subdict = []
    for key in d.keys():
        if type(d[key]) in [dict]:
            subdict.append(key)
        else:
            f.write('{0:20s}: {1}\n'.format(pre + key, d[key]))
    if subdict != []:
        for key_inner in subdict:
            f.write('{}\n'.format(pre + key_inner))
            prt_dict_tofile(d[key_inner], addfile=f, pre=pre + pre)


# Pipeline management


def list_steps():
    all_steps, pre_processing_steps, calibration_steps = list_of_steps()
    print('\npre_processing')
    for s in pre_processing_steps:
        print('    {}'.format(s))
    print('\ncalibration')
    for s in calibration_steps:
        print('    {}'.format(s))
    sys.exit()


def get_pipeline_version():
    return eMCP.__version__


def start_eMCP_dict(info_dir):
    try:
        eMCP = load_obj(info_dir + 'eMCP_info.pkl')
    except:
        eMCP = {}
        eMCP['steps'] = eMCP_info_start_steps()
        eMCP['img_stats'] = {}
    return eMCP


def list_of_steps():
    pre_processing_steps = [
        'run_importfits', 'flag_aoflagger', 'flag_apriori', 'flag_manual',
        'average', 'plot_data', 'save_flags'
    ]
    calibration_steps = [
        'restore_flags', 'flag_manual_avg', 'init_models', 'bandpass',
        'initial_gaincal', 'fluxscale', 'bandpass_final', 'gaincal_final',
        'applycal_all', 'flag_target', 'plot_corrected', 'first_images',
        'split_fields'
    ]
    all_steps = pre_processing_steps + calibration_steps
    return all_steps, pre_processing_steps, calibration_steps


def eMCP_info_start_steps():
    default_value = [0, 0, '']
    all_steps = list_of_steps()[0]

    steps = {}
    steps['start_pipeline'] = default_value
    for s in all_steps:
        steps[s] = default_value
    return steps


def check_pipeline_conflict(eMCP, pipeline_version):
    try:
        eMCP['pipeline_version']
        if eMCP['pipeline_version'] != pipeline_version:
            logger.warning(
                'The log shows that different versions of the pipeline'
                ' has been executed. Please verify versions')
            logger.warning('Previous version: {0}. Current version {1}'.format(
                eMCP['pipeline_version'], pipeline_version))
    except:
        pass


def read_inputs(inputs_file):
    config = configparser.ConfigParser()
    config.read(inputs_file)
    return config._sections['inputs']


def exit_pipeline(eMCP=''):
    if eMCP != '':
        logger.info('Something went wrong. Producing weblog before quiting')
    logger.info('Now quiting')
    sys.exit()


def find_run_steps(eMCP, run_steps, skip_steps=[]):
    if run_steps == '': run_steps = []
    if skip_steps == '': skip_steps = []
    logger.info('Step selection')
    logger.info('run_steps : {}'.format(run_steps))
    logger.info('skip_steps: {}'.format(skip_steps))

    all_steps, pre_processing_steps, calibration_steps = list_of_steps()

    # Populate list of steps selected
    step_list = []
    if 'pre_processing' in run_steps:
        step_list += pre_processing_steps
        run_steps.remove('pre_processing')
    if 'calibration' in run_steps:
        step_list += calibration_steps
        run_steps.remove('calibration')
    if 'all' in run_steps:
        step_list += all_steps
        run_steps.remove('all')
    step_list += run_steps

    # Check if all are valid steps:
    wrong_steps = [s for s in step_list if s not in all_steps]
    if wrong_steps != []:
        ws = ', '.join(wrong_steps)
        logger.critical('Not available step(s) to run: {0}'.format(ws))
        exit_pipeline(eMCP='')

    wrong_steps = [s for s in skip_steps if s not in all_steps]
    if wrong_steps != []:
        ws = ', '.join(wrong_steps)
        logger.critical('Not available step(s) to skip: {0}'.format(ws))
        exit_pipeline(eMCP='')

    # Remove skipped steps:
    for skip_step in skip_steps:
        if skip_step != '':
            step_list.remove(skip_step)

    # Define final step dictionary:
    logger.info('Sorted list of steps to execute:')
    input_steps = {}
    for s in all_steps:
        if s in step_list:
            logger.info('{0:16s}: {1}'.format(s,
                                              eMCP['defaults']['global'][s]))
            input_steps[s] = eMCP['defaults']['global'][s]
        elif s not in step_list:
            logger.info('{0:16s}: {1}'.format(s, 0))
            input_steps[s] = 0
        else:
            pass

    return input_steps


## CASACORE funcions


def read_keyword(infile, column, subtable=None):
    with casacore.tables.table(infile, ack=False) as maintable:
        if subtable != None:
            tb = casacore.tables.table(maintable.getkeyword(subtable),
                                       ack=False)
            maintable.close()
        else:
            tb = maintable
        column = tb.getcol(column)


#        tb.close()
    return column


def read_all_keywords(infile, subtable):
    with casacore.tables.table(infile, ack=False) as maintable:
        with casacore.tables.table(maintable.getkeyword(subtable),
                                   ack=False) as subtable:
            keywords = {
                column: subtable.getcol(column)
                for column in subtable.colnames()
            }
    return keywords


def read_caltable_data(caltable):
    with casacore.tables.table(caltable, ack=False) as maintable:
        colnames = maintable.colnames()
        dontread = ['WEIGHT', 'INTERVAL']
        data = {
            colname: maintable.getcol(colname)
            for colname in colnames if colname not in dontread
        }
    return data


def find_source_timerange(msfile, source):
    field_names = read_keyword(msfile, 'NAME', 'FIELD')
    source_id = [i for i, j in enumerate(field_names) if source == j][0]
    field_id = read_keyword(msfile, 'FIELD_ID')
    t = casacore.tables.table(msfile, ack=False)
    t1 = casacore.tables.taql('select from $t where FIELD_ID == $source_id')
    times = t1.getcol('TIME')
    t.close()
    t1.close()
    return times.min(), times.max()
