#!/usr/local/python
import os
import subprocess
import numpy as np
import socket
import pickle
import glob
import re
import itertools
import sys
import copy
import datetime
import dateutil
import shutil
from scipy.stats import mode
import logging
from astropy.coordinates import SkyCoord
import astropy.units as u
from casacore import tables as casacore_tables
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import configparser
from pathlib import Path

from astropy.io import fits

from ..plots import eMCP_plots as emplt
from ..weblog import eMCP_weblog as emwlog
from ..utils import eMCP_utils as emutils
from ..fluxscale import run_fluxscale
from ..flagstatistics import run_flagstats

from casatasks import mstransform, applycal, gaincal, flagmanager, flagdata, concat,\
    setjy, importfitsidi, listobs, vishead, visstat, fixvis, phaseshift, statwt, fluxscale, immath,\
    imstat, fringefit, smoothcal, bandpass, delmod, clearcal, initweights, ft, tclean

from casatools import table, msmetadata
from casaviewer import imview

tb = table()
msmd = msmetadata()

# Logging
logger = logging.getLogger('logger')

weblog_dir = './weblog/'
info_dir = './weblog/info/'
calib_dir = './weblog/calib/'
plots_dir = './weblog/plots/'
logs_dir = './logs/'
images_dir = './weblog/images/'

weblog_link = './'
info_link = './info/'
calib_link = './calib/'
plots_link = './plots/'
images_link = './images/'

line0 = '-' * 15


def backslash_check(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory


def read_inputs(inputs_file):
    config = configparser.ConfigParser()
    config.read(inputs_file)
    return config._sections['inputs']


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
    if wrong_steps:
        ws = ', '.join(wrong_steps)
        logger.critical('Not available step(s) to run: {0}'.format(ws))
        exit_pipeline(eMCP='')

    wrong_steps = [s for s in skip_steps if s not in all_steps]
    if wrong_steps:
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


def exit_pipeline(eMCP=''):
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    if eMCP != '':
        logger.info('Something went wrong. Producing weblog before quiting')
        emwlog.start_weblog(eMCP)
    logger.info('Now quiting')
    sys.exit()


def find_casa_problems():
    logger.debug('Checking casa_log for problems')
    # TODO: Here we need to make a list with all the casa log files
    # TODO: Then check the file that has the most recent date and
    # TODO: check for SEVERE words.
    pass


# Functions to save and load dictionaries
def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def add_step_time(step, eMCP, msg, t0, doweblog=True):
    t1 = datetime.datetime.utcnow()
    timestamp = t1.strftime('%Y-%m-%d %H:%M:%S')
    delta_t_min = (t1 - t0).total_seconds() / 60.
    eMCP['steps'][step] = [timestamp, delta_t_min, msg]
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    if doweblog:
        emwlog.start_weblog(eMCP)
    return eMCP


def check_pipeline_conflict(eMCP, pipeline_version):
    try:
        if eMCP['pipeline_version'] != pipeline_version:
            logger.warning(
                'The log shows that different versions of the pipeline'
                ' has been executed. Please verify versions')
            logger.warning('Previous version: {0}. Current version {1}'.format(
                eMCP['pipeline_version'], pipeline_version))
    except:
        pass


def update_mixed_mode(eMCP):
    # Default can be 'auto', 'none', 'force'
    default_mixed_mode = eMCP['defaults']['global']['is_mixed_mode']
    logger.debug('default_mixed_mode = {0}'.format(default_mixed_mode))

    msfile_sp = get_msfile_sp(eMCP)
    is_mixed_mode = False
    if default_mixed_mode == 'auto':
        if os.path.isdir(msfile_sp):
            logger.info('Found mixed mode file {0}'.format(msfile_sp))
            is_mixed_mode = True
        else:
            logger.info('No mixed mode file found')
            is_mixed_mode = False
    elif default_mixed_mode == 'none':
        is_mixed_mode = False
        logger.info('Forced no mixed mode')
    elif default_mixed_mode == 'force':
        logger.info('Forced mixed mode')
        if os.path.isdir(msfile_sp):
            is_mixed_mode = True
            logger.info('Found narrow band data {0}'.format(msfile_sp))
        else:
            logger.critical('Cannot narrow band data {}'.format(msfile_sp))
            is_mixed_mode = False
            exit_pipeline(eMCP)
    eMCP['is_mixed_mode'] = is_mixed_mode
    return eMCP


def check_band(eMCP, msfile):
    # Output example: 'C'
    # Take first frequency in the MS
    freq = emutils.read_keyword(
        msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW').mean() / 1e9
    band = ''
    if 1.2 < freq < 1.7:
        band = 'L'
    elif 4 < freq < 8:
        band = 'C'
    elif 17 < freq < 26:
        band = 'K'
    else:
        logger.critical('Cannot determine band from frequency {}'.format(freq))
        exit_pipeline(eMCP)
    logger.debug(f'Mean frequency in {msfile}; {freq}. Band {band}')
    return band


def get_baselines(msfile):
    antennas = emutils.read_keyword(msfile, 'NAME', subtable='ANTENNA')
    baselines = [
        f'{pair[0]}-{pair[1]}' for pair in itertools.combinations(antennas, 2)
    ]
    logger.debug(f'Baselines in {msfile}; {baselines}')
    return np.array(baselines)


def get_dates(d):
    t_mjd = d['axis_info']['time_axis']['MJDseconds'] / 60. / 60. / 24.
    t = np.array([mjdtodate(i) for i in t_mjd])
    return t_mjd, t


def mjdtodate(mjd):
    origin = datetime.datetime(1858, 11, 17)
    date = origin + datetime.timedelta(mjd)
    return date


def join_lists(x=[]):
    x1 = ','.join(x)
    x2 = set(x1.split(','))
    return ','.join(x2)


def user_sources(inputs):
    sources = {}
    sources['targets'] = inputs['targets']
    sources['phscals'] = inputs['phscals']
    sources['fluxcal'] = inputs['fluxcal']
    sources['bpcal'] = inputs['bpcal']
    sources['ptcal'] = inputs['ptcal']
    sources['calsources'] = join_lists([
        sources['phscals'], sources['fluxcal'], sources['bpcal'],
        sources['ptcal']
    ])
    sources['calsources'] = sources['calsources'].strip(',')
    sources['maincal'] = join_lists(
        [sources['fluxcal'], sources['bpcal'], sources['ptcal']])
    sources['allsources'] = join_lists(
        [sources['calsources'], sources['targets']])
    sources['allsources'] = sources['allsources'].strip(',')
    sources['no_fluxcal'] = sources['allsources'].replace(
        sources['fluxcal'], '').replace(',,', ',').strip(',')
    sources['cals_no_fluxcal'] = sources['calsources'].replace(
        sources['fluxcal'], '').replace(',,', ',').strip(',')
    sources['targets_phscals'] = join_lists(
        [sources['targets'], sources['phscals']])
    logger.debug('Targets:   {0}'.format(sources['targets']))
    logger.debug('Phasecals: {0}'.format(sources['phscals']))
    logger.debug('Fluxcal:   {0}'.format(sources['fluxcal']))
    logger.debug('Bandpass:  {0}'.format(sources['bpcal']))
    logger.debug('Pointcal:  {0}'.format(sources['ptcal']))
    return sources


def get_antennas(msfile):
    # Output example: ['Mk2', 'Pi', 'Da', 'Kn', 'De', 'Cm']
    # Antenna list
    antennas = emutils.read_keyword(msfile, 'NAME', subtable='ANTENNA')
    nice_order = ['Lo', 'Mk2', 'Pi', 'Da', 'Kn', 'De', 'Cm']
    antennas_sorted = [a for a in nice_order if a in antennas]
    external = [a for a in antennas if a not in antennas_sorted]
    antennas_all = antennas_sorted + external
    logger.debug('Antennas in MS {0}: {1}'.format(msfile, antennas_all))
    return antennas_all


def get_obstime(msfile):
    # returns datetime object of first and last times
    # Output example:  time objects: 2017-12-20 00:35:00, 2017-12-20 03:42:58
    times = emutils.read_keyword(msfile, 'TIME')
    t_ini = mjdtodate(times.min() / 60. / 60. / 24.)
    t_end = mjdtodate(times.max() / 60. / 60. / 24.)
    logger.debug(f'Ini and end times: {t_ini}, {t_end}')
    return t_ini, t_end


def get_obsfreq(msfile):
    # Returns freq of first channel, end chan, channel resolution
    # and number of channels (first spw) in GHz
    # Output example: (4.816125, 5.327875, 0.00025, 512)
    freq_ini = emutils.read_keyword(
        msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW').min() / 1e9
    freq_end = emutils.read_keyword(
        msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW').max() / 1e9
    chan_res = emutils.read_keyword(
        msfile, 'RESOLUTION', subtable='SPECTRAL_WINDOW').mean() / 1e9
    nchan = emutils.read_keyword(msfile,
                                 'CHAN_FREQ',
                                 subtable='SPECTRAL_WINDOW').shape[1]
    logger.debug(
        f'freq_ini, freq_end, chan_res, nchan: {freq_ini}, {freq_end}, {chan_res}, {nchan}'
    )
    return freq_ini, freq_end, chan_res, nchan


def find_mssources(msfile):
    # Output example: '1107-1226,1109-1235,1118-1232,1331+305,1407+284'
    fieldnames = emutils.read_keyword(msfile, 'NAME', subtable='FIELD')
    mssources = ','.join(np.sort(fieldnames))
    logger.debug('Sources in MS {0}: {1}'.format(msfile, mssources))
    return mssources


def find_source_intent(msinfo, cats=None):
    if cats is None:
        cats = ['targets', 'phscals', 'bpcal', 'fluxcal', 'ptcal']
    fields_ms = msinfo['sources']['mssources'].split(',')
    return {
        source: ','.join([
            cat for cat in cats if source in msinfo['sources'][cat].split(',')
        ])
        for source in fields_ms
    }


def find_source_timerange(msfile):
    fieldnames = emutils.read_keyword(msfile, 'NAME', subtable='FIELD')
    source_timerange_mjd = {}
    fact = 60. * 60. * 24.
    for source in fieldnames:
        mjd_min, mjd_max = emutils.find_source_timerange(msfile, source)
        source_timerange_mjd[source] = [mjd_min / fact, mjd_max / fact]
    return source_timerange_mjd


def get_project(msfile):
    # Output example: 'CY0000'
    project = emutils.read_keyword(msfile, 'PROJECT', subtable='OBSERVATION')
    logger.debug(f'Read project from {msfile}: {project}')
    return project


def get_polarization(msfile):
    # Output example: 'L, R'
    pol_types = np.unique(
        emutils.read_keyword(msfile, 'POLARIZATION_TYPE', 'FEED')['array'])
    polarization = ', '.join(pol_types)
    logger.debug(f'Polarization types for {msfile}: {polarization}')
    return polarization


def get_directions(msfile):
    directions = {}
    field_names = emutils.read_keyword(msfile, 'NAME', 'FIELD')
    phase_dir = emutils.read_keyword(msfile, 'PHASE_DIR', 'FIELD')
    for i, field in enumerate(field_names):
        ra = phase_dir[i][0][0] * u.rad
        dec = phase_dir[i][0][1] * u.rad
        # TODO: Frame need to be read from Measurement Set
        directions[field] = SkyCoord(ra, dec, frame='icrs')
    return directions


def get_distances(msfile, directions=''):
    if directions == '':
        directions = get_directions(msfile)
    field_names = emutils.read_keyword(msfile, 'NAME', 'FIELD')
    separations = {}
    # Write all separations in a txt file
    with open(os.path.join(info_dir, 'source_separations.txt'),
              'w') as sep_file:
        for (f1, f2) in itertools.combinations(field_names, 2):
            separations[f1 + '-' + f2] = directions[f1].separation(
                directions[f2]).to(u.deg)
            sep_file.write('{0:10} {1:10} {2:7.2f}\n'.format(
                f1, f2, separations[f1 + '-' + f2]))
    return separations


def get_integration_time(msfile):
    time_diff = np.diff(emutils.read_keyword(msfile, 'TIME'))
    intervals = time_diff[time_diff > 0]
    int_mode = mode(intervals)
    int_time = int_mode[0][0]
    perc = int_mode[1][0] / len(intervals) * 100.
    logger.debug(f'Mode time interval: {int_time}s (for {perc:4.1f}% of gaps)')
    return int_time


def get_msfile_sp(eMCP):
    ext_ms = {False: '.ms', True: '.mms'}
    do_ms2mms = eMCP['defaults']['import_eM']['ms2mms']
    msfile_sp = os.path.join(
        '/', eMCP['inputs']['inbase'] + '_sp' + ext_ms[do_ms2mms])
    return msfile_sp


def find_wide_narrow(spw_sp, cent_chan_sp, msfile):
    msmd.open(msfile)
    spws = msmd.spwfordatadesc()
    main_spw = []
    logger.debug('')
    logger.debug('Narrow sp {0}, central freq. {1:5.3f}'.format(
        spw_sp, cent_chan_sp / 1e9))
    for spw in spws:
        chan = msmd.chanfreqs(spw)
        freq_ini = chan[0]
        freq_end = chan[-1]
        logger.debug('Wide sp {0}: {1:5.3f} {2:5.3f}'.format(
            spw, freq_ini / 1e9, freq_end / 1e9))
        if cent_chan_sp > freq_ini and cent_chan_sp <= freq_end:
            main_spw.append(spw)
    logger.debug('Corresponding wide spw {0}'.format(main_spw))
    msmd.done()
    # print results
    if len(main_spw) == 0:
        logger.warning('Could not find wideband spw containing narrow band ' \
              'centred at {0}. Choose spwmap manually!'.format(cent_chan_sp*1e9))
    elif len(main_spw) > 1:
        logger.warning('Narrow spw {0} can be linked ' \
              'to different wide spw {1}. Choose spwmap manually!'.format(spw_sp, main_spw))
    else:
        pass
    return main_spw


def get_cent_freq(msfile):
    msmd.open(msfile)
    spws = msmd.spwfordatadesc()
    cent_freq = np.array([np.mean(msmd.chanfreqs(spw)) for spw in spws])
    msmd.done()
    return cent_freq


def get_chan_width(msfile):
    msmd.open(msfile)
    spws = msmd.spwfordatadesc()
    chan_width = np.array([np.unique(msmd.chanwidths(spw))[0] for spw in spws])
    msmd.done()
    return chan_width


def get_spwmap_sp(msfile, msfile_sp):
    msmd.open(msfile_sp)
    spws_sp = msmd.spwfordatadesc()
    msmd.done()
    cent_freq_sp = get_cent_freq(msfile_sp)
    spwmap_sp = []
    for (spw_sp, cent_chan_sp) in zip(spws_sp, cent_freq_sp):
        main_spw = find_wide_narrow(spw_sp, cent_chan_sp, msfile)
        spwmap_sp.append(main_spw[0])
    return spwmap_sp


def find_spwmap_sp(eMCP, msfile, msfile_sp):
    if eMCP['defaults']['import_eM']['spwmap_sp'] == []:
        spwmap_sp = get_spwmap_sp(msfile, msfile_sp)
    else:
        spwmap_sp = eMCP['defaults']['import_eM']['spwmap_sp']
    return spwmap_sp


def info_mixed_mode(eMCP, msinfo):
    eMCP = update_mixed_mode(eMCP)
    # Deal with mixed mode info
    if eMCP['is_mixed_mode']:
        logger.debug(eMCP['is_mixed_mode'])
        logger.debug('msfile: {0}'.format(msinfo['msfile']))
        logger.debug('msfile narrow: {0}'.format(msinfo['msfile_sp']))
        msinfo['spwmap_sp'] = find_spwmap_sp(eMCP, msinfo['msfile'],
                                             msinfo['msfile_sp'])
        logger.info('spwmap_sp = {0}'.format(msinfo['spwmap_sp']))
        msinfo['spwmap_sp_freq'] = get_cent_freq(msinfo['msfile_sp'])
        msinfo['spwmap_sp_width'] = get_chan_width(msinfo['msfile_sp'])
    else:
        msinfo['spwmap_sp'] = []
        msinfo['spwmap_sp_freq'] = []
        msinfo['spwmap_sp_width'] = []
    return eMCP, msinfo


def get_msinfo(eMCP, msfile, doprint=False):
    inputs = eMCP['inputs']
    logger.info('Found MS file: {0}'.format(msfile))
    logger.info('Reading ms file information for MS: {0}'.format(msfile))
    msinfo = {}
    msinfo['msfile'] = msfile
    msinfo['msfile_sp'] = get_msfile_sp(eMCP)
    msinfo['msfilename'] = os.path.splitext(msfile)[0].split('/')[-1]
    msinfo['project'] = get_project(msfile)
    msinfo['run'] = inputs['inbase']
    msinfo['sources'] = user_sources(inputs)
    msinfo['sources']['mssources'] = find_mssources(msfile)
    msinfo['sources']['source_intent'] = find_source_intent(msinfo)
    msinfo['sources']['source_timerange_mjd'] = find_source_timerange(
        msinfo['msfile'])
    msinfo['antennas'] = get_antennas(msfile)
    msinfo['band'] = check_band(eMCP, msfile)
    msinfo['baselines'] = get_baselines(msfile)
    #    msinfo['num_spw'] = len(casatasks.vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    msinfo['num_spw'] = len(
        emutils.read_keyword(msfile,
                             'MEAS_FREQ_REF',
                             subtable='SPECTRAL_WINDOW'))
    t_ini, t_end = get_obstime(msfile)
    freq_ini, freq_end, chan_res, nchan = get_obsfreq(msfile)
    msinfo['t_ini'] = t_ini
    msinfo['t_end'] = t_end
    msinfo['freq_ini'] = freq_ini
    msinfo['freq_end'] = freq_end
    msinfo['int_time'] = get_integration_time(msfile)
    msinfo['chan_res'] = chan_res
    msinfo['nchan'] = nchan
    msinfo['innerchan'] = '{0:.0f}~{1:.0f}'.format(
        0.1 * (nchan - nchan / 512.), 0.9 * (nchan - nchan / 512.))
    msinfo['polarizations'] = get_polarization(msfile)
    try:
        msinfo['refant'] = eMCP['msinfo']['refant']
    except:
        msinfo['refant'] = ''
    msinfo['directions'] = get_directions(msfile)
    msinfo['separations'] = get_distances(msfile,
                                          directions=msinfo['directions'])
    msinfo['Lo_dropout_scans'] = msinfo.get('Lo_dropout_scans', 'none')
    # Info related to mixed mode:
    eMCP, msinfo = info_mixed_mode(eMCP, msinfo)
    # Show summary
    logger.info('> Sources ({0}): {1}'.format(
        len(msinfo['sources']['mssources'].split(',')),
        msinfo['sources']['mssources']))
    logger.info('> Number of spw: {0}'.format(msinfo['num_spw']))
    logger.info('> Channels per spw: {0}'.format(msinfo['nchan']))
    logger.info('> Integration time {0:3.1f}s'.format(msinfo['int_time']))
    if doprint:
        emutils.prt_dict(msinfo)
    eMCP['msinfo'] = msinfo
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    return eMCP, msinfo, msfile


def get_unique_field(caltable):
    # If more than one field was used to generate the caltable but
    # combine='field' was used, only one field is included in the table. This
    # function will find which was actually used.
    field_id = np.unique(emutils.read_keyword(caltable, 'FIELD_ID'))[0]
    unique_field = emutils.read_keyword(caltable, 'NAME',
                                        subtable='FIELD')[field_id]
    return unique_field


def backup_table(caltable):
    try:
        emutils.rmdir(caltable + '_backup_missing')
    except OSError:
        pass
    shutil.copytree(caltable, caltable + '_backup_missing')
    logger.info('Backup of table {0} to {1}'.format(
        caltable, caltable + '_backup_missing'))


def remove_missing_scans(caltable, scans2flag):
    # Backup original table
    backup_table(caltable)
    antenna_names = emutils.read_keyword(caltable, 'NAME', subtable='ANTENNA')

    if 'Lo' in antenna_names:
        anten_Lo = np.argwhere(np.array(antenna_names) == 'Lo')[0][0]
    else:
        anten_Lo = None

    scan_number = emutils.read_keyword(caltable, 'SCAN_NUMBER')
    antenna1_Lo = emutils.read_keyword(caltable, 'ANTENNA1') == anten_Lo
    missing_rows = np.array([
        str(r) in scans2flag.replace(' ', '').split(',') for r in scan_number
    ])
    index_missing_rows = np.where(antenna1_Lo * missing_rows)[0]
    logger.info('Removing Lo solutions for dropout scans from {0}: {1}'.format(
        caltable, scans2flag))
    with casacore_tables.table(caltable, ack=False,
                               readonly=False) as maintable:
        maintable.removerows(index_missing_rows)


def run_listobs(msfile):
    outfile = os.path.join(info_dir, msfile.rstrip('/'), '.listobs.txt')
    listobs(vis=msfile, listfile=outfile, overwrite=True)
    find_casa_problems()
    logger.info('Listobs file in: {0}'.format(outfile))


def decide_hanning(import_eM, msfile):
    run_hanning = import_eM['run_hanning']
    if run_hanning == 'auto':
        logger.debug('run_hanning = "auto" means check if data are L band.')
        band = check_band(import_eM, msfile)
        if band == 'L':
            logger.info('L band dataset. Hanning smoothing will be executed.')
            apply_hanning = True
        else:
            logger.info('Dataset is not L band. Hanning smoothing not needed.')
            apply_hanning = False
    elif run_hanning == 1:
        apply_hanning = True
    elif run_hanning == 0:
        apply_hanning = False
    else:
        logger.warning('Not valid run_hanning parameter (use "auto", 0 or 1)')
        exit_pipeline()
    return apply_hanning


def mixed_mode(msfile):
    bw_spw = emutils.read_keyword(msfile,
                                  'TOTAL_BANDWIDTH',
                                  subtable='SPECTRAL_WINDOW')
    is_mixed_mode = len(np.unique(bw_spw)) != 1
    logger.info('Natural spw and bandwidth of imported data:')
    for i, bw_spw_i in enumerate(bw_spw):
        logger.info('{0}: {1}'.format(i, bw_spw_i))
    if is_mixed_mode:
        cont_spw = ','.join(
            np.where(bw_spw == np.max(np.unique(bw_spw)))[0].astype('str'))
        line_spw = ','.join(
            np.where(bw_spw != np.max(np.unique(bw_spw)))[0].astype('str'))
        spw_separation = [cont_spw, line_spw]  # Example ['0,1,2,3','4,5']
    else:
        spw_separation = ['', '']
    return is_mixed_mode, spw_separation


def plot_elev_uvcov(eMCP):
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    elevplot = plots_dir + 'plots_observation/{0}_elevation.png'.format(
        eMCP['msinfo']['msfilename'])
    if os.path.isfile(elevplot):
        logger.info('Elevation plot found.')
        logger.info('To regenerate elev and uvcov plots remove:')
        logger.info('{}.'.format(elevplot))
    else:
        emplt.make_elevation(msfile, msinfo)
        emplt.make_uvcov(msfile, msinfo)
    emwlog.start_weblog(eMCP)


def import_eMERLIN_fitsIDI(eMCP):
    import_eM = eMCP['defaults']['import_eM']
    logger.info('Start run_importfits')
    fits_path = backslash_check(eMCP['inputs']['fits_path'])
    msfile_name = eMCP['inputs']['inbase']
    msg = 'first execution'
    t0 = datetime.datetime.utcnow()
    eMCP = add_step_time('start_pipeline', eMCP, msg, t0, doweblog=False)

    # importfitsIDI
    constobsid = import_eM['constobsid']
    scanreindexgap_s = import_eM['scanreindexgap_s']
    fitsfiles = []
    for infile in os.listdir(fits_path):
        if infile.endswith('fits') or infile.endswith('FITS'):
            fitsfiles = fitsfiles + [fits_path + infile]
            logger.info('FITS file found to be imported: {0}'.format(infile))
    if fitsfiles == []:
        logger.critical('No fits files found in {}'.format(fits_path))
        exit_pipeline(eMCP='')

    emutils.rmdir(eMCP['inputs']['inbase'] + '.ms')
    emutils.rmdir(eMCP['inputs']['inbase'] + '.ms.flagversions')
    emutils.rmdir(eMCP['inputs']['inbase'] + '.mms')
    emutils.rmdir(eMCP['inputs']['inbase'] + '.mms.flagversions')
    emutils.rmdir(eMCP['inputs']['inbase'] + '_sp.ms')
    emutils.rmdir(eMCP['inputs']['inbase'] + '_sp.mms')

    logger.info('Running importfitsIDI')
    t0 = datetime.datetime.utcnow()
    msfile0 = msfile_name + '_imported.ms'
    emutils.rmdir(msfile0)

    importfitsidi(fitsidifile=fitsfiles,
                  vis=msfile0,
                  constobsid=constobsid,
                  scanreindexgap_s=scanreindexgap_s)
    find_casa_problems()
    logger.info('Created file: {}'.format(msfile0))
    logger.info('Removing initial correlator flags')
    flagdata(vis=msfile0, mode='unflag', flagbackup=False)
    find_casa_problems()
    logger.info('Finished importfitsIDI')
    msg = 'constobsid={0}, scanreindexgap_s={1}'.format(
        constobsid, scanreindexgap_s)
    eMCP['msfile'] = eMCP['inputs']['inbase'] + '.ms'
    msfile = msfile0
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    emwlog.start_weblog(eMCP)
    #    eMCP = add_step_time('importfitsIDI', eMCP, msg, t0, doweblog=True)

    # mstransform
    datacolumn = 'data'
    if import_eM['antenna'] == '':
        antenna = '*&*'
    else:
        antenna = '*&*;' + import_eM['antenna']
        logger.info('Splitting antennas: {}'.format(antenna))
    field = import_eM['field']
    timeaverage = import_eM['timeaverage']
    timebin = import_eM['timebin']
    chanaverage = import_eM['chanaverage']
    chanbin = import_eM['chanbin']
    usewtspectrum = import_eM['usewtspectrum']
    do_hanning = decide_hanning(import_eM, msfile0)
    do_ms2mms = import_eM['ms2mms']
    is_mixed_mode, spw_separation = mixed_mode(msfile0)
    if import_eM['spw_separation'] != ['', '']:
        spw_separation = import_eM['spw_separation']
        logger.info('Using user-defined spw_separation:')
    logger.info('Continuum:   {}'.format(spw_separation[0]))
    if is_mixed_mode:
        logger.info('Narrow (sp): {}'.format(spw_separation[1]))
    ext_ms = {False: '.ms', True: '.mms'}
    msfile1 = msfile_name + '_transformed' + ext_ms[do_ms2mms]
    if timeaverage:
        logger.info('Data will be averaged to {}'.format(timebin))
    if chanaverage:
        logger.info('Data will be averaged to {} chan/spw'.format(chanbin))
    if do_hanning:
        logger.info('Running Hanning smoothing')
    else:
        logger.info('Not running Hanning smoothing')
    if do_ms2mms:
        logger.info('Data will be converted to MMS')
    logger.info('Start mstransform')
    t0 = datetime.datetime.utcnow()
    emutils.rmdir(msfile1)
    if import_eM['fix_repeated_sources'] == True:
        fix_repeated_sources(msfile0, msfile1, datacolumn, antenna,
                             str(spw_separation[0]), timeaverage, timebin,
                             chanaverage, chanbin, usewtspectrum, do_hanning,
                             do_ms2mms)
    else:
        mstransform(vis=msfile0,
                    outputvis=msfile1,
                    keepflags=True,
                    datacolumn=datacolumn,
                    antenna=antenna,
                    spw=str(spw_separation[0]),
                    timeaverage=timeaverage,
                    timebin=timebin,
                    chanaverage=chanaverage,
                    chanbin=chanbin,
                    usewtspectrum=usewtspectrum,
                    hanning=do_hanning,
                    createmms=do_ms2mms)
        find_casa_problems()
    logger.info('Transformed: {0} into {1}'.format(msfile0, msfile1))
    if is_mixed_mode:  # No hanning and no channel average
        logger.info('Mixed mode data detected')
        msfile1_sp = msfile_name + '_transformed_sp' + ext_ms[do_ms2mms]
        logger.info('Running mstransform on spectral line data')
        emutils.rmdir(msfile1_sp)
        if import_eM['fix_repeated_sources'] == True:
            fix_repeated_sources(msfile0, msfile1_sp, datacolumn, antenna,
                                 str(spw_separation[1]), timeaverage, timebin,
                                 chanaverage, chanbin, usewtspectrum,
                                 do_hanning, do_ms2mms)
        else:
            mstransform(vis=msfile0,
                        outputvis=msfile1_sp,
                        keepflags=True,
                        datacolumn=datacolumn,
                        antenna=antenna,
                        spw=str(spw_separation[1]),
                        timeaverage=timeaverage,
                        timebin=timebin,
                        usewtspectrum=usewtspectrum,
                        createmms=do_ms2mms)
            find_casa_problems()
        logger.info('Transformed: {0} into {1}'.format(msfile0, msfile1_sp))
    if os.path.isdir(msfile1) == True:
        emutils.rmdir(msfile0)
    else:
        logger.critical('Problem generating {}. Stopping ' \
                        'pipeline'.format(msfile1))
        exit_pipeline(eMCP)
    logger.info('Finished mstransform')
    msg += '. Hanning={0}, createmms={1}'.format(do_hanning, do_ms2mms)
    if timeaverage:
        msg += ', timebin={0}'.format(timebin)
    if chanaverage:
        msg += ', chanbin={0}'.format(chanbin)
    if import_eM['antenna'] != '':
        msg += ', antenna="{}"'.format(antenna)
    msfile = msfile1
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    #    eMCP = add_step_time('mstransform', eMCP, msg, t0, doweblog=True)

    # Phase shift
    msfile = eMCP['inputs']['inbase'] + ext_ms[do_ms2mms]
    logger.info('Start phase shift')
    #    t0 = datetime.datetime.utcnow()
    fixvis(vis=msfile1, outputvis=msfile, reuse=False)
    find_casa_problems()
    logger.info('Fixed {0} into {1}'.format(msfile1, msfile))
    run_listobs(msfile)
    if os.path.isdir(msfile) == True:
        emutils.rmdir(msfile1)
    else:
        logger.critical('Problem generating {}. Stopping ' \
                        'pipeline'.format(msfile))
        exit_pipeline(eMCP)
    if is_mixed_mode:
        msfile_sp = get_msfile_sp(eMCP)
        logger.info('Running phaseshift on spectral line data')
        fixvis(vis=msfile1_sp, outputvis=msfile_sp, reuse=False)
        find_casa_problems()
        logger.info('Fixed {0} into {1}'.format(msfile1_sp, msfile_sp))
        run_listobs(msfile_sp)
        if os.path.isdir(msfile_sp) == True:
            emutils.rmdir(msfile1_sp)
        else:
            logger.critical(
                'Problem generating phase shifted ms. Stopping pipeline')
            exit_pipeline(eMCP)
    logger.info('Finished phaseshift')
    msg = ''
    logger.info('End run_importfits')
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    eMCP = add_step_time('run_importfits', eMCP, msg, t0, doweblog=True)
    flag_statistics(eMCP, step='run_importfits')
    return eMCP


def fix_repeated_sources(msfile0, msfile1, datacolumn, antenna, spw_separation,
                         timeaverage, timebin, chanaverage, chanbin,
                         usewtspectrum, do_hanning, do_ms2mms):
    """ In case a source with the same name and phase centre is loaded
    with different field ID, despite being the same, we can split each
    field ID individually, and then concatenate. The concatenation will
    check names and positions and merge them if they are the same"""
    mssources = find_mssources(msfile0)
    num_sources = len(mssources.split(','))
    tmp_files = []
    for f in range(num_sources):
        outputmsfile_tmp = 'field_{0:02.0f}'.format(f)
        emutils.rmdir(outputmsfile_tmp)
        logger.info('Creating ... {0}'.format(outputmsfile_tmp))
        commands = {}
        commands['mstransform'] = {
            'vis': msfile0,
            'outputvis': outputmsfile_tmp,
            'keepflags': True,
            'datacolumn': datacolumn,
            'antenna': antenna,
            'spw': spw_separation,
            'timeaverage': timeaverage,
            'timebin': timebin,
            'chanaverage': chanaverage,
            'chanbin': chanbin,
            'usewtspectrum': usewtspectrum,
            'hanning': do_hanning,
            'createmms': do_ms2mms
        }
        mstransform(vis=msfile0,
                    outputvis=outputmsfile_tmp,
                    field=str(f),
                    keepflags=True,
                    datacolumn=datacolumn,
                    antenna=antenna,
                    spw=spw_separation,
                    timeaverage=timeaverage,
                    timebin=timebin,
                    chanaverage=chanaverage,
                    chanbin=chanbin,
                    usewtspectrum=usewtspectrum,
                    hanning=do_hanning,
                    createmms=do_ms2mms)
        find_casa_problems()
        tmp_files.append(outputmsfile_tmp)
    logger.info('Concat individual fields')
    concat(vis=tmp_files, concatvis=msfile1, timesort=True)
    find_casa_problems()
    for f in range(num_sources):
        outputmsfile_tmp = 'field_{0:02.0f}'.format(f)
        logger.info('Removing... {0}'.format(outputmsfile_tmp))
        emutils.rmdir(outputmsfile_tmp)
    return


def run_aoflagger_fields(eMCP):
    """
    This version of the autoflagger iterates through the ms within the mms structure selecting individual fields.
    It uses pre-defined strategies. The np.unique in the loop below is needed for single source files.
    The mms thinks there are many fields (one per mms). I think it is a bug from virtualconcat.
    """
    fields = eMCP['defaults']['aoflagger']['fields']
    separate_bands = eMCP['defaults']['aoflagger']['separate_bands']
    pipeline_path = eMCP['pipeline_path']
    msfile = eMCP['msinfo']['msfile']
    run_input = eMCP['defaults']['aoflagger']['run']
    if run_input == 'auto':
        logger.info('aoflagger mode "auto"')
        if eMCP['msinfo']['band'] == 'L':
            logger.info('L-band data, aoflagger will be executed')
            run_aoflagger = True
        else:
            run_aoflagger = False
            logger.info('Not L band, aoflagger will NOT be executed')
            return eMCP
    elif run_input == '':
        run_aoflagger = True

    logger.info('Start run_aoflagger_fields')
    t0 = datetime.datetime.utcnow()
    if separate_bands:
        logger.info('Bands will be processed separately.')
    elif not separate_bands:
        logger.info('Bands will be processed all together.')
    else:
        logger.warning('separate_bands can only be True or False')
        exit_pipeline(eMCP)
    # Check if aoflagger is available:
    aoflagger_available = check_command('aoflagger')
    if not aoflagger_available:
        logger.critical('aoflagger requested but not available.')
        logger.warning('Exiting pipeline.')
        exit_pipeline(eMCP)
    # Check that version is at least 2.9+
    old_aoflagger = check_aoflagger_version()
    if old_aoflagger:
        logger.critical('aoflagger version <2.9 does not work correctly.')
        logger.warning('Exiting pipeline.')
        exit_pipeline(eMCP)
    vis_fields = emutils.read_keyword(msfile, 'NAME', 'FIELD')

    fields_num = {f: i for i, f in enumerate(vis_fields)}
    if fields == 'all':
        fields = vis_fields
    else:
        fields = np.atleast_1d(fields)
    for field in np.unique(fields):
        # First, check if user has a new strategy for this field in the local folder.
        # If not, check if user has produced a new strategy for this field in the pipeline folder
        # (for typical sources, etc).
        # If not, check for default strategies for this field
        # If nothing is found, just use the default strategy
        ao_strategy = ""
        aoflagger_strategies_path = emutils.get_project_root(
        ) / "aoflagger_strategies"
        aoflagger_strategies_user = Path("./") / "user"
        aoflagger_strategies_default = aoflagger_strategies_path / "default"

        ao_strategy_path_user = aoflagger_strategies_user / "{0}.lua".format(
            field)
        ao_strategy_path_default = aoflagger_strategies_default / "{0}.lua".format(
            field)
        ao_strategy_path_default_faint = aoflagger_strategies_default / "{0}.lua".format(
            'default_faint')

        logger.info('Looking for file {0}'.format(str(ao_strategy_path_user)))
        if ao_strategy_path_user.is_file():
            ao_strategy = str(ao_strategy_path_user)
        elif ao_strategy_path_default.is_file():
            logger.info('File {0} not found'.format(
                str(ao_strategy_path_user)))
            logger.info('Looking for file {0}'.format(
                str(ao_strategy_path_default)))
            ao_strategy = str(ao_strategy_path_default)
        else:
            logger.info('Files {0} and {1} not found'.format(
                str(ao_strategy_path_user), str(ao_strategy_path_default)))
            logger.info('Looking for file {0}'.format(
                str(ao_strategy_path_default_faint)))
            ao_strategy = str(ao_strategy_path_default_faint)

        logger.info(
            'Running AOFLagger for field {0} ({1}) using strategy {2}'.format(
                field, fields_num[field], ao_strategy))
        if separate_bands:
            num_spw = len(
                vishead(msfile, mode='list',
                        listitems=['spw_name'])['spw_name'][0])
            for b in range(num_spw):
                logger.info('Processing source {0}, band {1}'.format(field, b))
                flag_command = 'aoflagger -fields {2} -bands {3} -strategy {0} {1}'.format(
                    ao_strategy, msfile, fields_num[field], b)
                os.system(flag_command)
                logger.info('Last AOFlagger command: {}'.format(flag_command))
        else:
            logger.info('Processing source {0}, all bands'.format(field))
            flag_command = 'aoflagger -fields {2} -strategy {0} {1}'.format(
                ao_strategy, msfile, fields_num[field])
            os.system(flag_command + ' | tee -a pre-cal_flag_stats.txt')
            logger.info('Last AOFlagger command: {}'.format(flag_command))

    flag_statistics(eMCP, step='flag_aoflagger')
    logger.info('End flag_aoflagger')
    msg = ''
    eMCP = add_step_time('flag_aoflagger', eMCP, msg, t0)
    return eMCP


def check_command(command):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([command], stdout=devnull,
                         stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


def check_aoflagger_version():
    logger.info('Checking AOflagger version')
    from subprocess import Popen, PIPE
    try:
        process = Popen(['aoflagger --version'], stdout=PIPE)
        (output, err) = process.communicate()
        exit_code = process.wait()
        version = output.split()[1]
        version_list = version.split('.')
        if (version_list[0] == '2') and (int(version_list[1]) < 9):
            old_aoflagger = True
        else:
            old_aoflagger = False
        logger.info('AOflagger version is {0}'.format(version))
    except:
        logger.info('AOflagger not available in this computer.')
        old_aoflagger = False
    return old_aoflagger


def find_quacktime(msinfo, s1, s2):
    separations = msinfo['separations']
    if s1 not in msinfo['sources']['mssources'].split(','):
        logger.warning('{} not in MS'.format(s1))
        separation = 'Not in MS'
        return 0
    elif s2 not in msinfo['sources']['mssources'].split(','):
        logger.warning('{} not in MS'.format(s2))
        separation = 'Not in MS'
        return 0
    else:
        try:
            separation = separations[s1 + '-' + s2]
        except KeyError:
            try:
                separation = separations[s2 + '-' + s1]
            except:
                separation = 0.0

    if separation < 1.0 * u.deg:
        quack_time = 20.
    elif 1.0 * u.deg <= separation < 2.0 * u.deg:
        quack_time = 25.
    elif 2.0 * u.deg <= separation < 3.5 * u.deg:
        quack_time = 30.
    elif separation >= 3.5 * u.deg:
        quack_time = 35.
    else:
        quack_time = 0
    return quack_time


def quack_estimating(eMCP):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    # Main calibrators quack
    standard_cal_list = [
        '1331+305', '1407+284', '0319+415', 'J1331+305', 'J1407+284',
        'J0319+415', '1331+3030', '1407+2728', '0319+4130', 'J1331+3030',
        'J1407+2728', 'J0319+4130'
    ]
    bright_cal = join_lists([
        si for si in standard_cal_list
        if si in msinfo['sources']['mssources'].split(',')
    ])
    if bright_cal != '':
        std_cal_quack = eMCP['defaults']['flag_apriori']['std_cal_quack']
        logger.info(
            'Flagging {0}s from bright calibrators'.format(std_cal_quack))
        flagdata(vis=msfile,
                 field=bright_cal,
                 mode='quack',
                 quackinterval=std_cal_quack,
                 flagbackup=False)
        find_casa_problems()
    else:
        logger.warning(
            'No main calibrators (3C84, 3C286, OQ208) found in data set')
    for s1, s2 in zip(msinfo['sources']['targets'].split(','),
                      msinfo['sources']['phscals'].split(',')):
        sources_in_ms = msinfo['sources']['mssources'].split(',')
        missing_sources = [si for si in [s1, s2] if si not in sources_in_ms]
        if not missing_sources:
            quacktime = find_quacktime(msinfo, s1, s2)
            if s1 != '':
                logger.info(
                    'Flagging first {0} sec of target {1} and phasecal {2}'.
                    format(quacktime, s1, s2))
                flagdata(vis=msfile,
                         field=','.join([s1, s2]),
                         mode='quack',
                         quackinterval=quacktime,
                         flagbackup=False)
                find_casa_problems()
        else:
            logger.warning(
                'Warning, source(s) {} not present in MS, will not flag this pair'
                .format(','.join(missing_sources)))


def read_flag_database(logfile, t0, t1):
    """ Returns a numpy array with four columns and N entries going
    from the starting time to the end time requested"""
    column_names = ['timestamp', 'antenna', 'status']
    data_all = np.genfromtxt(logfile,
                             delimiter=',',
                             dtype=None,
                             names=column_names,
                             encoding=None)
    # Filter by desired time range
    timestamps = np.array(
        [dateutil.parser.parse(d) for d in data_all['timestamp']])
    on_timerange = (timestamps > t0) & (timestamps < t1)
    data = data_all[on_timerange]
    return data


def write_flags_antenna(data, antenna, t0, t1):
    data_ant = data[data['antenna'] == antenna]
    if len(data_ant) == 0:
        print('No times found for antenna {} in this time interval'.format(
            antenna))
        return ''
    flag_commands = []
    # First a flag from t0 to the first entry
    if data_ant[0]['status'] == 'ONSOURCE':
        time_ini = t_casaformat(t0)
        time_end = t_casaformat(data_ant[0]['timestamp'])
        flag_commands.append(casa_flagcommand(antenna, time_ini, time_end))
    for i, d in enumerate(data_ant):
        try:
            if data_ant[i]['status'] == 'OFFSOURCE':  # and \
                #            data_ant[i+1]['status'] == 'ONSOURCE':
                time_ini = t_casaformat(data_ant[i]['timestamp'])
                time_end = t_casaformat(data_ant[i + 1]['timestamp'])
                flag_command = casa_flagcommand(antenna, time_ini, time_end)
                flag_commands.append(flag_command)
        except IndexError:  # Last flag until t1
            if data_ant[i]['status'] == 'OFFSOURCE':
                time_ini = t_casaformat(data_ant[i]['timestamp'])
                time_end = t_casaformat(t1)
                flag_command = casa_flagcommand(antenna, time_ini, time_end)
                flag_commands.append(flag_command)

    logger.debug(f"flag_commands: {flag_commands}")
    all_commands_ant = ''.join(flag_commands)
    return all_commands_ant


def t_casaformat(t):
    try:
        timestamp = dateutil.parser.parse(t).strftime(
            format='%Y/%m/%d/%H:%M:%S.%f')
    except:
        timestamp = t.strftime(format='%Y/%m/%d/%H:%M:%S.%f')
    return timestamp


def casa_flagcommand(antenna, time_ini, time_end):
    command = "mode='manual' antenna='{0}' timerange='{1}~{2}'\n".format(
        antenna, time_ini, time_end)
    return command


def search_observatory_flags(eMCP):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    antennas = msinfo['antennas']
    flagfile = 'observatory.flags'
    if os.path.exists(flagfile):
        logger.info('Observatory flags already present, '\
                    'will not overwrite file {}'.format(flagfile))
        logger.info('Applying observatory flags in {}'.format(flagfile))
        flagdata(vis=msfile, mode='list', inpfile=flagfile, flagbackup=False)
        find_casa_problems()
        finished_autoflag = True
    else:
        t0, t1 = get_obstime(msfile)
        logger.info('Initial obs time {0}'.format(t0))
        logger.info('Final   obs time {0}'.format(t1))
        s = socket.gethostname()
        if t0 > datetime.datetime(2019, 1, 25, 17, 15, 0) and s == 'pipeline':
            # Retrieve database (Only locally at pipeline machine JBO)
            # Only available from 2019-01-25 at 17:15:00
            logger.info('Trying to retrieve observatory flags (locally)')
            s1 = 'emproc1'
            s2 = 'ast.man.ac.uk'
            loc = '/home/emerlin/jmoldon/otcx/antenna_monitor.log'
            try:
                os.system(
                    'scp -pr {0}.{1}:{2} /pipeline1/emerlin/files/'.format(
                        s1, s2, loc))
            except:
                pass
            logfile = '/pipeline1/emerlin/files/antenna_monitor.log'
            data = read_flag_database(logfile, t0, t1)
            flag_commands = ''
            logger.info('Generating flags for this dataset')
            with open(flagfile, 'w') as f_file:
                for ant in antennas:
                    logger.debug(f"Antenna: {ant}")
                    all_commands_ant = write_flags_antenna(data, ant, t0, t1)
                    f_file.write(all_commands_ant)
            logger.info('Applying observatory flags in {}'.format(flagfile))
            flagdata(vis=msfile,
                     mode='list',
                     inpfile=flagfile,
                     flagbackup=False)
            find_casa_problems()
            finished_autoflag = True
        else:
            logger.info('Database only after 25 January 2019')
            finished_autoflag = False
    return finished_autoflag


def select_first_last_chan(msfile, spw_frac, nchan):
    first_spw = emutils.read_keyword(msfile, 'DATA_DESC_ID').min()
    last_spw = emutils.read_keyword(msfile, 'DATA_DESC_ID').max()
    chan_fact = (nchan - nchan / 512.)
    first_chan = '{0}:{1:.0f}~{2:.0f}'.format(first_spw, 0. * chan_fact,
                                              spw_frac * chan_fact)
    last_chan = '{0}:{1:.0f}~{2:.0f}'.format(last_spw,
                                             (1. - spw_frac) * chan_fact,
                                             1.0 * chan_fact)
    return first_chan, last_chan


def flagdata1_apriori(eMCP):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    sources = msinfo['sources']
    logger.info('Start flagdata1_apriori')
    t0 = datetime.datetime.utcnow()
    antennas = get_antennas(msfile)
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    # Find number of channels in MS:
    nchan = emutils.read_keyword(msfile,
                                 'CHAN_FREQ',
                                 subtable='SPECTRAL_WINDOW').shape[1]
    msg = ''
    # Remove pure zeros
    logger.info('Flagging zeros')
    flagdata(vis=msfile, mode='clip', clipzeros=True, flagbackup=False)
    find_casa_problems()
    msg += 'Clip zeros. '
    if eMCP['is_mixed_mode']:
        msfile_sp = eMCP['msinfo']['msfile_sp']
        logger.info('Flagging zeros from narrow: {0}'.format(msfile_sp))
        flagdata(vis=msfile_sp, mode='clip', clipzeros=True, flagbackup=False)
        find_casa_problems()
    # Subband edges
    channels_to_flag = '*:0~{0};{1}~{2}'.format(nchan / 128 - 1,
                                                nchan - nchan / 128, nchan - 1)
    logger.info('MS has {} channels/spw'.format(nchan))
    logger.info('Flagging edge channels {0}'.format(channels_to_flag))

    flagdata(vis=msfile, mode='manual', spw=channels_to_flag, flagbackup=False)
    find_casa_problems()
    msg += 'Subband edges {}. '.format(channels_to_flag)
    # 5% of first and last channels of MS
    spw_frac = eMCP['defaults']['flag_apriori']['border_chan_perc'] / 100.
    first_chan, last_chan = select_first_last_chan(msfile, spw_frac, nchan)
    logger.info('Flagging initial/end channels: {0} and {1}'.format(
        first_chan, last_chan))

    flagdata(vis=msfile, mode='manual', spw=first_chan, flagbackup=False)
    find_casa_problems()

    flagdata(vis=msfile, mode='manual', spw=last_chan, flagbackup=False)
    find_casa_problems()
    msg += 'channels {0} {1}. '.format(first_chan, last_chan)
    # Remove 4 sec from everywhere:
    all_quack = eMCP['defaults']['flag_apriori']['all_quack']
    logger.info('Flagging first {} seconds from all scans'.format(all_quack))

    flagdata(vis=msfile,
             field='',
             mode='quack',
             quackinterval=all_quack,
             flagbackup=False)
    find_casa_problems()
    # Flag slewing times based on observatory records
    observatory_flags = eMCP['defaults']['flag_apriori']['observatory_flags']
    do_estimated_quack = eMCP['defaults']['flag_apriori']['do_estimated_quack']
    if observatory_flags:
        logger.debug('Observatory_flags selected')
        finished_autoflag = search_observatory_flags(eMCP)
        logger.debug(f"finished_autoflag: {finished_autoflag}")
        if not finished_autoflag and do_estimated_quack == 'auto':
            logger.warning(
                'Observatory flags failed. Starting ad-hoc procedure')
            quack_estimating(eMCP)
            msg += 'Estimated quack. '
        else:
            logger.warning('observatory.flags')
            msg += 'Observatory flags: observatory.flags. '
    else:
        logger.info('No obsevatory flagging selected')
    # Flag using estimated quacking interval
    if do_estimated_quack == 'force':
        logger.info('Starting ad-hoc procedure to estimate quack flagging.')
        quack_estimating(eMCP)
        msg += 'Estimated quack. '

    # Flag Lo-Mk2
    if 'Lo' in antennas and 'Mk2' in antennas:
        if eMCP['defaults']['flag_apriori']['flag_Lo-Mk2']:
            logger.info('Flagging Lo-Mk2 baseline')

            flagdata(vis=msfile,
                     mode='manual',
                     antenna='Lo*&Mk2*',
                     flagbackup=False)
            find_casa_problems()
            msg += '. Lo-Mk2 flagged'
        else:
            logger.info('Flagging Lo-Mk2 baseline is disabled')
    flag_statistics(eMCP, step='flag_apriori')
    logger.info('End flag_apriori')
    eMCP = add_step_time('flag_apriori', eMCP, msg, t0)
    return eMCP


def log_manual_flags(inpfile):
    with open(inpfile, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
        lines = [line for line in lines if line[0] != '#']
    logger.info('Contents of {}'.format(inpfile))
    max_lines = 30
    if len(lines) > max_lines:
        logger.info('Flagfile has {0} lines. '\
                    'Only showing the first {1}'.format(len(lines),
                                                        max_lines))
    for l in lines[:max_lines]:
        logger.info(l)
    if len(lines) > max_lines:
        logger.info('...')
    if len(lines) > 0:
        are_there_flags = True
    else:
        are_there_flags = False
    return are_there_flags


def flagdata_manual(eMCP, run_name='flag_manual'):
    logger.info(line0)
    msfile = eMCP['msinfo']['msfile']
    if run_name == 'flag_manual':
        inpfile = './manual.flags'
    elif run_name == 'flag_manual_avg':
        inpfile = './manual_avg.flags'
    else:
        logger.warning('Wrong run_name specified')
        inpfile = ''
    logger.info('Start {}'.format(run_name))
    t0 = datetime.datetime.utcnow()
    if os.path.isfile(inpfile) == True:
        logger.info('Applying manual flags from file: {0}'.format(inpfile))
        are_there_flags = log_manual_flags(inpfile)
        if are_there_flags:
            flagdata(vis=msfile,
                     mode='list',
                     inpfile=inpfile,
                     flagbackup=False)
            find_casa_problems()
        else:
            logger.warning('Flagfile is empty')
        msg = '{0}. '.format(inpfile)
    else:
        logger.info('No manual flag file selected: {0}'.format(inpfile))
        msg = 'No flagging file. '
    if run_name == 'flag_manual_avg':
        # Compute Lo_dropouts if needed and flag data
        eMCP, msg_i = flag_Lo_dropouts(eMCP)
        msg += msg_i
        # Apply manual_narrow.flags if available
        inpfile_narrow = './manual_narrow.flags'
        if eMCP['is_mixed_mode'] and os.path.isfile(inpfile_narrow) == True:
            logger.info(
                'Applying manual flags from file: {0}'.format(inpfile_narrow))
            are_there_flags = log_manual_flags(inpfile_narrow)
            if are_there_flags:
                flagdata(vis=eMCP['msinfo']['msfile_sp'],
                         mode='list',
                         inpfile=inpfile_narrow,
                         flagbackup=False)
                find_casa_problems()
            else:
                logger.warning('Narrow (sp) flagfile is empty')
            msg += 'Narrow (sp): {0}'.format(inpfile_narrow)
    flag_statistics(eMCP, step=run_name)
    logger.info('End {}'.format(run_name))
    eM = add_step_time(run_name, eMCP, msg, t0)
    return eMCP


def flagdata_tfcrop(eMCP, defaults):
    logger.info(line0)
    if defaults == 'flag_target':
        t0 = datetime.datetime.utcnow()
        logger.info('Start flag_target')
    msinfo = eMCP['msinfo']
    msfile = eMCP['msinfo']['msfile']
    tfcrop = eMCP['defaults'][defaults]['tfcrop']
    logger.info("Running flagdata, mode = '{}'".format(tfcrop['mode']))
    logger.info("correlation='{}'".format(tfcrop['correlation']))
    logger.info("ntime='{0}', combinescans={1}, datacolumn='{2}'".format(
        tfcrop['ntime'], tfcrop['combinescans'], tfcrop['datacolumn']))
    logger.info(
        "winsize={0}, timecutoff={1}, freqcutoff={2}, maxnpieces={3}".format(
            tfcrop['winsize'], tfcrop['timecutoff'], tfcrop['freqcutoff'],
            tfcrop['maxnpieces']))
    logger.info("usewindowstats='{0}', halfwin={1}, extendflags={2}".format(
        tfcrop['uwstats'], tfcrop['halfwin'], tfcrop['extendflags']))

    flagdata(vis=msfile,
             mode='tfcrop',
             field=eMCP['msinfo']['sources'][tfcrop['sources']],
             antenna=tfcrop['antenna'],
             scan=tfcrop['scan'],
             spw=tfcrop['spw'],
             correlation=tfcrop['correlation'],
             ntime=tfcrop['ntime'],
             combinescans=tfcrop['combinescans'],
             datacolumn=tfcrop['datacolumn'],
             winsize=tfcrop['winsize'],
             timecutoff=tfcrop['timecutoff'],
             freqcutoff=tfcrop['freqcutoff'],
             maxnpieces=tfcrop['maxnpieces'],
             usewindowstats=tfcrop['uwstats'],
             halfwin=tfcrop['halfwin'],
             extendflags=tfcrop['extendflags'],
             action=tfcrop['action'],
             display=tfcrop['display'],
             flagbackup=tfcrop['flagbackup'])
    find_casa_problems()
    if defaults == 'flag_target':
        flag_statistics(eMCP, step='flag_target')
        msg = "mode={0}, maxnpieces={1}, timecutoff={2}, freqcutoff={3}".format(
            'tfcrop', tfcrop['maxnpieces'], tfcrop['timecutoff'],
            tfcrop['freqcutoff'])
        eMCP = add_step_time('flag_target', eMCP, msg, t0)
        logger.info('End flag_target')
    return eMCP


def flagdata_rflag(eMCP, defaults):
    logger.info(line0)
    if defaults == 'flag_target':
        t0 = datetime.datetime.utcnow()
        logger.info('Start flag_target')
    msinfo = eMCP['msinfo']
    msfile = eMCP['msinfo']['msfile']
    rflag = eMCP['defaults'][defaults]['rflag']
    logger.info("Running flagdata, mode = '{}'".format(rflag['mode']))
    logger.info("ntime='{0}', combinescans={1}, datacolumn='{2}'".format(
        rflag['ntime'], rflag['combinescans'], rflag['datacolumn']))
    logger.info("timedevscale = {0}, freqdevscale = {1}".format(
        rflag['timedevscale'], rflag['freqdevscale']))

    flagdata(vis=msfile,
             mode=rflag['mode'],
             field=eMCP['msinfo']['sources'][rflag['sources']],
             antenna=rflag['antenna'],
             scan=rflag['scan'],
             spw=rflag['spw'],
             correlation=rflag['correlation'],
             ntime=rflag['ntime'],
             combinescans=rflag['combinescans'],
             datacolumn=rflag['datacolumn'],
             timedevscale=rflag['timedevscale'],
             freqdevscale=rflag['freqdevscale'],
             action=rflag['action'],
             display=rflag['display'],
             flagbackup=rflag['flagbackup'])
    find_casa_problems()
    if defaults == 'flag_target':
        flag_statistics(eMCP, step='flag_target')
        msg = "mode={0}, ntime={1}, timedevscale={2}, freqdevscale={3}".format(
            'rflag', rflag['ntime'], rflag['timedevscale'],
            rflag['freqdevscale'])
        eMCP = add_step_time('flag_target', eMCP, msg, t0)
        logger.info('End flag_target')
    return eMCP


def run_flag_target(eMCP):
    mode = eMCP['defaults']['flag_target']['mode_to_run']
    if mode == 'rflag':
        eMCP = flagdata_rflag(eMCP, 'flag_target')
    elif mode == 'tfcrop':
        eMCP = flagdata_tfcrop(eMCP, 'flag_target')


def define_refant(eMCP):
    msfile = eMCP['msinfo']['msfile']
    default_refant = eMCP['defaults']['global']['refant']
    if default_refant == 'compute':
        recompute = True
        logger.info('Forcing recompute of refant')
    elif default_refant == '':
        try:
            refant = eMCP['msinfo']['refant']
            if refant == '':
                logger.info('No refant specified. Will compute optimal')
                recompute = True
            else:
                logger.info(
                    'Refant already in eMCP["msinfo"], will not recompute')
                recompute = False
        except:
            logger.info('No refant specified. Will compute optimal')
            recompute = True
    else:
        recompute = False
        logger.info('Refant from inputs file')
        refant = default_refant
    if recompute:
        logger.info('Recomputing best reference antenna')
        refant0 = default_refant
        refant_user = refant0.replace(' ', '').split(',')
        antennas = get_antennas(msfile)
        refant_in_ms = (np.array([ri in antennas for ri in refant_user])).all()
        if not refant_in_ms:
            if refant0 != '' and default_refant != 'compute':
                logger.warning(
                    'Selected reference antenna(s) {0} not in MS! User selection will be ignored'
                    .format(refant0))
            # Finding best antennas for refant
            logger.info(
                'To avoid this, set a reference antenna in the default parameters file.'
            )
            refant = find_refant(msfile,
                                 field=eMCP['msinfo']['sources']['calsources'])
        else:
            refant = ','.join(refant_user)
    logger.info('Refant in eMCP: {}'.format(refant))
    return refant


def find_refant(msfile, field):
    # Find phase solutions per scan:
    tablename = os.path.join(calib_dir, 'find_refant.phase')
    logger.info('Phase calibration per scan to count number of good solutions')

    gaincal(vis=msfile,
            caltable=tablename,
            field=field,
            refantmode='flex',
            solint='inf',
            minblperant=2,
            gaintype='G',
            calmode='p')
    find_casa_problems()
    # Read solutions (phases):
    antenna_names = emutils.read_keyword(tablename, 'NAME', subtable='ANTENNA')
    antenna_ids = emutils.read_keyword(tablename, 'ANTENNA1')
    flags = emutils.read_keyword(tablename, 'FLAG')
    phases = np.angle(emutils.read_keyword(tablename, 'CPARAM'))
    snrs = emutils.read_keyword(tablename, 'SNR')

    # Analyse number of good solutions:
    good_frac = []
    good_snrs = []
    for i, ant_id in enumerate(np.unique(antenna_ids)):
        cond = antenna_ids == ant_id
        #t = times[cond]
        f = flags[:, 0, 0][cond]
        p = phases[:, 0, 0][cond]
        snr = snrs[:, 0, 0][cond]
        frac = 1.0 * np.count_nonzero(~f) / len(f) * 100.
        snr_mean = np.nanmean(snr[~f])
        good_frac.append(frac)
        good_snrs.append(snr_mean)
    sort_idx = np.argsort(good_frac)[::-1]
    logger.info('Antennas sorted by % of good solutions:')
    for i in sort_idx:
        logger.info('{0:3}: {1:4.1f}, <SNR> = {2:4.1f}'.format(
            antenna_names[i], good_frac[i], good_snrs[i]))
    if good_frac[sort_idx[0]] < 90:
        logger.warning(
            'Small fraction of good solutions with selected refant!')
        logger.warning('Please inspect antennas to select optimal refant')
        logger.warning(
            'You may want to use refantmode= flex" in default_params')
    pref_ant = np.array(antenna_names)[sort_idx]
    if 'Lo' in antenna_names:
        priorities = ['Pi', 'Da', 'Kn', 'De', 'Cm']
    else:
        priorities = ['Mk2', 'Pi', 'Da', 'Kn', 'Cm', 'De']
    refant = ','.join([a for a in pref_ant if a in priorities])
    return refant


def saveflagstatus(eMCP):
    msinfo = eMCP['msinfo']
    logger.info('Starting saveflagstatus')
    t0 = datetime.datetime.utcnow()
    logger.info(
        'Saving current flagging status to versionname=\'initialize_flags\'')

    flagmanager(
        msinfo['msfile'],
        mode='save',
        versionname='initialize_flags',
        comment=
        'Restore this version to restart calibration without the flags produced by the calibration',
        merge='replace')
    find_casa_problems()
    if eMCP['is_mixed_mode']:
        flagmanager(
            msinfo['msfile_sp'],
            mode='save',
            versionname='initialize_flags',
            comment=
            'Restore this version to restart calibration without the flags produced by the calibration',
            merge='replace')
        find_casa_problems()
    msg = 'versionname=initialize_flags'
    eMCP = add_step_time('save_flags', eMCP, msg, t0)
    return eMCP


def restoreflagstatus(eMCP):
    msinfo = eMCP['msinfo']
    t0 = datetime.datetime.utcnow()
    logger.info('Starting restoreflagstatus')
    logger.info(
        'Restoring flagging status in versionname=\'initialize_flags\'')
    flagmanager(msinfo['msfile'],
                mode='restore',
                versionname='initialize_flags',
                merge='replace')
    find_casa_problems()
    if eMCP['is_mixed_mode']:
        flagmanager(msinfo['msfile_sp'],
                    mode='restore',
                    versionname='initialize_flags',
                    merge='replace')
        find_casa_problems()

    flag_statistics(eMCP, step='restore_flags')
    msg = 'versionname=initialize_flags'
    eMCP = add_step_time('restore_flags', eMCP, msg, t0)
    return eMCP


def check_sources_in_ms(eMCP):
    msfile = eMCP['msinfo']['msfile']
    sources = eMCP['msinfo']['sources']
    mssources = find_mssources(msfile)
    targets = sources['targets']
    phscals = sources['phscals']
    # Check that targets/phscals are not empty:
    if targets == '' or phscals == '':
        logger.critical('Targets or phase calibrators not specified')
        logger.warning('Stopping pipeline at this step')
        exit_pipeline(eMCP)
    # Check that all sources are in the MS:
    sources_not_in_msfile = [
        s for s in sources['allsources'].split(',')
        if s not in mssources.split(',')
    ]
    if len(sources_not_in_msfile) > 0:
        fields = mssources
        logger.critical('Fields {} not present in MS but listed in ' \
                       'inputs file.'.format(','.join(sources_not_in_msfile)))
        logger.warning('Stopping pipeline at this step')
        exit_pipeline(eMCP)


def check_table_exists(caltables, tablename):
    logger.debug('Try existence: {0}'.format(caltables[tablename]['table']))
    try:
        if os.path.isdir(caltables[tablename]['table']):
            logger.debug('OK')
            it_exists = True
        else:
            logger.debug('Not OK')
            it_exists = False
    except:
        it_exists = False
    if not it_exists:
        logger.critical('Calibration table {} required but ' \
                        'not available'.format(tablename))
        logger.warning('Stopping pipeline at this step')
        exit_pipeline()


### Run CASA calibration functions
def list_sources_out(eMCP):
    mssources = find_mssources(eMCP['msinfo']['msfile']).split(',')
    allsources = eMCP['msinfo']['sources']['allsources'].split(',')
    ignored_sources = [s for s in mssources if s not in allsources]
    if ignored_sources:
        logger.info('Sources ignored for not being in inputs file: {0}'.format(
            ','.join(ignored_sources)))


def run_average(eMCP):
    logger.info(line0)
    logger.info('Start average')
    t0 = datetime.datetime.utcnow()
    chanbin = eMCP['defaults']['average']['chanbin']
    timebin = eMCP['defaults']['average']['timebin']
    msfile = eMCP['msinfo']['msfile']
    if timebin == '1s':
        timeaverage = False
    else:
        timeaverage = True
    if chanbin == 1:
        chanaverage = False
    else:
        chanaverage = True
    datacolumn = eMCP['defaults']['average']['datacolumn']
    scan = eMCP['defaults']['average']['scan']
    antenna = eMCP['defaults']['average']['antenna']
    timerange = eMCP['defaults']['average']['timerange']
    fields = eMCP['msinfo']['sources']['allsources']
    name, exte = os.path.splitext(msfile)
    outputmsfile = name + '_avg' + exte
    emutils.rmdir(outputmsfile)
    emutils.rmdir(outputmsfile + '.flagversions')
    logger.info('Input MS: {0}'.format(msfile))
    logger.info('Output MS: {0}'.format(outputmsfile))
    logger.info('chanbin={0}, timebin={1}'.format(chanbin, timebin))
    logger.info('Fields: {0}'.format(fields))
    list_sources_out(eMCP)
    logger.info('Data column: {0}'.format(datacolumn))

    mstransform(vis=msfile,
                outputvis=outputmsfile,
                field=fields,
                timeaverage=timeaverage,
                chanaverage=chanaverage,
                timerange=timerange,
                scan=scan,
                antenna=antenna,
                timebin=timebin,
                chanbin=chanbin,
                datacolumn=datacolumn,
                keepflags=True)
    find_casa_problems()
    run_listobs(outputmsfile)
    msg = 'chanbin={0}, timebin={1}, datacolumn={2}'.format(
        chanbin, timebin, datacolumn)
    logger.info('End average')

    if eMCP['defaults']['average']['shift_phasecenter']:
        shift_all_positions(eMCP)
        msg += '. Shifted position'
    eMCP = add_step_time('average', eMCP, msg, t0)
    return eMCP


def load_3C286_model(eMCP):
    fluxcal = eMCP['msinfo']['sources']['fluxcal']
    msfile = eMCP['msinfo']['msfile']
    init_models = eMCP['defaults']['init_models']
    models_path = emutils.get_project_root() / "calibrator_models"
    # Check dataset frequency:
    band = check_band(eMCP, msfile)
    if band == 'C':
        model_file_path = models_path / '3C286_C.clean.model.tt0'
        logger.info('Dataset is band C. Using C band model of 3C286')
    elif band == 'L':
        model_file_path = models_path / '1331+305.clean.model.tt0'
        logger.info('Dataset is band L. Using L band model of 3C286')
    else:
        logger.warning('No 3C286 model available!')
        model_file_path = ""
    model_file_path = str(model_file_path)
    logger.info('Initializing model for 3C286')
    logger.info('Model {0}'.format(model_file_path))

    setjy(vis=msfile,
          field=fluxcal,
          model=model_file_path,
          scalebychan=True,
          usescratch=True)
    find_casa_problems()
    if eMCP['is_mixed_mode']:
        msfile_sp = eMCP['msinfo']['msfile_sp']
        logger.info('Initializing model for 3C286 for {0}'.format(msfile_sp))

        setjy(vis=msfile_sp,
              field=fluxcal,
              model=model_file_path,
              scalebychan=True,
              usescratch=True)
        find_casa_problems()
    return


def run_initialize_models(eMCP):
    logger.info(line0)
    logger.info('Start init_models')
    t0 = datetime.datetime.utcnow()
    init_models = eMCP['defaults']['init_models']
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    msfile = eMCP['msinfo']['msfile']
    fluxcal = eMCP['msinfo']['sources']['fluxcal']
    logger.info('Resetting corrected column with clearcal')

    clearcal(vis=msfile)
    find_casa_problems()
    logger.info('Resetting weights')
    logger.info('wtmode={0}, dowtsp={1}'.format(init_models['wtmode'],
                                                init_models['dowtsp']))

    initweights(vis=msfile,
                wtmode=init_models['wtmode'],
                dowtsp=init_models['dowtsp'])
    find_casa_problems()

    if eMCP['is_mixed_mode']:
        msfile_sp = eMCP['msinfo']['msfile_sp']
        logger.info('Resetting from {0}'.format(msfile_sp))
        clearcal(vis=msfile_sp)
        find_casa_problems()
    logger.info('Deleting model of all sources')
    delmod(vis=msfile, otf=True, scr=True)  #scr to delete MODEL column
    find_casa_problems()
    if eMCP['is_mixed_mode']:
        msfile_sp = eMCP['msinfo']['msfile_sp']
        logger.info('Deleting model from {0}'.format(msfile_sp))

        delmod(vis=msfile_sp, otf=True, scr=True)
        find_casa_problems()
    # Set flux density of flux calibrator:
    if fluxcal in ['1331+305', '1331+3030', 'J1331+305', 'J1331+3030']:
        load_3C286_model(eMCP)
    elif fluxcal == '':
        logger.warning('No flux calibrator selected')
        pass
    else:
        logger.warning('Using a non-standard flux calibrator.')
        if init_models['manual_fluxcal']:
            logger.info('Manual flux density for: {0}'.format(fluxcal))
            logger.info('Flux {0}, spix {1}, reffreq {2}'.format(
                init_models['fluxcal_flux'], init_models['fluxcal_spix'],
                init_models['fluxcal_reffreq']))

            setjy(vis=msfile,
                  field=fluxcal,
                  standard='manual',
                  fluxdensity=init_models['fluxcal_flux'],
                  spix=init_models['fluxcal_spix'],
                  reffreq=init_models['fluxcal_reffreq'],
                  usescratch=True)
            find_casa_problems()
        else:
            logger.warning('Cannot initialize non-standard flux calibrator'\
                           ' {}'.format(fluxcal))
            logger.warning('Please provide manual flux density, spectral '\
                           'index and reffreq in defaults file')
            exit_pipeline(eMCP)
    logger.info('End init_models')
    msg = ''
    eMCP = add_step_time('init_models', eMCP, msg, t0)
    return eMCP


def initialize_cal_dict(eMCP):
    # All the calibration steps will be saved in the dictionary caltables.pkl
    # located in the calib directory. If it does not exist a new one is created.
    # Reference antenna(s)
    try:
        msinfo = eMCP['msinfo']
    except:
        logger.warning('No msinfo found')
        sys.exit()
    any_calsteps = [
        'flag_manual_avg', 'bandpass', 'initial_gaincal', 'fluxscale',
        'bandpass_final', 'gaincal_final', 'applycal_all'
    ]
    if np.array([eMCP['input_steps'][cal] > 0 for cal in any_calsteps]).any():
        try:
            caltables = load_obj(os.path.join(calib_dir, 'caltables.pkl'))
            logger.info('Loaded previous calibration tables from: {0}'.format(
                calib_dir + 'caltables.pkl'))
        except:
            caltables = {}
            caltables['inbase'] = eMCP['inputs']['inbase']
            caltables['plots_dir'] = plots_dir
            caltables['calib_dir'] = calib_dir
            caltables['num_spw'] = msinfo['num_spw']
            logger.info(
                'New caltables dictionary created. Saved to: {0}'.format(
                    calib_dir + 'caltables.pkl'))
        # Refant
        eMCP['msinfo']['refant'] = define_refant(eMCP)
        save_obj(eMCP, os.path.join(info_dir, 'eMCP_info.pkl'))
        caltables['refant'] = eMCP['msinfo']['refant']
        caltables['Lo_dropout_scans'] = eMCP['msinfo']['Lo_dropout_scans']
        caltables['refantmode'] = eMCP['defaults']['global']['refantmode']
        save_obj(caltables, os.path.join(calib_dir, 'caltables.pkl'))
        return caltables


def read_gaincal_solutions():
    """ Compile solution statistics from the CASA log
        Needs to be executed right after gaincal """
    # TODO: read casa logs to find solution statistics
    pass


def read_applycal_flags():
    """ Compile flag statistics from the CASA log
        Needs to be executed right after applycal """
    #'#    casalog_file = casa_log.logfile()
    # TODO: read casa logs to find flag statistics
    pass


def run_gaincal(msfile, caltables, caltable_name):
    logger.info(line0)
    emutils.rmdir(caltables[caltable_name]['table'])
    logger.info('Running gaincal to generate: {0}'.format(
        caltables[caltable_name]['name']))
    logger.info('Field(s) = {0}, gaintype = {1}, calmode = {2}'.format(
        caltables[caltable_name]['field'],
        caltables[caltable_name]['gaintype'],
        caltables[caltable_name]['calmode']))
    logger.info('solint = {0}, spw = {1},  combine = {2}'.format(
        caltables[caltable_name]['solint'], caltables[caltable_name]['spw'],
        caltables[caltable_name]['combine']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    spwmap = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(
        caltables[caltable_name]['table']))
    # Run CASA task gaincal
    gaincal(vis=msfile,
            caltable=caltables[caltable_name]['table'],
            field=caltables[caltable_name]['field'],
            gaintype=caltables[caltable_name]['gaintype'],
            calmode=caltables[caltable_name]['calmode'],
            solint=caltables[caltable_name]['solint'],
            combine=caltables[caltable_name]['combine'],
            spw=caltables[caltable_name]['spw'],
            refant=caltables['refant'],
            refantmode=caltables['refantmode'],
            gaintable=gaintable,
            gainfield=gainfield,
            interp=interp,
            spwmap=spwmap,
            minblperant=caltables[caltable_name]['minblperant'],
            minsnr=caltables[caltable_name]['minsnr'])
    find_casa_problems()
    read_gaincal_solutions()
    logger.info('caltable {0} in {1}'.format(
        caltables[caltable_name]['name'], caltables[caltable_name]['table']))


def convert_spwmap(spwmap, spwmap_sp):
    if spwmap == []:
        spwmap_out = spwmap_sp
    else:
        spwmap_out = [spwmap[s] for s in spwmap_sp]
    logger.debug('Wide spwmap {}'.format(spwmap))
    logger.debug('Narrow spwmap {}'.format(spwmap_sp))
    logger.debug('To use spwmap {}'.format(spwmap_out))
    return spwmap_out


def run_gaincal_narrow(msfile_sp, caltables, caltable_name, spwmap_sp):
    logger.info(line0)
    emutils.rmdir(caltables[caltable_name]['table'])
    logger.info('Running gaincal to generate: {0}'.format(
        caltables[caltable_name]['name']))
    logger.info('Field(s) = {0}, gaintype = {1}, calmode = {2}'.format(
        caltables[caltable_name]['field'],
        caltables[caltable_name]['gaintype'],
        caltables[caltable_name]['calmode']))
    logger.info('solint = {0}, spw = {1},  combine = {2}'.format(
        caltables[caltable_name]['solint'], caltables[caltable_name]['spw'],
        caltables[caltable_name]['combine']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    spwmap = []
    for p in previous_cal:
        logger.debug('Previous cal: {}'.format(p))
        spwmap.append(convert_spwmap(caltables[p]['spwmap'], spwmap_sp))
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(
        caltables[caltable_name]['table']))
    # Run CASA task gaincal
    gaincal(vis=msfile_sp,
            caltable=caltables[caltable_name]['table'],
            field=caltables[caltable_name]['field'],
            gaintype=caltables[caltable_name]['gaintype'],
            calmode=caltables[caltable_name]['calmode'],
            solint=caltables[caltable_name]['solint'],
            combine=caltables[caltable_name]['combine'],
            spw=caltables[caltable_name]['spw'],
            refant=caltables['refant'],
            refantmode=caltables['refantmode'],
            gaintable=gaintable,
            gainfield=gainfield,
            interp=interp,
            spwmap=spwmap,
            minblperant=caltables[caltable_name]['minblperant'],
            minsnr=caltables[caltable_name]['minsnr'])
    find_casa_problems()
    read_gaincal_solutions()
    logger.info('caltable {0} in {1}'.format(
        caltables[caltable_name]['name'], caltables[caltable_name]['table']))


def run_bandpass(msfile, caltables, caltable_name, minblperant=3, minsnr=2):
    logger.info(line0)
    emutils.rmdir(caltables[caltable_name]['table'])
    logger.info('Running bandpass to generate: {0}'.format(
        caltables[caltable_name]['name']))
    logger.info(
        'Field(s) {0}, solint = {1}, spw = {2}, combine = {3}, solnorm = {4}'.
        format(caltables[caltable_name]['field'],
               caltables[caltable_name]['solint'],
               caltables[caltable_name]['spw'],
               caltables[caltable_name]['combine'],
               caltables[caltable_name]['solnorm']))
    logger.debug('uvrange = {0}'.format(caltables[caltable_name]['uvrange']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    spwmap = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(
        caltables[caltable_name]['table']))
    # Run CASA task bandpass

    bandpass(vis=msfile,
             caltable=caltables[caltable_name]['table'],
             field=caltables[caltable_name]['field'],
             solint=caltables[caltable_name]['solint'],
             combine=caltables[caltable_name]['combine'],
             spw=caltables[caltable_name]['spw'],
             solnorm=caltables[caltable_name]['solnorm'],
             uvrange=caltables[caltable_name]['uvrange'],
             fillgaps=caltables[caltable_name]['fillgaps'],
             refant=caltables['refant'],
             gaintable=gaintable,
             gainfield=gainfield,
             interp=interp,
             spwmap=spwmap,
             minblperant=minblperant,
             minsnr=minsnr)
    find_casa_problems()
    logger.info('caltable {0} in {1}'.format(
        caltables[caltable_name]['name'], caltables[caltable_name]['table']))


def run_bandpass_narrow(msfile_sp,
                        caltables,
                        caltable_name,
                        spwmap_sp,
                        minblperant=3,
                        minsnr=2):
    logger.info(line0)
    emutils.rmdir(caltables[caltable_name]['table'])
    logger.info('Running bandpass to generate: {0}'.format(
        caltables[caltable_name]['name']))
    logger.info(
        'Field(s) {0}, solint = {1}, spw = {2}, combine = {3}, solnorm = {4}'.
        format(caltables[caltable_name]['field'],
               caltables[caltable_name]['solint'],
               caltables[caltable_name]['spw'],
               caltables[caltable_name]['combine'],
               caltables[caltable_name]['solnorm']))
    logger.debug('uvrange = {0}'.format(caltables[caltable_name]['uvrange']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    spwmap = []
    for p in previous_cal:
        if p == 'narrow_p_offset.G3':
            spwmap.append([])
        else:
            logger.debug('Previous cal: {}'.format(p))
            spwmap.append(convert_spwmap(caltables[p]['spwmap'], spwmap_sp))
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(
        caltables[caltable_name]['table']))
    # Run CASA task bandpass
    bandpass(vis=msfile_sp,
             caltable=caltables[caltable_name]['table'],
             field=caltables[caltable_name]['field'],
             solint=caltables[caltable_name]['solint'],
             combine=caltables[caltable_name]['combine'],
             spw=caltables[caltable_name]['spw'],
             solnorm=caltables[caltable_name]['solnorm'],
             uvrange=caltables[caltable_name]['uvrange'],
             fillgaps=caltables[caltable_name]['fillgaps'],
             refant=caltables['refant'],
             gaintable=gaintable,
             gainfield=gainfield,
             interp=interp,
             spwmap=spwmap,
             minblperant=minblperant,
             minsnr=minsnr)
    find_casa_problems()
    logger.info('caltable {0} in {1}'.format(
        caltables[caltable_name]['name'], caltables[caltable_name]['table']))


def smooth_caltable(msfile,
                    tablein,
                    plotdir,
                    caltable='',
                    field='',
                    smoothtype='median',
                    smoothtime=120.):
    logger.info(
        'Smoothing table: {0}, field {1}, smoothtype {2}, smoothtime {3}'.
        format(tablein, field, smoothtype, smoothtime))
    basename = os.path.basename(tablein)
    caltableplot_phs = plotdir + basename + '_phs.png'
    caltableplot_amp = plotdir + basename + '_amp.png'
    if caltable == '':
        os.system('mv {0} {1}'.format(
            caltableplot_phs, plotdir + basename + '_phs_pre_smooth.png'))
        os.system('mv {0} {1}'.format(
            caltableplot_amp, plotdir + basename + '_amp_pre_smooth.png'))
    smoothcal(vis=msfile,
              tablein=tablein,
              caltable=tablein + 'smooth',
              field='',
              smoothtype='median',
              smoothtime=60 * 20.)
    logger.info('Pre-smoothing table saved to: {0}'.format(tablein +
                                                           '_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein, tablein + '_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein + 'smooth', tablein))
    return


def run_applycal(eMCP, caltables, step, insources=''):
    logger.info(line0)
    logger.info('Start applycal')
    logger.info('Applying calibration up to step: {}'.format(step))
    previous_cal = eMCP['defaults'][step]['apply_calibrators']
    previous_cal_targets = eMCP['defaults'][step]['apply_targets']
    applymode = eMCP['defaults']['global']['applymode']
    msfile = eMCP['msinfo']['msfile']
    sources = eMCP['msinfo']['sources']
    # 1 correct non-target sources:
    logger.info('Applying calibration to calibrator sources')
    if insources == '':
        fields = sources['calsources']
    else:
        fields = insources
    logger.info('Fields: {0}'.format(fields))
    # Check if tables exist:
    for table_i in previous_cal:
        check_table_exists(caltables, table_i)
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    spwmap = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]

    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    applycal(vis=msfile,
             field=fields,
             gaintable=gaintable,
             gainfield=gainfield,
             interp=interp,
             spwmap=spwmap,
             applymode=applymode,
             flagbackup=False)
    find_casa_problems()
    read_applycal_flags()
    applycal_dict = {}
    for source in fields.split(','):
        applycal_dict[source] = {}
        for i, p in enumerate(previous_cal):
            applycal_dict[source][p] = [gainfield[i], spwmap[i], interp[i]]

    # 2 correct targets
    if previous_cal_targets:
        logger.info('Applying calibration to target sources')
        logger.info('Target fields: {0}'.format(sources['targets']))
        if previous_cal_targets == '':
            previous_cal_targets = previous_cal
        for i, s in enumerate(sources['targets'].split(',')):
            phscal = sources['phscals'].split(',')[i]
            if s != '' and s != phscal:
                # Check if tables exist:
                for table_i in previous_cal_targets:
                    check_table_exists(caltables, table_i)
                # Previous calibration
                gaintable = [
                    caltables[p]['table'] for p in previous_cal_targets
                ]
                interp = [caltables[p]['interp'] for p in previous_cal_targets]
                spwmap = [caltables[p]['spwmap'] for p in previous_cal_targets]
                gainfield = [
                    caltables[p]['gainfield']
                    if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))
                    < 2 else phscal for p in previous_cal_targets
                ]
                logger.info('Target: {0}. Phase calibrator: {1}'.format(
                    s, phscal))
                logger.info('Previous calibration applied: {0}'.format(
                    str(previous_cal_targets)))
                logger.info('Previous calibration gainfield: {0}'.format(
                    str(gainfield)))
                logger.info('Previous calibration spwmap: {0}'.format(
                    str(spwmap)))
                logger.info('Previous calibration interp: {0}'.format(
                    str(interp)))

                applycal(vis=msfile,
                         field=s,
                         gaintable=gaintable,
                         gainfield=gainfield,
                         interp=interp,
                         spwmap=spwmap,
                         applymode=applymode,
                         flagbackup=False)
                find_casa_problems()
                read_applycal_flags()
                applycal_dict[s] = {}
                for j, p in enumerate(previous_cal_targets):
                    applycal_dict[s][p] = [gainfield[j], spwmap[j], interp[j]]
            else:
                logger.warning('Source {} is not phase-referenced'.format(s))
    else:
        logger.info('Not applying calibration to target sources at this stage')

    eMCP['applycal_dict'] = applycal_dict
    logger.info('End applycal')


def run_applycal_narrow(eMCP, caltables, step, spwmap_sp, insources=''):
    logger.info(line0)
    logger.info('Start applycal')
    logger.info('Applying narrow band calibration up to step: {}'.format(step))
    previous_cal = eMCP['defaults'][step]['apply_narrow_calibrators']
    previous_cal_targets = eMCP['defaults'][step]['apply_narrow_targets']
    msfile = eMCP['msinfo']['msfile']
    msfile_sp = eMCP['msinfo']['msfile_sp']
    sources = eMCP['msinfo']['sources']
    # 1 correct non-target sources:
    logger.info('Applying calibration to calibrator sources')
    if insources == '':
        fields = sources['calsources']
    else:
        fields = insources
    logger.info('Fields: {0}'.format(fields))
    # Check if tables exist:
    for table_i in previous_cal:
        check_table_exists(caltables, table_i)
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    spwmap = []
    for p in previous_cal:
        if p == 'narrow_p_offset.G3' or p == 'narrow_bpcal.BP2':
            spwmap.append([])
        else:
            logger.debug('Previous cal: {}'.format(p))
            spwmap.append(convert_spwmap(caltables[p]['spwmap'], spwmap_sp))

    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))

    applycal(vis=msfile_sp,
             field=fields,
             gaintable=gaintable,
             gainfield=gainfield,
             interp=interp,
             spwmap=spwmap,
             flagbackup=False)
    find_casa_problems()
    read_applycal_flags()
    applycal_dict_sp = {}
    for source in fields.split(','):
        applycal_dict_sp[source] = {}
        for i, p in enumerate(previous_cal):
            applycal_dict_sp[source][p] = [gainfield[i], spwmap[i], interp[i]]

    # 2 correct targets
    if previous_cal_targets:
        logger.info('Applying calibration to target sources')
        logger.info('Target fields: {0}'.format(sources['targets']))
        if previous_cal_targets == '':
            previous_cal_targets = previous_cal
        for i, s in enumerate(sources['targets'].split(',')):
            phscal = sources['phscals'].split(',')[i]
            if s != '' and s != phscal:
                # Check if tables exist:
                for table_i in previous_cal_targets:
                    check_table_exists(caltables, table_i)
                # Previous calibration
                gaintable = [
                    caltables[p]['table'] for p in previous_cal_targets
                ]
                interp = [caltables[p]['interp'] for p in previous_cal_targets]
                gainfield = [
                    caltables[p]['gainfield']
                    if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))
                    < 2 else phscal for p in previous_cal_targets
                ]
                spwmap = []
                for p in previous_cal_targets:
                    if p == 'narrow_p_offset.G3' or p == 'narrow_bpcal.BP2':
                        spwmap.append([])
                    else:
                        logger.debug('Previous cal: {}'.format(p))
                        spwmap.append(
                            convert_spwmap(caltables[p]['spwmap'], spwmap_sp))
                logger.info('Target: {0}. Phase calibrator: {1}'.format(
                    s, phscal))
                logger.info('Previous calibration applied: {0}'.format(
                    str(previous_cal_targets)))
                logger.info('Previous calibration gainfield: {0}'.format(
                    str(gainfield)))
                logger.info('Previous calibration spwmap: {0}'.format(
                    str(spwmap)))
                logger.info('Previous calibration interp: {0}'.format(
                    str(interp)))

                applycal(vis=msfile_sp,
                         field=s,
                         gaintable=gaintable,
                         gainfield=gainfield,
                         interp=interp,
                         spwmap=spwmap,
                         flagbackup=False)
                find_casa_problems()
                read_applycal_flags()
                applycal_dict_sp[s] = {}
                for j, p in enumerate(previous_cal_targets):
                    applycal_dict_sp[s][p] = [
                        gainfield[j], spwmap[j], interp[j]
                    ]
            else:
                logger.warning('Source {} is not phase-referenced'.format(s))
    else:
        logger.info('Not applying calibration to target sources at this stage')
    eMCP['applycal_dict_sp'] = applycal_dict_sp

    logger.info('End applycal')


### Calibration steps


def initial_bp_cal(eMCP, caltables):
    logger.info('Start initial_bpcal')
    t0 = datetime.datetime.utcnow()
    # Pass 1
    logger.info('Starting pass 1 of initial_bpcal')
    eMCP, caltables = run_bpcal(eMCP, caltables, doplots=False)

    # Apply solutions
    run_applycal(eMCP,
                 caltables,
                 step='bp_apply_mid',
                 insources=eMCP['msinfo']['sources']['bpcal'])

    # Flagging
    if eMCP['defaults']['bandpass']['run_flag']:
        eMCP = flagdata_tfcrop(eMCP, defaults='bandpass')

    # Pass 2
    logger.info('Starting pass 2 of initial_bpcal')
    eMCP, caltables = run_bpcal(eMCP, caltables)

    # Apply calibration if requested:
    if eMCP['input_steps']['bandpass'] == 2:
        run_applycal(eMCP, caltables, step='bandpass')

    flag_statistics(eMCP, step='bandpass')
    save_obj(caltables, os.path.join(caltables['calib_dir'], 'caltables.pkl'))
    bp = eMCP['defaults']['bandpass']
    msg = 'field={0}, combine={1}, solint={2}'.format(
        eMCP['msinfo']['sources']['bpcal'], bp['bp_combine'], bp['bp_solint'])
    eMCP = add_step_time('bandpass', eMCP, msg, t0)
    logger.info('End bandpass')
    return eMCP, caltables


def run_bpcal(eMCP, caltables, doplots=True):
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    bp = eMCP['defaults']['bandpass']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # 0 Delay calibration of bpcal
    caltable_name = bp['delay_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = os.path.join(
        caltables['calib_dir'], caltables['inbase'] + '_' + caltable_name)
    caltables[caltable_name]['previous_cal'] = bp['delay_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = bp['delay_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['delay_spw'])
    caltables[caltable_name]['combine'] = bp['delay_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     bp['delay_combine'])
    caltables[caltable_name]['interp'] = bp['delay_interp']
    caltables[caltable_name]['minblperant'] = bp['delay_minblperant']
    caltables[caltable_name]['minsnr'] = bp['delay_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])

# Plots
    if doplots:
        caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_1.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot,
                            gaintype='K')

        logger.info('Delay plot in: {0}'.format(caltableplot))

    # 1 Phase calibration
    caltable_name = bp['phase_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = bp['phase_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = bp['phase_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['phase_spw'])
    caltables[caltable_name]['combine'] = bp['phase_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     bp['phase_combine'])
    caltables[caltable_name]['interp'] = bp['phase_interp']
    caltables[caltable_name]['minblperant'] = bp['phase_minblperant']
    caltables[caltable_name]['minsnr'] = bp['phase_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])

# Plots
    if doplots:
        caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot,
                            gaintype='G',
                            calmode='p')

        logger.info('BP0_p phase plot: {0}'.format(caltableplot))
    logger.info(
        'Apply phase calibration flags to bandpass, applymode=flagonly')

    applycal(vis=msfile,
             gaintable=caltables[caltable_name]['table'],
             field=msinfo['sources']['bpcal'],
             applymode='flagonly',
             flagbackup=False)
    find_casa_problems()
    read_applycal_flags()

    # 2 Amplitude calibration
    caltable_name = bp['ap_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = bp['ap_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['solint'] = bp['ap_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['ap_spw'])
    caltables[caltable_name]['combine'] = bp['ap_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     bp['ap_combine'])
    caltables[caltable_name]['interp'] = bp['ap_interp']
    caltables[caltable_name]['minblperant'] = bp['ap_minblperant']
    caltables[caltable_name]['minsnr'] = bp['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


# Plots
    if doplots:
        caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot_phs,
                            gaintype='G',
                            calmode='p')

        logger.info('BP0_ap phase plot: {0}'.format(caltableplot_phs))
        caltableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_amp.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot_amp,
                            gaintype='G',
                            calmode='ap')

        logger.info('BP0_ap amp plot: {0}'.format(caltableplot_amp))

    # 3 Bandpass calibration
    caltable_name = bp['bp_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = bp['bp_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['solint'] = bp['bp_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['bp_spw'])
    caltables[caltable_name]['combine'] = bp['bp_combine']
    caltables[caltable_name]['uvrange'] = bp['bp_uvrange']
    caltables[caltable_name]['fillgaps'] = bp['bp_fillgaps']
    caltables[caltable_name]['solnorm'] = bp['bp_solnorm']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     bp['bp_combine'])
    caltables[caltable_name]['interp'] = bp['bp_interp']
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name)
    caltables[caltable_name]['gainfield'] = get_unique_field(
        caltables[caltable_name]['table'])
    #    logger.info('Bandpass0 {0}: {1}'.format(caltable_name,bptable))
    # Plots
    if doplots:
        bptableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        bptableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_amp.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            bptableplot_phs,
                            gaintype='B',
                            calmode='p')

        logger.info('BP0 phase plot: {0}'.format(bptableplot_phs))
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            bptableplot_amp,
                            gaintype='B',
                            calmode='ap')

        logger.info('BP0 amplitude plot: {0}'.format(bptableplot_amp))
    return eMCP, caltables


def initial_gaincal(eMCP, caltables):
    logger.info('Start initial_gaincal')
    t0 = datetime.datetime.utcnow()
    # Pass 1
    logger.info('Starting pass 1 of initial_gaincal')
    # Delay calibration #
    if not eMCP['defaults']['initial_gaincal']['delay']['use_fringefit']:
        eMCP, caltables = solve_delays(eMCP, caltables, doplots=False)
    else:
        logger.info('Full fringe fit selected.')
        eMCP, caltables = delay_fringefit(eMCP, caltables)

    # Initial gain pass 1
    eMCP, caltables = gain_p_ap(eMCP, caltables, doplots=False)
    run_applycal(eMCP, caltables, step='initial_gaincal')

    # Flagging steps:
    flagmode = ''
    if eMCP['defaults']['initial_gaincal']['flagmode'] == 'tfcrop':
        eMCP = flagdata_tfcrop(eMCP, 'initial_gaincal')
        flagmode += 'tfcrop'
    elif eMCP['defaults']['initial_gaincal']['flagmode'] == 'rflag':
        eMCP = flagdata_rflag(eMCP, 'initial_gaincal')
        flagmode += 'rflag'
    else:
        logger.info('No flagging selected')
    ### Pass2
    logger.info('Starting pass 2 of initial_gaincal')
    # Delay calibration #
    if not eMCP['defaults']['initial_gaincal']['delay']['use_fringefit']:
        eMCP, caltables = solve_delays(eMCP, caltables)
    else:
        logger.info('Full fringe fit selected.')
        eMCP, caltables = delay_fringefit(eMCP, caltables)

    # Initial gain pass 2
    eMCP, caltables = gain_p_ap(eMCP, caltables)

    # Apply calibration if requested:
    if eMCP['input_steps']['initial_gaincal'] == 2:
        run_applycal(eMCP, caltables, step='initial_gaincal')
    flag_statistics(eMCP, step='initial_gaincal')
    save_obj(caltables, caltables['calib_dir'] + 'caltables.pkl')
    ini_gaincal = eMCP['defaults']['initial_gaincal']
    msg = 'delay solint={0}, combine={1}, flagmode={2}, p_solint={3}, ' \
          'ap_solint={4}'.format(ini_gaincal['delay']['solint'],
                                 ini_gaincal['delay']['combine'],
                                 flagmode,
                                 ini_gaincal['p_solint'],
                                 ini_gaincal['ap_solint'])
    eMCP = add_step_time('initial_gaincal', eMCP, msg, t0)
    logger.info('End initial_gaincal')
    return eMCP, caltables


def select_calibrator(param, eMCP, add_main=True):
    # Default is using all the calsources. User can specify sources in params
    msinfo = eMCP['msinfo']
    if param == 'default':
        calsources = msinfo['sources']['calsources']
        field = ','.join(np.unique(calsources.split(',')))
    elif param == 'phscals':
        calsources = msinfo['sources']['phscals']
        field = ','.join(np.unique(calsources.split(',')))
    elif param != '':
        field = param
    else:
        logger.critical('No valid delay calibrator. Exiting pipeline.')
        exit_pipeline(eMCP)
    # Check if all selected sources are in mssources
    close_pipeline = False
    for f in field.split(','):
        if f not in msinfo['sources']['mssources'].split(','):
            logger.critical('Selected calibrator {} not in MS!'.format(f))
            close_pipeline = True
    if close_pipeline:
        exit_pipeline(eMCP)
    # Include maincals because needed in ap_G1 fluxscale
    if add_main:
        field = ','.join(
            np.unique(
                np.concatenate([
                    field.split(','), msinfo['sources']['maincal'].split(',')
                ])))
    return field


def solve_delays(eMCP, caltables, doplots=True):
    logger.info('Starting solve_delays')
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    delay = eMCP['defaults']['initial_gaincal']['delay']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    caltable_name = delay['tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = delay['prev_cal']
    caltables[caltable_name]['field'] = select_calibrator(
        delay['delay_cal'], eMCP)
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = delay['solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, delay['spw'])
    caltables[caltable_name]['combine'] = delay['combine']
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     delay['combine'])
    caltables[caltable_name]['interp'] = delay['interp']
    caltables[caltable_name]['minblperant'] = delay['minblperant']
    caltables[caltable_name]['minsnr'] = delay['minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


#    logger.info('Delay calibration {0}: {1}'.format(caltable_name, caltable))
# Plots
    if doplots:
        # 1 No range
        caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_1.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot,
                            gaintype='K')

        logger.info('Delay plot: {0}'.format(caltableplot))
        logger.info('Finished solve_delays')
    return eMCP, caltables


def run_fringefit(msfile,
                  caltables,
                  caltable_name,
                  minblperant=3,
                  minsnr=2,
                  smodel=[]):
    previous_cal = caltables[caltable_name]['previous_cal']
    emutils.rmdir(caltables[caltable_name]['table'])
    logger.info('Running fringefit to generate: {0}'.format(
        caltables[caltable_name]['name']))
    logger.info('Field(s) = {0}, zerorates = {1}'.format(
        caltables[caltable_name]['field'],
        caltables[caltable_name]['zerorates']))
    logger.info('solint = {0}, spw = {1},  combine = {2}'.format(
        caltables[caltable_name]['solint'], caltables[caltable_name]['spw'],
        caltables[caltable_name]['combine']))
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp = [caltables[p]['interp'] for p in previous_cal]
    spwmap = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [
        caltables[p]['gainfield']
        if len(np.atleast_1d(caltables[p]['gainfield'].split(','))) < 2 else
        'nearest' for p in previous_cal
    ]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(
        caltables[caltable_name]['table']))
    # Run CASA task fringefit

    fringefit(vis=msfile,
              caltable=caltables[caltable_name]['table'],
              field=caltables[caltable_name]['field'],
              solint=caltables[caltable_name]['solint'],
              combine=caltables[caltable_name]['combine'],
              spw=caltables[caltable_name]['spw'],
              refant=caltables['refant'],
              antenna='*&*',
              gaintable=gaintable,
              gainfield=gainfield,
              interp=interp,
              spwmap=spwmap,
              zerorates=caltables[caltable_name]['zerorates'])
    find_casa_problems()
    logger.info('caltable {0} in {1}'.format(
        caltables[caltable_name]['name'], caltables[caltable_name]['table']))


def delay_fringefit(eMCP, caltables, doplots=True):
    # Currently does not combine spws because that is not correctly implemented
    # in the pre-release version of task fringefit
    logger.info('Starting delay_fringefit')
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    delay = eMCP['defaults']['initial_gaincal']['delay']
    calsources = eMCP['msinfo']['sources']['calsources']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    caltable_name = delay['tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = delay['prev_cal']
    caltables[caltable_name]['field'] = ','.join(
        np.unique(calsources.split(',')))
    caltables[caltable_name]['solint'] = delay['solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, delay['spw'])
    caltables[caltable_name]['combine'] = delay['combine']
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     delay['combine'])
    caltables[caltable_name]['interp'] = delay['interp']
    caltables[caltable_name]['zerorates'] = delay['zerorates']
    caltables[caltable_name]['minblperant'] = delay['minblperant']
    caltables[caltable_name]['minsnr'] = delay['minsnr']
    caltable = caltables[caltable_name]['table']

    # Calibration
    run_fringefit(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


# Plots
# 1 Phases
    if doplots:
        caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot,
                            gaintype='G',
                            calmode='p')
        logger.info(
            'Fringe calibration (phases) plot in: {0}'.format(caltableplot))
    # 2 Delays
    caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_dela.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot,
                        gaintype='K')
    logger.info(
        'Fringe calibration (delays) plot in: {0}'.format(caltableplot))
    caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_dela2.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot,
                        gaintype='K')
    logger.info(
        'Fringe calibration (delays range) plot in: {0}'.format(caltableplot))
    logger.info('Delay calibration plot: {0}'.format(caltableplot))
    logger.info('End delay_fringefit')
    return eMCP, caltables


def make_spwmap(caltables, combine):
    if 'spw' in combine:
        spwmap = [0] * caltables['num_spw']
    else:
        spwmap = []
    return spwmap


def make_spw(msinfo, spw_list):
    spws = spw_list[0]
    chans = spw_list[1]
    if chans == 'innerchan':
        chans = ':' + msinfo['innerchan']
    elif chans == '':
        chans = ''
    else:
        chans = ':' + chans
    return '{0}{1}'.format(spws, chans)


def gain_p_ap(eMCP, caltables, doplots=True):
    logger.info('Running gain_p_ap')
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    gain_p_ap = eMCP['defaults']['initial_gaincal']
    calsources = eMCP['msinfo']['sources']['calsources']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']

    # 1 Phase calibration
    caltable_name = gain_p_ap['p_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_p_ap['p_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = ','.join(
        np.unique(calsources.split(',')))
    caltables[caltable_name]['solint'] = gain_p_ap['p_solint']
    caltables[caltable_name]['combine'] = gain_p_ap['p_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_p_ap['p_spw'])
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     gain_p_ap['p_combine'])
    caltables[caltable_name]['interp'] = gain_p_ap['p_interp']
    caltables[caltable_name]['minblperant'] = gain_p_ap['p_minblperant']
    caltables[caltable_name]['minsnr'] = gain_p_ap['p_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('G0 phase {0}: {1}'.format(caltable_name,caltable))
# Plots
    if doplots:
        caltableplot = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot,
                            gaintype='G',
                            calmode='p')

        logger.info('G0 phase plot: {0}'.format(caltableplot))
    logger.info(
        'Apply phase calibration flags to calibrators, applymode=flagonly')

    applycal(vis=msfile,
             gaintable=caltables[caltable_name]['table'],
             field=calsources,
             applymode='flagonly',
             flagbackup=False)
    find_casa_problems()
    read_applycal_flags()

    # 2 Amplitude calibration
    caltable_name = gain_p_ap['ap_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_p_ap['ap_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = select_calibrator(
        gain_p_ap['ap_calibrator'], eMCP)
    caltables[caltable_name]['solint'] = gain_p_ap['ap_solint']
    caltables[caltable_name]['combine'] = gain_p_ap['ap_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_p_ap['ap_spw'])
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     gain_p_ap['ap_combine'])
    caltables[caltable_name]['interp'] = gain_p_ap['ap_interp']
    caltables[caltable_name]['minblperant'] = gain_p_ap['ap_minblperant']
    caltables[caltable_name]['minsnr'] = gain_p_ap['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


# Plots
    if doplots:
        caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot_phs,
                            gaintype='G',
                            calmode='p')

        logger.info('G1_ap phase plot: {0}'.format(caltableplot_phs))
        caltableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_amp.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot_amp,
                            gaintype='G',
                            calmode='ap')

        logger.info('G1_ap amp plot: {0}'.format(caltableplot_amp))

    logger.info('Finished gain_p_ap')
    return eMCP, caltables


def find_anten_fluxscale(antennas):
    if 'Lo' in antennas:
        anten_for_flux = ['!Lo']
    else:
        anten_for_flux = ['']
    return anten_for_flux


def read_source_model(model, field, msinfo, eMcalflux):
    logger.info('Found tt0 model for source {0}: {1}'.format(field, model[0]))
    logger.info('Found tt1 model for source {0}: {1}'.format(field, model[1]))
    model_path, modelfilename = os.path.split(model[0])
    model_name = os.path.splitext(modelfilename)[0]
    scaled_model_base = '{0}/{1}_{2}.'.format(model_path, msinfo['msfilename'],
                                              model_name)
    scaled_model = [scaled_model_base + ext for ext in ['tt0', 'tt1']]
    emutils.rmdir(scaled_model[0])
    emutils.rmdir(scaled_model[1])
    flux_in_model = imstat(model[0])
    find_casa_problems()
    logger.info(
        'Flux in model tt0 (sum): {0:5.3g}, rms: {1:5.3g}, mean: {2:5.3g}'.
        format(flux_in_model['sum'][0], flux_in_model['rms'][0],
               flux_in_model['mean'][0]))
    logger.info('Scaling model {0} to eMcalflux {1} Jy to create {2}'.format(
        model[0], eMcalflux, scaled_model[0]))
    logger.info('Scaling model {0} to eMcalflux {1} Jy to create {2}'.format(
        model[1], eMcalflux, scaled_model[1]))
    logger.info('The multiplying factor is: {0:5.3f}'.format(
        eMcalflux / flux_in_model['sum'][0]))
    immath(imagename=model[0],
           mode='evalexpr',
           expr='IM0*{}'.format(eMcalflux / flux_in_model['sum'][0]),
           outfile=scaled_model[0])
    immath(imagename=model[1],
           mode='evalexpr',
           expr='IM0*{}'.format(eMcalflux / flux_in_model['sum'][0]),
           outfile=scaled_model[1])
    return scaled_model


def eM_fluxscale(eMCP, caltables):
    logger.info('Start eM_fluxscale')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    flux = eMCP['defaults']['fluxscale']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    sources = eMCP['msinfo']['sources']
    antennas = eMCP['msinfo']['antennas']
    ampcal_table = flux['ampcal_table']
    anten_for_flux = find_anten_fluxscale(antennas)
    cals_to_scale = select_calibrator(
        eMCP['defaults']['initial_gaincal']['ap_calibrator'], eMCP)
    fluxcal = sources['fluxcal']
    caltable_name = flux['tablename']
    # Remove previous results:
    emutils.rmfile(calib_dir + 'allcal_ap.G1_fluxes.txt')
    emutils.rmfile(calib_dir +
                   '{0}_fluxscale.png'.format(msinfo['msfilename']))
    # Check if table allcal_ap.G1 exists:
    check_table_exists(caltables, ampcal_table)
    caltables[caltable_name] = copy.copy(caltables[ampcal_table])
    caltables[caltable_name][
        'table'] = caltables[ampcal_table]['table'] + '_fluxscaled'
    fluxes_txt = calib_dir + caltables[caltable_name]['name'] + '_fluxes.txt'
    logger.info('Flux density scale from: {0}'.format(fluxcal))
    logger.info('Transfered to: {0}'.format(cals_to_scale))
    logger.info('Input caltable: {0}'.format(caltables[ampcal_table]['table']))
    logger.info('Antennas used to scale: {0}'.format(anten_for_flux))
    emutils.rmdir(caltables[caltable_name]['table'])

    run_fluxscale(msfile, fluxcal, cals_to_scale, ','.join(anten_for_flux),
                  caltables[ampcal_table]['table'],
                  caltables[caltable_name]['table'], fluxes_txt)
    calfluxes = load_obj(calib_dir + 'calfluxes.pkl')

    logger.info('Modified caltable: {0}'.format(
        caltables[caltable_name]['table']))
    logger.info(
        'Spectrum information: {0}'.format(caltables[caltable_name]['table'] +
                                           '_fluxes.txt'))
    if calfluxes is None:
        logger.critical('Something went wrong with fluxscale')
        logger.warning('This probably means that necessary data are missing:')
        logger.warning('Required sources: {0}'.format(cals_to_scale))
        logger.warning('Required antennas: {0}'.format(
            ','.join(anten_for_flux)))
        exit_pipeline(eMCP)
    # Compute correction to scale the flux density of 3C286 according to
    # resolution provided by the shortest available baseline of e-MERLIN
    eMfactor = calc_eMfactor(msfile, field=fluxcal)
    # Include a note in the fluxes.txt file warning that the values in that
    # file should be corrected by eMfactor
    with open(fluxes_txt, 'a') as file:
        file.write(
            '# WARNING: All flux densities in this file need to be multiplied by eMfactor={0:6.4f} to match the corrections that have been applied to the data.'
            .format(eMfactor))
    # Get fitted flux density and spectral index, correctly scaled for e-MERLIN
    eMcalfluxes = {}
    for k in calfluxes.keys():
        if type(calfluxes[k]) is dict:
            try:
                a = []
                a.append(calfluxes[k]['fitFluxd'] * eMfactor)
                a.append(calfluxes[k]['spidx'][1])
                a.append(calfluxes[k]['fitRefFreq'])
                eMcalfluxes[calfluxes[k]['fieldName']] = a
                logger.info(
                    'Spectrum for {0:>9s}: Flux density = {1:6.3f} +/-{2:6.3f}, spidx ={3:5.2f}+/-{4:5.2f}'
                    .format(calfluxes[k]['fieldName'],
                            calfluxes[k]['fitFluxd'] * eMfactor,
                            calfluxes[k]['fitFluxdErr'] * eMfactor,
                            calfluxes[k]['spidx'][1],
                            calfluxes[k]['spidxerr'][1]))
            except:
                pass
    logger.info('Plotting fluxscale models')
    emplt.fluxscale_models(calfluxes, eMfactor, msinfo)
    # Scale and fill model column:
    for field in eMcalfluxes.keys():
        # Check if there are image models (tt0, tt1) for each field:
        model = [
            './source_models/{0}.model.tt0'.format(field),
            './source_models/{0}.model.tt1'.format(field)
        ]
        if os.path.exists(model[0]) and os.path.exists(model[1]):
            scaled_model = read_source_model(model, field, msinfo,
                                             eMcalfluxes[field][0])
            logger.info('New model for this observation: {0}, {1}'.format(
                scaled_model[0], scaled_model[1]))
            ft(vis=msfile,
               field=field,
               model=scaled_model,
               nterms=2,
               usescratch=True)
            logger.info('Model for {0} included in MODEL column in {1}'.format(
                field, msfile))
            if eMCP['is_mixed_mode']:
                msfile_sp = eMCP['msinfo']['msfile_sp']
                logger.info('Filling sp: {0}'.format(msfile_sp))
                ft(vis=msfile_sp,
                   field=field,
                   model=scaled_model,
                   nterms=2,
                   usescratch=True)
        else:
            # If there is no model, just assume point source:
            logger.info(
                'Filling model column for point-like source: {0}'.format(
                    field))
            logger.info(
                'Model: flux={0:6.3g}Jy, spix={1:6.3g}, reffreq={2:6.3g}Hz'.
                format(eMcalfluxes[field][0], eMcalfluxes[field][1],
                       eMcalfluxes[field][2]))
            setjy(vis=msfile,
                  field=field,
                  standard='manual',
                  fluxdensity=eMcalfluxes[field][0],
                  spix=eMcalfluxes[field][1],
                  reffreq=str(eMcalfluxes[field][2]) + 'Hz',
                  usescratch=True)
            find_casa_problems()
            if eMCP['is_mixed_mode']:
                msfile_sp = eMCP['msinfo']['msfile_sp']
                logger.info('Filling model for sp: {0}'.format(msfile_sp))

                setjy(vis=msfile_sp,
                      field=field,
                      standard='manual',
                      fluxdensity=eMcalfluxes[field][0],
                      spix=eMcalfluxes[field][1],
                      reffreq=str(eMcalfluxes[field][2]) + 'Hz',
                      usescratch=True)
                find_casa_problems()

    # Apply calibration if requested:
    if eMCP['input_steps']['fluxscale'] == 2:
        run_applycal(eMCP, caltables, step='fluxscale')
    logger.info('End fluxscale')
    save_obj(caltables, caltables['calib_dir'] + 'caltables.pkl')
    msg = ''
    eMCP = add_step_time('fluxscale', eMCP, msg, t0)
    return eMCP, caltables


def compile_delays(tablename, outname):
    tb.open(tablename + '/ANTENNA')
    antennas = tb.getcol('NAME')
    tb.close()
    tb.open(tablename)
    a = tb.getcol('ANTENNA1')
    times = tb.getcol('TIME')
    delays = tb.getcol('FPARAM')
    tb.close()
    delay_stats = []
    for i in range(len(times)):
        delay_stats.append([antennas[a[i]], 'RR', delays[0, 0, i]])
        delay_stats.append([antennas[a[i]], 'LL', delays[1, 0, i]])
    delay_stats = np.asarray(delay_stats)
    np.save('delay_' + outname, delay_stats)
    logger.info('Delay statistics saved to: {0}'.format(outname))


def calc_eMfactor(msfile, field='1331+305'):
    logger.info('Computing eMfactor')
    if field not in ['1331+305', '1331+3030', 'J1331+305', 'J1331+3030']:
        logger.warning(
            'Scaling flux assuming 3C286 is the flux calibrator. Your flux calibrator is: {}. Scaling could wrong.'
            .format(field))
        logger.info('Assuming eMfactor = 1')
        eMfactor = 1.0
        return eMfactor
    tb.open(msfile + '/FIELD')
    names = tb.getcol('NAME')
    field_id = np.argwhere(names == field)[0][0]
    tb.close()

    tb.open(msfile + '/ANTENNA')
    anten = tb.getcol('NAME')
    tb.close()
    try:
        Lo_id = np.argwhere(anten == 'Lo')[0][0]
    except:
        Lo_id = -1

    tb.open(msfile)
    uvw = tb.getcol('UVW')
    a1 = tb.getcol('ANTENNA1')
    a2 = tb.getcol('ANTENNA2')
    field = tb.getcol('FIELD_ID')
    tb.close()

    uvdist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    mask = (uvdist == 0) + (field != field_id) + (anten[a1] == 'Lo')
    uvdist_nonzero = np.ma.array(uvdist, mask=mask)

    # To exclude completely flagged data:
    # I comment this out because for some data sets tb.getcol returns a flatten
    # array without shape

    n = np.argmin(uvdist_nonzero)

    tb.open(msfile + '/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    tb.close()

    shortest_baseline = uvdist_nonzero[n]  # Shortest baseline in m
    center_freq = np.mean(chan_freq) / 1e6  # Center frequency in MHz
    dfluxpy_output = dfluxpy(center_freq, shortest_baseline)
    eMfactor = dfluxpy_output[1] / dfluxpy_output[0]
    logger.info(
        'Shortest projected baseline: {0} [{1}] - {2} [{3}] {4:10.2f}m'.format(
            anten[a1[n]], a1[n], anten[a2[n]], a2[n], uvdist_nonzero[n]))
    logger.info('Central frequency of the MS: {0} MHz'.format(center_freq))
    logger.info('eMfactor: {0:6.4f}'.format(eMfactor))
    return eMfactor


def dfluxpy(freq, baseline):
    import math
    #######
    # Python version of 3C286 flux calculation program (original author unknown)
    # ..............................
    # Author DMF       20/10/2011
    # ..............................
    #
    # Update to use Perley & Butler 2012 coefficients
    # 10/04/2013
    # DMF
    ########

    # Reworked to use the 1999 VLA flux formula, and a 2nd formula to give a continuous
    # estimate of the resolved fraction,
    # by Ian Stewart, JBO, 8 Aug 2007.
    # Minor changes by amsr, 8 Aug 2007

    lowest_freq = 300.0
    highest_freq = 50000.0
    if freq < lowest_freq or freq > highest_freq:
        print(
            "Frequency must be between $lowest_freq and $highest_freq MHz. \n")

    # Perley & Butler 2012 values
    A = 1.2515
    B = -0.4605
    C = -0.1715
    D = 0.0336

    log10f = (math.log(freq) / 2.3025851) - 3.0
    # Why the -3? Because freq has to be GHz for the formula to work.
    log_flux = A + B * log10f + C * log10f**2 + D * log10f**3
    vlaflux = math.pow(10.0, log_flux)

    # The VLA flux must now be corrected to account for the higher resolving power of merlin.
    # The formula used was obtained with the help of Peter Thomasson. If we assume that 3C286 is represented
    # by a gaussian of angular size theta_s, and represent the resolving power of the given baseline as a
    # function of frequency f and antenna separation L by theta_b(f,L), then the reduction in central
    # flux density A(0) due to the finite theta_s is given by
    #
    #                               1
    #                -----------------------------------
    #     A'(0)       2 pi (theta_b(f,L)^2 + theta_s^2)
    #    ------- = --------------------------------------- ,
    #     A(0)                      1
    #                       ---------------------
    #                        2 pi theta_b(f,L)^2
    #
    #                   1
    #            = -------------- ,
    #               1 + rho(f,L)
    #
    # where the resolved fraction rho(f,L) is given by
    #
    #                  theta_s^2
    #    rho(f,L) = ---------------- .
    #                theta_b(f,L)^2
    #
    # Use of theta_b(f,L) = k/(fL) allows this to be written
    #
    #               (   f*L     )^2
    #    rho(f,L) = (-----------)   * rho_ref .
    #               (f_ref*L_ref)
    #
    # The reference value of rho is fixed at 0.04 for the MK-TA baseline at 5 GHz (Peter Thomasson).

    ref_bl_length = 11236.79  # MK-TA separation in metres.
    ref_freq = 5000.0
    ref_rho = 0.04
    thisbl = "this baseline (Mk-Ta)"

    bl_length = baseline

    frac = (freq / ref_freq) * (bl_length / ref_bl_length)
    rho = frac**2 * ref_rho
    merlinflux = vlaflux / (1.0 + rho)

    # Another useful quantity is the resolved percentage:
    resolved_percent = 100.0 * rho / (1.0 + rho)
    caution_res_pc = 10.0

    return vlaflux, merlinflux, resolved_percent, caution_res_pc, thisbl


def plot_image(eMCP, imagename, center, ext='.tt0', dozoom=False):
    imstat_residual = imstat(imagename + '.residual' + ext)
    find_casa_problems()
    imstat_image = imstat(imagename + '.image' + ext)
    find_casa_problems()
    noise = imstat_residual['rms'][0]
    peak = imstat_image['max'][0]
    scaling = np.min([0, -np.log(1.0 * peak / noise) + 4])
    imgpar = eMCP['defaults']['first_images']
    level0 = imgpar['level0']
    levels = level0 * np.sqrt(3**np.arange(20))
    zoom_range = 150
    for extension in ['image' + ext]:
        filename = '{0}.{1}'.format(imagename, extension)
        logger.info('Creating png for image: {0}'.format(filename))
        logger.info('Peak: {0:5.2e} mJy, noise: {1:5.2e} mJy, ' \
                    'scaling: {2:3.1f}'.format(peak*1000.,
                                               noise*1000.,
                                               scaling))
        imview(raster={
            'file': filename,
            'scaling': float(scaling),
            'colorwedge': True
        },
               contour={
                   'file': filename,
                   'levels': list(levels),
                   'base': 0,
                   'unit': float(noise) * 2.
               },
               out=filename + '.png')
        find_casa_problems()
        if dozoom:
            imview(raster={
                'file': filename,
                'scaling': float(scaling),
                'colorwedge': True
            },
                   contour={
                       'file': filename,
                       'levels': list(levels),
                       'base': 0,
                       'unit': float(noise) * 2.
                   },
                   zoom={
                       'blc': [center - zoom_range, center - zoom_range],
                       'trc': [center + zoom_range, center + zoom_range]
                   },
                   out=filename + '_zoom.png')
            find_casa_problems()
    return peak, noise, scaling


def plot_image_add(imagename, center, zoom_range, ext='.tt0', dozoom=False):
    filename = imagename
    imview(raster={
        'file': filename + '.residual' + ext,
        'colorwedge': True
    },
           contour={
               'file': filename + '.mask',
               'levels': [1]
           },
           out=filename + '.residual' + ext + '.png')
    find_casa_problems()
    if dozoom:
        imview(raster={
            'file': filename + '.residual' + ext,
            'colorwedge': True
        },
               contour={
                   'file': filename + '.mask',
                   'levels': [1]
               },
               zoom={
                   'blc': [center - zoom_range, center - zoom_range],
                   'trc': [center + zoom_range, center + zoom_range]
               },
               out=filename + '.residual' + ext + '_zoom.png')
        find_casa_problems()
    return


def write_wsclean_command(msfile, config_wsclean):
    logger.debug('config_wsclean')
    # Duplicate size if needed:
    if type(config_wsclean['-size']) == int:
        size_int = config_wsclean['-size']
        config_wsclean['-size'] = '{0} {0}'.format(size_int)


#    # Skip casa-mask if no mask specified
#    if config_wsclean['-casa-mask'] == '':
#        config_wsclean.pop('-casa-mask')
# Only keep keys starting with - that will be passed to wsclean
    for key in config_wsclean.keys():
        if key[0] != '-':
            config_wsclean = config_wsclean.pop(key)

    logger.debug(config_wsclean)
    wsclean_params = ' '.join(
        ['{0} {1}'.format(k, v) for (k, v) in config_wsclean.items()])
    wsclean_command = '{0} {1} {2}'.format('wsclean', wsclean_params, msfile)
    return wsclean_command


def single_tclean(eMCP, s, num):
    msinfo = eMCP['msinfo']
    logger.info('Producing tclean images for {0}, field: {1}'.format(
        msinfo['msfile'], s))
    emutils.makedir(images_dir + '{}'.format(s))
    imagename = images_dir + '{0}/{1}_{0}_img{2:02d}'.format(
        s, msinfo['msfilename'], num)
    prev_images = glob.glob(imagename + '*')
    for prev_image in prev_images:
        emutils.rmdir(prev_image)
    cellsize = {'C': '0.008arcsec', 'L': '0.02arcsec', 'K': '0.002arcsec'}
    imgpar = eMCP['defaults']['first_images']
    imsize = imgpar['imsize']
    cell = cellsize[msinfo['band']]
    niter = imgpar['niter']
    gain = imgpar['gain']
    uvtaper = imgpar['uvtaper']
    uvrange = imgpar['uvrange']
    restoringbeam = imgpar['restoringbeam']
    deconvolver = imgpar['deconvolver']
    nterms = imgpar['nterms']
    scales = imgpar['scales']
    weighting = imgpar['weighting']
    robust = imgpar['robust']
    usemask = 'auto-multithresh'
    nsigma = imgpar['nsigma']
    sidelobethreshold = imgpar['sidelobethreshold']
    noisethreshold = imgpar['noisethreshold']
    lownoisethreshold = imgpar['lownoisethreshold']
    minbeamfrac = imgpar['minbeamfrac']
    growiterations = imgpar['growiterations']
    parallel = False

    logger.info('imsize = {0}, cell = {1}, niter = {2}'.format(
        imsize, cell, niter))
    logger.info('weighting = {0}, robust = {1}'.format(weighting, robust))
    logger.info('usemask = {0}, growiterations = {1}'.format(
        usemask, growiterations))
    tclean(vis=msinfo['msfile'],
           field=s,
           datacolumn='corrected',
           imagename=imagename,
           imsize=imsize,
           cell=cell,
           deconvolver=deconvolver,
           gain=gain,
           nterms=nterms,
           scales=scales,
           uvtaper=uvtaper,
           uvrange=uvrange,
           restoringbeam=restoringbeam,
           weighting=weighting,
           robust=robust,
           niter=niter,
           usemask=usemask,
           nsigma=nsigma,
           sidelobethreshold=sidelobethreshold,
           noisethreshold=noisethreshold,
           lownoisethreshold=lownoisethreshold,
           minbeamfrac=minbeamfrac,
           growiterations=growiterations,
           savemodel='none',
           parallel=parallel)
    find_casa_problems()
    if nterms > 1:
        ext = '.tt0'
    else:
        ext = ''
    peak, noise, scaling = plot_image(eMCP,
                                      imagename,
                                      center=int(imsize / 2.0),
                                      ext=ext,
                                      dozoom=True)
    zoom_range = eMCP['defaults']['first_images']['zoom_range_pix']
    plot_image_add(imagename,
                   center=int(imsize / 2.0),
                   zoom_range=zoom_range,
                   ext=ext,
                   dozoom=True)
    eMCP['img_stats'][s] = [peak, noise, scaling]
    return eMCP


def single_wsclean(eMCP, s, field_id):
    # TODO: Here we should call either wsclean or tclean
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    # Check if wsclean is available:
    wsclean_available = check_command('wsclean')
    if not wsclean_available:
        logger.critical('wsclean not available.')
        logger.warning('Exiting pipeline.')
        exit_pipeline(eMCP)
    logger.info(f"Producing wsclean images for {msfile}, field: {s}")
    emutils.makedir(os.path.join(images_dir, '{}'.format(s)))
    num = 0
    imagename = os.path.join(images_dir, s,
                             f"{msinfo['msfilename']}_{s}_img{num:02d}")
    logger.debug(f"Imagename: {imagename}")
    prev_images = glob.glob(imagename + '*')
    for prev_image in prev_images:
        emutils.rmdir(prev_image)
    config_wsclean = eMCP['defaults']['first_images']['wsclean']
    config_wsclean['-name'] = imagename
    config_wsclean['-field'] = field_id
    # Set up cellsize:
    cellsize = {'C': '0.008asec', 'L': '0.02asec', 'K': '0.002asec'}
    if config_wsclean['-scale'] == 'auto':
        config_wsclean['-scale'] = cellsize[msinfo['band']]
    logger.info(
        f"imsize = {config_wsclean['-size']}, scale = {config_wsclean['-scale']}, niter = {config_wsclean['-niter']}"
    )
    logger.info(f"weight = {config_wsclean['-weight']}")
    wsclean_command = write_wsclean_command(msfile, config_wsclean)
    logger.info(f'Full wsclean command:\n{wsclean_command}')
    logger.info(f'{wsclean_command.split()}')

    with open('stdouterr.log', 'a') as f:
        subprocess.run(wsclean_command.split(),
                       stdout=f,
                       stderr=subprocess.STDOUT)


#    if nterms > 1:
#        ext = '.tt0'
#    else:
#        ext = ''
    fitsfile = imagename + '-image.fits'
    # Image statistics
    imstats_img = get_image_stats(fitsfile)
    imstats_res = get_image_stats(fitsfile.replace('-image', '-residual'))
    scaling = np.min(
        [0, -np.log(1.0 * imstats_img['max'] / imstats_res['rms']) + 4])
    eMCP['img_stats'][s] = [imstats_img['max'], imstats_res['rms'], scaling]
    logger.debug(imstats_img['max'])
    logger.debug(imstats_res['rms'])
    logger.debug(scaling)
    # Convert to png
    emplt.fits2png(fitsfile,
                   rms=imstats_res['rms'],
                   scaling=scaling,
                   contour=False)
    emplt.fits2png(fitsfile,
                   rms=imstats_res['rms'],
                   scaling=scaling,
                   zoom=True)
    emplt.fits2png(fitsfile.replace('-image', '-residual'),
                   scaling=scaling,
                   rms=imstats_res['rms'],
                   contour=False)
    emplt.fits2png(fitsfile.replace('-image', '-residual'),
                   scaling=scaling,
                   rms=imstats_res['rms'],
                   contour=False,
                   zoom=True)
    return eMCP


def get_image_stats(image_file):
    logger.info(f"Getting statistics from {image_file}")
    with fits.open(image_file) as hdu_list:
        image_data = hdu_list[0].data
    imstats = {}
    imstats['min'] = np.min(image_data)
    imstats['max'] = np.max(image_data)
    imstats['mean'] = np.mean(image_data)
    imstats['stdev'] = np.std(image_data)
    imstats['rms'] = np.sqrt(np.mean(np.square(image_data)))
    logger.info(f"Min: {imstats['min']}, Max: {imstats['max']}")
    return imstats


def run_first_images(eMCP):
    msinfo = eMCP['msinfo']
    logger.info(line0)
    logger.info('Start run_first_images')
    t0 = datetime.datetime.utcnow()

    eMCP['img_stats'] = {}
    field_names = np.array(
        emutils.read_keyword(msinfo['msfile'], 'NAME', subtable='FIELD'))
    for s in msinfo['sources']['targets_phscals'].split(','):
        field_id = np.argwhere(field_names == s)[0][0]
        #        eMCP = single_tclean(eMCP, s, num)
        eMCP = single_wsclean(eMCP, s, field_id)
    logger.info('End first_images')
    msg = ''
    eMCP = add_step_time('first_images', eMCP, msg, t0)
    return eMCP


def run_split_fields(eMCP):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    logger.info(line0)
    logger.info('Start split_fields')
    t0 = datetime.datetime.utcnow()
    split_fields_defaults = eMCP['defaults']['split_fields']
    fields_to_split = split_fields_defaults['fields']
    datacolumn = split_fields_defaults['datacolumn']
    createmms = split_fields_defaults['createmms']
    output_dir = split_fields_defaults['output_dir']
    emutils.makedir(output_dir)
    if fields_to_split in msinfo['sources'].keys():
        fields_str = msinfo['sources'][fields_to_split]
    elif fields_to_split == '':
        logger.critical('No fields selected to split')
        exit_pipeline(eMCP)
    else:
        fields_str = fields_to_split
    fields = sorted(fields_str.split(','))
    logger.info('Selected fields: {}'.format(fields))
    for i, field in enumerate(fields):
        logger.info('Processing field: {}'.format(field))
        chanbin = split_fields_defaults['chanbin']
        timebin = split_fields_defaults['timebin']
        chanaverage = split_fields_defaults['chanaverage']
        timeaverage = split_fields_defaults['timeaverage']
        # If a list was provided, use corresponding value
        if type(timebin) is list:
            logger.info('timebin is a list: {}'.format(timebin))
            timebin = timebin[i]
            timeaverage = timeaverage[i]
        if type(chanbin) is list:
            logger.info('chanbin is a list: {}'.format(chanbin))
            chanbin = chanbin[i]
            chanaverage = chanaverage[i]
        # Don't use channel average if it only 1 channel selected
        if chanaverage is True and chanbin == 1:
            chanaverage = False
        # Define output name
        msfile_name = '{0}_{1}.ms'.format(msinfo['run'], field)
        outputmsfile = os.path.join(output_dir, msfile_name)
        emutils.rmdir(outputmsfile)
        emutils.rmdir(outputmsfile + '.flagversions')
        logger.info('Input MS: {0}'.format(msfile))
        logger.info('Output MS: {0}'.format(outputmsfile))
        if chanaverage:
            logger.info('chanbin={0}'.format(chanbin))
        if timeaverage:
            logger.info('timebin={0}'.format(timebin))
        logger.info('Datacolumn: {0}'.format(datacolumn))
        if createmms:
            logger.info('Create MMS: {0}'.format(createmms))
        mstransform(vis=msfile,
                    outputvis=outputmsfile,
                    field=field,
                    timeaverage=timeaverage,
                    chanaverage=chanaverage,
                    timebin=str(timebin),
                    chanbin=chanbin,
                    datacolumn=datacolumn,
                    keepflags=True)
        find_casa_problems()
        flagtable_info = 'after_split'
        current_time = datetime.datetime.utcnow()
        flagmanager(vis=outputmsfile,
                    mode='save',
                    versionname=flagtable_info,
                    comment='After splitting the source at {0}'.format(
                        str(current_time)))
        find_casa_problems()
        if not os.path.exists(outputmsfile):
            logger.critical('Could not split field')
            exit_pipeline(eMCP)
    logger.info('End split_fields')
    msg = ''
    eMCP = add_step_time('split_fields', eMCP, msg, t0)
    return eMCP


def shift_field_position(eMCP, msfile, shift):
    field = shift['field']
    new_pos = shift['new_position']
    position_name = shift['new_field_name']
    logger.info('Field {0} will be shifted to {1} on {2}'.format(
        field, position_name, new_pos))
    msfile_split = '{0}_{1}'.format(msfile, position_name)
    mssources = vishead(msfile, mode='list', listitems='field')['field'][0]
    if field not in mssources:
        logger.critical(
            'Requested field to shift: {} not in MS! Closing '.format(field))
        exit_pipeline(eMCP)
    emutils.rmdir(msfile_split)
    # Split
    logger.info('Splitting field: {}'.format(field))
    mstransform(msfile, outputvis=msfile_split, field=field, datacolumn='data')
    find_casa_problems()
    logger.info('Changing phase center to: {}'.format(new_pos))
    phaseshift(vis=msfile_split,
               field=field,
               outputvis=msfile_split,
               phasecenter=new_pos,
               datacolumn='data')
    find_casa_problems()
    # Change field name
    tb.open(msfile_split + '/FIELD', nomodify=False)
    st = tb.selectrows(0)
    st.putcol('NAME', '{0}'.format(position_name))
    st.done()
    tb.close()
    # Average individual field
    chanbin = eMCP['defaults']['average']['chanbin']
    timebin = eMCP['defaults']['average']['timebin']
    if timebin == '1s':
        timeaverage = False
    else:
        timeaverage = True
    if chanbin == 1:
        chanaverage = False
    else:
        chanaverage = True
    datacolumn = eMCP['defaults']['average']['datacolumn']
    scan = eMCP['defaults']['average']['scan']
    antenna = eMCP['defaults']['average']['antenna']
    timerange = eMCP['defaults']['average']['timerange']
    emutils.rmdir(msfile_split + '_avg')
    mstransform(vis=msfile_split,
                outputvis=msfile_split + '_avg',
                timeaverage=timeaverage,
                chanaverage=chanaverage,
                timerange=timerange,
                scan=scan,
                antenna=antenna,
                timebin=timebin,
                chanbin=chanbin,
                datacolumn=datacolumn,
                keepflags=True)
    find_casa_problems()
    if os.path.isdir(msfile_split + '_avg'):
        emutils.rmdir(msfile_split)
    # Concatenate again
    s = msfile.split('.')
    sn, se = s[:-1], s[-1]
    msfile0 = '.'.join(sn) + '_avg.' + se
    logger.info('Concatenating {0} into {1}'.format(msfile_split + '_avg',
                                                    msfile0))
    concat(vis=msfile_split + '_avg', concatvis=msfile0)
    emutils.rmdir(msfile_split + '_avg')
    logger.info('Updating listobs for MS: {}'.format(msfile0))
    run_listobs(msfile0)
    logger.warning(
        'New field: {} in MS. Make sure you include it in the inputs file'.
        format(position_name))


def read_shifts_file(shifts_file):
    shifts_list = []
    with open(shifts_file, 'r') as shifts_txt:
        lines = shifts_txt.readlines()
        for line in lines:
            if line.strip() != '':
                shift = {
                    'field': line.split(',')[0].strip(),
                    'new_field_name': line.split(',')[1].strip(),
                    'new_position': line.split(',')[2].strip()
                }
                shifts_list.append(shift)
    return shifts_list


def shift_all_positions(eMCP):
    logger.info(line0)
    msfile = eMCP['msinfo']['msfile']
    logger.info('Running shift_all_pos')
    #    t0 = datetime.datetime.utcnow()
    shifts_file = './shift_phasecenter.txt'
    try:
        shifts_list = read_shifts_file(shifts_file)
        logger.info('Reading shifts from {0}'.format(shifts_file))
    except:
        logger.critical(
            'Unable to open {0} with new position information'.format(
                shifts_file))
        exit_pipeline(eMCP)
    listobs_file = info_dir + msfile
    emutils.mvdir(listobs_file + '.listobs.txt',
                  listobs_file + 'preshift_listobs.txt')
    logger.info(
        'Found {0} shifts to apply. {0} new fields will be added'.format(
            len(shifts_list)))
    for shift in shifts_list:
        shift_field_position(eMCP, msfile, shift)
    logger.info('Finished shift_all_pos')

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


def find_fields_scans(msfile):
    msmd.open(msfile)
    scans = msmd.scannumbers()
    dict_scans = msmd.fieldsforscans(scans,
                                     True,
                                     asmap=True,
                                     obsid=0,
                                     arrayid=0)
    msmd.done()
    field_for_scan = np.array([dict_scans[str(scan)][0] for scan in scans])
    return scans, field_for_scan


def flag_Lo_dropouts(eMCP):
    # Search for Lo dropouts
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    antennas = get_antennas(msfile)
    msg_out = ''
    if 'Lo' in antennas:
        Lo_dropout_scans = eMCP['defaults']['flag_manual_avg']['Lo_dropout']
        if Lo_dropout_scans == 'none':
            eMCP['msinfo']['Lo_dropout_scans'] = 'none'
        elif Lo_dropout_scans == '':
            if msinfo['sources']['phscals'] != '':
                logger.info('Searching for Lo dropout scans')
                logger.info(
                    "To avoid this, set deafult Lo_dropout_scans = 'none' "
                    "when rerun flag_manual_avg")
                phscals = msinfo['sources']['phscals'].split(',')
                Lo_drop_list = find_Lo_drops(msfile, phscals, eMCP)
                if not Lo_drop_list:
                    eMCP['msinfo']['Lo_dropout_scans'] = 'none'
                else:
                    eMCP['msinfo']['Lo_dropout_scans'] = ','.join(
                        Lo_drop_list.astype('str'))
                    logger.info('Flagging Lo dropout scans: '
                                '{0}'.format(
                                    eMCP['msinfo']['Lo_dropout_scans']))
                    flagdata(vis=msfile,
                             antenna='Lo',
                             scan=eMCP['msinfo']['Lo_dropout_scans'],
                             flagbackup=False)
                    find_casa_problems()
                msg_out = ' Searched Lo dropouts'
        else:
            logger.info('Lo_dropout_scans manually selected: {0}'.format(
                Lo_dropout_scans))
            eMCP['msinfo']['Lo_dropout_scans'] = Lo_dropout_scans
            logger.info('Now flagging')
            flagdata(vis=msfile,
                     antenna='Lo',
                     scan=eMCP['msinfo']['Lo_dropout_scans'],
                     flagbackup=False)
            find_casa_problems()
    else:
        eMCP['msinfo']['Lo_dropout_scans'] = 'none'
    return eMCP, msg_out


def find_Lo_amp_spw(msfile, phscal, phscal_scans, spw, eMCP):
    amp_mean = np.ones_like(phscal_scans) * np.nan
    amp_std = np.ones_like(phscal_scans) * np.nan
    Lo_defaults = eMCP['defaults']['flag_manual_avg']
    results_tmp = visstat(msfile,
                          scan='',
                          field=phscal,
                          antenna='Lo&*',
                          spw=str(spw + ':' + eMCP['msinfo']['innerchan']),
                          datacolumn=Lo_defaults['Lo_datacolumn'],
                          useflags=Lo_defaults['Lo_useflags'],
                          correlation='RR,LL',
                          timeaverage=True,
                          timebin='9999999s',
                          timespan='')
    for key in results_tmp.keys():
        scan = key.split(',')[1].split('=')[-1]
        scan_idx = np.where(phscal_scans == int(scan))[0][0]
        amp_mean[scan_idx] = results_tmp[key]['median']
        amp_std[scan_idx] = results_tmp[key]['stddev']
    return amp_mean, amp_std


def find_Lo_amp(msfile, phscal, phscal_scans, eMCP, spws):
    amp_means = np.zeros((len(spws), len(phscal_scans)))
    amp_stds = np.zeros((len(spws), len(phscal_scans)))
    logger.info('Analysing {0} scans for phscal: {1}'.format(
        len(phscal_scans), phscal))
    for i, spw in enumerate(spws):
        logger.info('Processing spw: {}'.format(spw))
        amp_means_i, amp_std_i = find_Lo_amp_spw(msfile, phscal, phscal_scans,
                                                 spw, eMCP)
        amp_means[i] = amp_means_i
        amp_stds[i] = amp_std_i
    amp_mean = np.average(amp_means, weights=1.0 / amp_stds, axis=0)
    amp_std = np.average(amp_stds, axis=0)
    return amp_mean, amp_std


def plot_Lo_drops(msfile, phscal_scans, amp_mean, lo_dropout_scans, phscal,
                  eMCP):
    msinfo = eMCP['msinfo']
    drops = np.array([scan in lo_dropout_scans for scan in phscal_scans])
    fig = plt.figure(figsize=(30, 8))
    ax1 = fig.add_subplot(111)

    scans, field_for_scan = find_fields_scans(msfile)
    ax1.bar(scans - 0.5,
            np.ones_like(scans) * np.max(amp_mean) * 1.2,
            alpha=0.2,
            color='0.5',
            width=1),

    ax1.bar(phscal_scans - 0.5,
            amp_mean,
            alpha=1.0,
            color='0.5',
            width=1,
            label='{0}'.format(phscal))
    ax1.bar(phscal_scans[drops] - 0.5,
            amp_mean[drops],
            alpha=1.0,
            color='r',
            width=1,
            label='{0} Lo dropouts'.format(phscal))

    ax1.legend(loc=0)
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xlim(np.min(phscal_scans) - 0.5, np.max(phscal_scans) + 0.5)
    ax1.set_ylim(0, np.max(amp_mean) * 1.2)
    ax1.set_xlabel('Scan number')
    ax1.set_ylabel('Mean spw Lo raw amplitude')

    plots_obs_dir = './weblog/plots/plots_flagstats/'
    plot_file_Lo = plots_obs_dir + '{0}_Lo_dropout_scans{1}.png'.format(
        msinfo['msfilename'], phscal)
    fig.savefig(plot_file_Lo, bbox_inches='tight')


def bandpass_final(eMCP, caltables):
    logger.info('Start bandpass_final')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    bp_final = eMCP['defaults']['bandpass_final']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # Bandpass calibration
    caltable_name = bp_final['bp_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = bp_final['bp_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['solint'] = bp_final['bp_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp_final['bp_spw'])
    caltables[caltable_name]['combine'] = bp_final['bp_combine']
    caltables[caltable_name]['uvrange'] = bp_final['bp_uvrange']
    caltables[caltable_name]['fillgaps'] = bp_final['bp_fillgaps']
    caltables[caltable_name]['solnorm'] = bp_final['bp_solnorm']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     bp_final['bp_combine'])
    caltables[caltable_name]['interp'] = bp_final['bp_interp']
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name)
    caltables[caltable_name]['gainfield'] = get_unique_field(
        caltables[caltable_name]['table'])
    logger.info('Bandpass_final BP {0}: {1}'.format(caltable_name, bptable))
    # Plots
    bptableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    bptableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_amp.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        bptableplot_phs,
                        gaintype='B',
                        calmode='p')

    logger.info('Bandpass_final BP phase plot: {0}'.format(bptableplot_phs))
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        bptableplot_amp,
                        gaintype='B',
                        calmode='ap')

    logger.info(
        'Bandpass_final BP amplitude plot: {0}'.format(bptableplot_amp))
    logger.info('End bandpass_final')
    # Apply calibration if requested:
    if eMCP['input_steps']['bandpass_final'] == 2:
        run_applycal(eMCP, caltables, step='bandpass_final')

    save_obj(caltables, caltables['calib_dir'] + 'caltables.pkl')
    msg = 'field={0}, combine={1}, solint={2}'.format(
        msinfo['sources']['bpcal'], bp_final['bp_combine'],
        bp_final['bp_solint'])
    eMCP = add_step_time('bandpass_final', eMCP, msg, t0)
    return eMCP, caltables


def gaincal_final(eMCP, caltables):
    logger.info('Start gaincal_final')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS:
    check_sources_in_ms(eMCP)
    gain_final = eMCP['defaults']['gaincal_final']

    # Compute p, ap for calibrators
    eMCP, caltables = gaincal_final_cals(eMCP, caltables)

    # Compute p, ap for targets (per scan)
    eMCP, caltables = gaincal_final_scan(eMCP, caltables)

    if eMCP['is_mixed_mode']:
        eMCP, caltables = gaincal_narrow(eMCP, caltables)
    # Apply calibration if requested:
    if eMCP['input_steps']['gaincal_final'] == 2:
        run_applycal(eMCP, caltables, step='gaincal_final')
    save_obj(caltables, caltables['calib_dir'] + 'caltables.pkl')
    logger.info('End gaincal_final')
    msg = 'p_solint={0}, ap_solint={1}'.format(gain_final['p_solint'],
                                               gain_final['ap_solint'])
    eMCP = add_step_time('gaincal_final', eMCP, msg, t0)
    return eMCP, caltables


def gaincal_narrow(eMCP, caltables, doplots=True):
    gain_final = eMCP['defaults']['gaincal_final']
    msinfo = eMCP['msinfo']
    msfile = eMCP['msinfo']['msfile']
    msfile_sp = eMCP['msinfo']['msfile_sp']
    spwmap_sp = eMCP['msinfo']['spwmap_sp']

    # 1 Phase offset wide-narrow
    caltable_name = gain_final['p_offset_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['p_offset_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = msinfo['sources']['ptcal']
    caltables[caltable_name]['solint'] = gain_final['p_offset_solint']
    caltables[caltable_name]['combine'] = gain_final['p_offset_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo,
                                               gain_final['p_offset_spw'])
    caltables[caltable_name]['spwmap'] = make_spwmap(
        caltables, gain_final['p_offset_combine'])
    caltables[caltable_name]['interp'] = gain_final['p_offset_interp']
    caltables[caltable_name]['minblperant'] = gain_final[
        'p_offset_minblperant']
    caltables[caltable_name]['minsnr'] = gain_final['p_offset_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal_narrow(msfile_sp, caltables, caltable_name, spwmap_sp)
    caltables[caltable_name]['gainfield'] = get_unique_field(
        caltables[caltable_name]['table'])
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
    # Plots
    if doplots:
        caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
            'inbase'] + '_' + caltable_name + '_phs.png'
        emplt.plot_caltable(caltables[caltable_name]['table'],
                            caltableplot_phs,
                            gaintype='G',
                            calmode='p')

        logger.info('{0} phase plot: {1}'.format(caltable_name,
                                                 caltableplot_phs))

    # 2 Narrow bandpass final calibration
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # Bandpass calibration
    caltable_name = gain_final['narrow_bp_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['narrow_bp_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['ptcal']
    caltables[caltable_name]['solint'] = gain_final['narrow_bp_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo,
                                               gain_final['narrow_bp_spw'])
    caltables[caltable_name]['combine'] = gain_final['narrow_bp_combine']
    caltables[caltable_name]['uvrange'] = gain_final['narrow_bp_uvrange']
    caltables[caltable_name]['fillgaps'] = gain_final['narrow_bp_fillgaps']
    caltables[caltable_name]['solnorm'] = gain_final['narrow_bp_solnorm']
    caltables[caltable_name]['spwmap'] = make_spwmap(
        caltables, gain_final['narrow_bp_combine'])
    caltables[caltable_name]['interp'] = gain_final['narrow_bp_interp']
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass_narrow(msfile_sp, caltables, caltable_name, spwmap_sp)
    caltables[caltable_name]['gainfield'] = get_unique_field(
        caltables[caltable_name]['table'])
    logger.info('Bandpass_final BP {0}: {1}'.format(caltable_name, bptable))
    # Plots
    bptableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    bptableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_amp.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        bptableplot_phs,
                        gaintype='B',
                        calmode='p')

    logger.info('Bandpass_final BP phase plot: {0}'.format(bptableplot_phs))
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        bptableplot_amp,
                        gaintype='B',
                        calmode='ap')

    logger.info(
        'Bandpass_final BP amplitude plot: {0}'.format(bptableplot_amp))
    # Apply calibration if requested:
    if eMCP['input_steps']['bandpass_final'] == 2:
        run_applycal(eMCP, caltables, step='bandpass_final')
    logger.info('End bandpass_final')
    return eMCP, caltables


def gaincal_final_cals(eMCP, caltables):
    gain_final = eMCP['defaults']['gaincal_final']
    calsources = eMCP['msinfo']['sources']['calsources']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # 1 Phase calibration
    caltable_name = gain_final['p_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['p_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = ','.join(
        np.unique(calsources.split(',')))
    caltables[caltable_name]['solint'] = gain_final['p_solint']
    caltables[caltable_name]['combine'] = gain_final['p_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_final['p_spw'])
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     gain_final['p_combine'])
    caltables[caltable_name]['interp'] = gain_final['p_interp']
    caltables[caltable_name]['minblperant'] = gain_final['p_minblperant']
    caltables[caltable_name]['minsnr'] = gain_final['p_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
    # Plots
    caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_phs,
                        gaintype='G',
                        calmode='p')

    logger.info('{0} phase plot: {1}'.format(caltable_name, caltableplot_phs))

    # 1 Amplitude calibration
    caltable_name = gain_final['ap_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['ap_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = select_calibrator(
        gain_final['ap_calibrator'], eMCP)
    caltables[caltable_name]['solint'] = gain_final['ap_solint']
    caltables[caltable_name]['combine'] = gain_final['ap_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_final['ap_spw'])
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,
                                                     gain_final['ap_combine'])
    caltables[caltable_name]['interp'] = gain_final['ap_interp']
    caltables[caltable_name]['minblperant'] = gain_final['ap_minblperant']
    caltables[caltable_name]['minsnr'] = gain_final['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


# Plots
    caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_phs,
                        gaintype='G',
                        calmode='p')

    logger.info('{0} phase plot: {1}'.format(caltable_name, caltableplot_phs))
    caltableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_amp.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_amp,
                        gaintype='G',
                        calmode='ap')

    logger.info('{0} amp plot: {1}'.format(caltable_name, caltableplot_amp))
    return eMCP, caltables


def gaincal_final_scan(eMCP, caltables):
    gain_final = eMCP['defaults']['gaincal_final']
    phscals = eMCP['msinfo']['sources']['phscals']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']

    # 1 Phase calibration on phasecal: Scan-averaged phase calibration
    caltable_name = gain_final['p_scan_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['p_scan_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = ','.join(np.unique(phscals.split(',')))
    caltables[caltable_name]['solint'] = gain_final['p_scan_solint']
    caltables[caltable_name]['combine'] = gain_final['p_scan_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo,
                                               gain_final['p_scan_spw'])
    caltables[caltable_name]['gainfield'] = caltables[caltable_name]['field']
    caltables[caltable_name]['spwmap'] = make_spwmap(
        caltables, gain_final['p_scan_combine'])
    caltables[caltable_name]['interp'] = gain_final['p_scan_interp']
    caltables[caltable_name]['minblperant'] = gain_final['p_scan_minblperant']
    caltables[caltable_name]['minsnr'] = gain_final['p_scan_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])

# Plots
    caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_phs,
                        gaintype='G',
                        calmode='p')

    # 2 Amplitude calibration on phasecal: scan-averaged amplitude solutions
    caltable_name = gain_final['ap_scan_tablename']
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir'] + caltables[
        'inbase'] + '_' + caltable_name
    caltables[caltable_name]['previous_cal'] = gain_final['ap_scan_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = select_calibrator(
        gain_final['ap_scan_calibrator'], eMCP, add_main=False)
    caltables[caltable_name]['solint'] = gain_final['ap_scan_solint']
    caltables[caltable_name]['combine'] = gain_final['ap_scan_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo,
                                               gain_final['ap_scan_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['phscals']
    caltables[caltable_name]['spwmap'] = make_spwmap(
        caltables, gain_final['ap_scan_combine'])
    caltables[caltable_name]['interp'] = gain_final['ap_scan_interp']
    caltables[caltable_name]['minblperant'] = gain_final['ap_scan_minblperant']
    caltables[caltable_name]['minsnr'] = gain_final['ap_scan_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '' and caltables[
            'Lo_dropout_scans'] != 'none':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])


# Plots
    caltableplot_phs = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_phs.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_phs,
                        gaintype='G',
                        calmode='p')

    logger.info('{0} phase plot: {1}'.format(caltable_name, caltableplot_phs))
    caltableplot_amp = caltables['plots_dir'] + 'caltables/' + caltables[
        'inbase'] + '_' + caltable_name + '_amp.png'
    emplt.plot_caltable(caltables[caltable_name]['table'],
                        caltableplot_amp,
                        gaintype='G',
                        calmode='ap')

    logger.info('{0} amp plot: {1}'.format(caltable_name, caltableplot_amp))
    return eMCP, caltables


def applycal_all(eMCP, caltables):
    t0 = datetime.datetime.utcnow()
    run_applycal(eMCP, caltables, step='applycal_all')
    if eMCP['is_mixed_mode']:
        spwmap_sp = eMCP['msinfo']['spwmap_sp']
        run_applycal_narrow(eMCP,
                            caltables,
                            step='applycal_all',
                            spwmap_sp=spwmap_sp)

    flag_statistics(eMCP, step='applycal_all')
    # Run statwt if requested
    if eMCP['defaults']['applycal_all']['run_statwt']:
        timebin = eMCP['defaults']['applycal_all']['statwt_timebin']
        msfile = eMCP['msinfo']['msfile']
        logger.info('Running statwt on {}'.format(msfile))
        logger.info('timebin: {}'.format(timebin))
        statwt(vis=msfile, timebin=timebin)
        find_casa_problems()
        if eMCP['is_mixed_mode']:
            msfile_sp = eMCP['msinfo']['msfile_sp']
            logger.info('Running statwt on narrow data {}'.format(msfile_sp))

            statwt(vis=msfile_sp, timebin=timebin)
            find_casa_problems()
    else:
        logger.info('statwt not selected')
    msg = ''
    eMCP = add_step_time('applycal_all', eMCP, msg, t0)
    return eMCP


def calc_Lo_drops(amp_mean, phscal_scans, threshold=0.5):
    # Sort amplitude values:
    a = np.sort(amp_mean[~np.isnan(amp_mean)])
    # Iterate to create two groups. Compute sum of std of both groups
    # That sum will be minimum when the two groups are divided by the
    # amplitude intersecting the two real distributions
    astd = np.array([0.0])
    for i in range(1, len(a)):
        astd = np.append(astd, np.std(a[:i]) + np.std(a[i:]))

    if astd[1:].min() > threshold * astd[1:].max():
        logger.info('No evidence of bimodality. Lo_drop_scans will be empty')
        lo_dropout_scans = []
    else:
        astd[0] = astd.max()
        idxmin = np.argmin(astd)
        separation_amp = a[idxmin]
        drops = amp_mean < separation_amp
        lo_dropout_scans = phscal_scans[drops]
    return lo_dropout_scans


def find_Lo_drops(msfile, phscals, eMCP):
    spws = eMCP['defaults']['flag_manual_avg']['Lo_spws']
    logger.info('Searching for possible Lo dropout scans on phscal scans')
    lo_dropout_scans = np.array([], dtype='int')
    scans, field_for_scan = find_fields_scans(msfile)
    for phscal in phscals:
        logger.info('Now searching for phasecal: {}'.format(phscal))
        phscal_scans = scans[field_for_scan == phscal]
        amp_mean, amp_std = find_Lo_amp(msfile, phscal, phscal_scans, eMCP,
                                        spws)
        threshold = eMCP['defaults']['flag_manual_avg']['Lo_threshold']
        lo_dropout_scans_i = calc_Lo_drops(amp_mean, phscal_scans, threshold)
        emplt.plot_Lo_drops(msfile, phscal_scans, scans, amp_mean,
                            lo_dropout_scans_i, phscal, eMCP)
        if len(lo_dropout_scans_i) > 0:
            logger.info('Potential dropout scans: '
                        '{0}'.format(','.join(
                            lo_dropout_scans_i.astype('str'))))
            min_scans = eMCP['defaults']['flag_manual_avg']['Lo_min_scans']
            if min_scans == '':
                min_scans = int(len(phscal_scans) * 0.4)
            if len(lo_dropout_scans_i) >= min_scans:
                lo_dropout_scans = np.hstack(
                    [lo_dropout_scans, lo_dropout_scans_i])
            else:
                logger.info('Less than {0}/{1} dropout scans, not considered '
                            'persistent drops'.format(min_scans,
                                                      len(phscal_scans)))
                lo_dropout_scans = []
        else:
            lo_dropout_scans = []
    return lo_dropout_scans


def remove_flagversion(msfile, versionname):
    flagmanager(vis=msfile, mode="delete", versionname=versionname)
    find_casa_problems()


def flag_statistics(eMCP, step):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    logger.info(line0)
    logger.info('Start flagstatistics')
    plots_obs_dir = './weblog/plots/plots_flagstats/'
    emutils.makedir(plots_obs_dir)
    logger.info('Running flagdata on {0}'.format(step))
    logger.info('mode="summary", action="calculate", antenna="*&*"'.format(
        msinfo['msfile']))
    run_flagstats(msfile, step)
    versionname = 'eMCP_{}'.format(step)
    remove_flagversion(msfile, versionname)
    logger.info('Saving flagtable in versionname="{0}"'.format(versionname))
    current_time = datetime.datetime.utcnow()

    flagmanager(vis=msfile,
                mode="save",
                versionname=versionname,
                comment=str(current_time))

    find_casa_problems()
    outfile = weblog_dir + 'plots/plots_flagstats/flagstats_{}.pkl'.format(
        step)
    logger.info('flagstats file saved to: {}'.format(outfile))
    logger.info('Flag statistics ready. Now plotting.')
    flag_stats = load_obj(outfile)
    emplt.plot_flagstatistics(flag_stats, msinfo, step)
    logger.info('End flagstatistics')
