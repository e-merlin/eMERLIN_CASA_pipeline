#!/usr/local/python
import os, subprocess
import numpy as np
import pickle
import glob
from Tkinter import *
import tkMessageBox
import sys, shutil
import copy
import getopt
import datetime
import collections
from eMERLIN_CASA_GUI import GUI_pipeline
from scipy import stats
import logging

import functions.weblog as emwlog
import functions.eMERLIN_CASA_plots as emplt

# CASA imports
from taskinit import *
from tasks import *
from recipes.setOrder import setToCasaOrder
from casac import casac
msmd = casac.msmetadata()

# Logging
logger = logging.getLogger('logger')

weblog_dir = './weblog/'
info_dir   = './weblog/info/'
calib_dir  = './weblog/calib/'
plots_dir  = './weblog/plots/'
logs_dir   = './logs/'
images_dir = './weblog/images/'

weblog_link= './'
info_link  = './info/'
calib_link = './calib/'
plots_link = './plots/'
images_link= './images/'

line0 = '-'*10

def check_in(pipeline_path):
    try:
        opts, arg = getopt.getopt(sys.argv[1:],'i:c:hg',['help','input=','gui'])
        logger.debug(sys.argv[1:])
    except getopt.GetoptError as err:
        logger.error(err)
        sys.exit(2)
    for o,a in opts:
        logger.debug('{0} {1}'.format(o,a))
        if o in ('-i','--input'):
            logger.info('Inputs from file: {}'.format(a))
            inputs = headless(a) ## read input file
#            inputs['quit'] = 0 ##needed to add to be compatible with GUI
        elif o in ('-g','--gui'):
            inputs = GUI_pipeline(pipeline_path).confirm_parameters() ## read input file
            logger.info('inputs from GUI: {}'.format(inputs))
        elif o in ('-h','--help'):
            logger.debug('help will be written soon')
            sys.exit('Closing pipeline eMCP')
        elif o == '-c':
            logger.debug('Executing!')
        else:
            assert False, "rerun with either headless -i or gui" #if none are specifed run GUI
    return inputs

def backslash_check(directory):
    if directory[-1] != '/':
        return directory+'/'
    else:
        return directory


def headless(inputfile):
    ''' Parse the list of inputs given in the specified file. (Modified from evn_funcs.py)'''
    logger.info('Parameters in inputs file:')
    INPUTFILE = open(inputfile, "r")
    control = collections.OrderedDict()
    # a few useful regular expressions
    newline = re.compile(r'\n')
    space = re.compile(r'\s')
    char = re.compile(r'\w')
    comment = re.compile(r'#.*')
    # parse the input file assuming '=' is used to separate names from values
    for line in INPUTFILE:
        if char.match(line):
            line = comment.sub(r'', line)
            line = line.replace("'", '')
            (param, value) = line.split('=')
            param = newline.sub(r'', param)
            param = param.strip()
            param = space.sub(r'', param)
            value = newline.sub(r'', value)
            value = value.replace(' ','').strip()
            valuelist = value.split(',')
            if len(valuelist) == 1:
                if valuelist[0] == '0' or valuelist[0]=='1' or valuelist[0]=='2':
                    valueout = int(valuelist[0])
                else:
                    valueout = str(valuelist[0])
            else:
                valueout = ','.join(valuelist)
            control[param] = valueout
            logger.info('{0:16s}: {1}'.format(param, valueout))
    return control

def Tkinter_select():
    root = Tkinter.Tk()
    root.withdraw()
    file = tkFileDialog.askdirectory(parent=root,mode='rb',title='Choose a file')
    if file != None:
            print file
    return file

def check_history(vis):
    tb.open(vis+'/HISTORY')
    x = tb.getcol('MESSAGE')
    y = [i for i, item in enumerate(x) if 'eMER_CASA_Pipeline:' in item]
    if len(y) == 0:
        print 'Measurement set has not been processed \n'
    else:
        print 'WARNING: Some pipeline processes have already been run'
        for i in range(len(y)):
            print x[y[i]]

def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass

def rmdir(pathdir,message='Deleted:'):
    if os.path.exists(pathdir):
        try:
            shutil.rmtree(pathdir)
            logger.info('{0} {1}'.format(message, pathdir))
        except:
            logger.debug('Could not delete: {0} {1}'.format(message, pathdir))
            pass

def rmfile(pathdir,message='Deleted:'):
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


def exit_pipeline(msinfo=''):
    if msinfo != '':
        logger.info('Something went wrong. Producing weblog and quiting')
        emwlog.start_weblog()
    sys.exit()


# Functions to save and load dictionaries
def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def prt_dict(d, pre=''):
    subdict = []
    for key in d.keys():
        if type(d[key]) == dict:
            subdict.append(key)
        else:
            print('{0:20s}: {1}'.format(pre+key, d[key]))
    if subdict != []:
        for key_inner in subdict:
            print(pre+key_inner)
            prt_dict(d[key_inner], pre=pre+'   ')

def prt_dict_tofile(d, tofilename=None, addfile='', pre=' '):
    if tofilename != None:
        f = open(tofilename, 'wb')
    else:
        f = addfile
    subdict = []
    for key in d.keys():
        if type(d[key]) in [collections.OrderedDict, dict]:
            subdict.append(key)
        else:
            f.write('{0:20s}: {1}\n'.format(pre+key, d[key]))
    if subdict != []:
        for key_inner in subdict:
            f.write('{}\n'.format(pre+key_inner))
            prt_dict_tofile(d[key_inner], addfile=f, pre=pre+pre)

def add_step_time(step, eMCP, msg, t0, doweblog=True):
    t1 = datetime.datetime.utcnow()
    timestamp = t1.strftime('%Y-%m-%d %H:%M:%S')
    delta_t_min = (t1-t0).total_seconds()/60.
    eMCP['steps'][step] = [timestamp, delta_t_min, msg]
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    os.system('cp casa_eMCP.log {}casa_eMCP.log.txt'.format(info_dir))
    if doweblog:
        emwlog.start_weblog(eMCP)
    return eMCP

def check_pipeline_conflict(eMCP, pipeline_version):
    try:
        eMCP['pipeline_version']
        if ((~new_run) and (eMCP['pipeline_version'] != pipeline_version)):
            logger.warning(
            'The log shows that different versions of the pipeline'
            ' has been executed. Please verify versions')
            logger.warning('Previous version: {0}. Current version {1}'.format(
            eMCP['pipeline_version'], pipeline_version))
    except:
        pass

def check_mixed_mode(vis,mode):
    logger.info('Check for mixed mode')
    tb.open(vis + '/SPECTRAL_WINDOW')
    bw_spw = np.array(tb.getcol('TOTAL_BANDWIDTH'))
    tb.close()
    if len(np.unique(bw_spw)) != 1:
        if mode == 'split':
            logger.info('Splitting continuum from spectral line')
            cont_spw = np.where(bw_spw==np.max(np.unique(bw_spw)))[0]
            print np.array2string(cont_spw, separator=',')[1:-1]
            split(vis=vis, outputvis=vis+'.continuum', spw=np.array2string(cont_spw, separator=',')[1:-1], datacolumn='data')
            spec_line = np.delete(bw_spw, cont_spw)
            logger.info('Splitting spectral line')
            for i in range(len(np.unique(spec_line))):
                spec_line_spw = np.where(bw_spw==np.unique(spec_line)[i])[0]
                split(vis=vis, outputvis=vis+'.sp{0}'.format(i), spw=np.array2string(spec_line_spw, separator=',')[1:-1],datacolumn='data')
                ms.writehistory(message='eMER_CASA_Pipeline: Spectral line split from {0}'.format(vis),msname=vis+'.sp{0}'.format(i))
            ms.writehistory(message='eMER_CASA_Pipeline: Spectral lines split from this ms',msname=vis)
            os.system('mv {0} {1}'.format(vis, vis+'.original'))
            os.system('mv {0} {1}'.format(vis+'.continuum', vis))
            logger.info('Will continue with continuum, original data is {0}'.format(vis+'.original'))
            return_variable = ''
        if mode == 'check':
            logger.info('MS is mixed mode. Please split')
            return_variable = True
    else:
        if mode == 'split':
            logger.info('Not mixed mode, continuing')
            return_variable = ''
        if mode == 'check':
            return_variable = False
    return return_variable

def check_band(msfile):
    # Take first frequency in the MS
    msmd.open(msfile)
    freq = msmd.chanfreqs(0)[0]/1e9
    msmd.done()
    band = ''
    if (freq > 1.2) and (freq < 1.7):
        band = 'L'
    elif (freq > 4) and (freq < 8):
        band = 'C'
    elif (freq > 22) and (freq < 24):
        band = 'K'
    else:
        logger.critical('Cannot determine band from frequency {}'.format(freq))
        exit_pipeline(eMCP)
    return band

def get_baselines(msfile):
    msmd.open(msfile)
    antennas0 = msmd.antennanames()
    baselines0 = msmd.baselines()
    msmd.close()
    baselines = []
    for i, a in enumerate(antennas0):
        for j, b in enumerate(antennas0):
            if j > i:
                baselines.append('{0}-{1}'.format(a, b))
    return np.array(baselines)

def get_dates(d):
    t_mjd   = d['axis_info']['time_axis']['MJDseconds']/60./60./24.
    t = np.array([mjdtodate(timei) for timei in t_mjd])
    return t_mjd, t


def mjdtodate(mjd):
    origin = datetime.datetime(1858,11,17)
    date = origin + datetime.timedelta(mjd)
    return date

def join_lists(x=[]):
    x1 = ','.join(x)
    x2 = set(x1.split(','))
    return ','.join(x2)

def user_sources(inputs):
    sources = collections.OrderedDict()
    sources['targets'] = inputs['targets']
    sources['phscals'] = inputs['phscals']
    sources['fluxcal'] = inputs['fluxcal']
    sources['bpcal']   = inputs['bpcal']
    sources['ptcal']   = inputs['ptcal']
    sources['calsources'] = join_lists([sources['phscals'], sources['fluxcal'],sources['bpcal'], sources['ptcal']])
    sources['calsources'] = sources['calsources'].strip(',')
    sources['maincal'] = join_lists([sources['fluxcal'],sources['bpcal'], sources['ptcal']])
    sources['allsources'] = join_lists([sources['calsources'], sources['targets']])
    sources['allsources'] = sources['allsources'].strip(',')
    sources['no_fluxcal'] = sources['allsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
    sources['cals_no_fluxcal'] = sources['calsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
    sources['targets_phscals'] = join_lists([sources['targets'],sources['phscals']])
    #logger.info('Targets:   {0}'.format(sources['targets']))
    #logger.info('Phasecals: {0}'.format(sources['phscals']))
    #logger.info('Fluxcal:   {0}'.format(sources['fluxcal']))
    #logger.info('Bandpass:  {0}'.format(sources['bpcal']))
    #logger.info('Pointcal:  {0}'.format(sources['ptcal']))
    return sources

def get_antennas(msfile):
    # Antenna list 
    msmd.open(msfile)
    antennas = msmd.antennanames()
    msmd.close()
    nice_order = ['Lo', 'Mk2', 'Pi', 'Da', 'Kn', 'De', 'Cm']
    antennas = [a for a in nice_order if a in antennas]
    #logger.info('Antennas in MS {0}: {1}'.format(msfile, antennas))
    return antennas

def get_obstime(msfile):
    # returns datetime object of first and last times in obs_id 0
    msmd.open(msfile)
    t_ini = mjdtodate(msmd.timerangeforobs(0)['begin']['m0']['value'])
    t_end = mjdtodate(msmd.timerangeforobs(0)['end']['m0']['value'])
    msmd.done()
    return t_ini, t_end

def get_obsfreq(msfile):
    # Returns freq of first channel, end chan, channel resolution
    # and number of channels (first spw) in GHz
    msmd.open(msfile)
    nspw = msmd.nspw()
    freq_ini = msmd.chanfreqs(0)[0]/1e9
    freq_end = msmd.chanfreqs(nspw-1)[-1]/1e9
    chan_res = msmd.chanwidths(0)[0]/1e9
    nchan = len(msmd.chanwidths(0))
    msmd.done()
    return freq_ini, freq_end, chan_res, nchan

def find_mssources(msfile):
    #mssources = ','.join(vishead(msfile,mode='list',listitems='field')['field'][0])
    msmd.open(msfile)
    mssources = ','.join(msmd.fieldnames())
    msmd.done()
    #logger.info('Sources in MS {0}: {1}'.format(msfile, mssources))
    return mssources

def find_source_intent(msinfo, cats=['targets', 'phscals', 'bpcal', 'fluxcal', 'ptcal']):
    fields_ms = msinfo['sources']['mssources'].split(',')
    return {source: ','.join([cat for cat in cats if source in msinfo['sources'][cat].split(',')]) for source in fields_ms}

def get_project(msfile):
    tb.open(msfile+'/OBSERVATION')
    project = tb.getcol('PROJECT')
    tb.close()
    return project[0]

def get_polarization(msfile):
    tb.open(msfile+'/FEED')
    polarization = tb.getcol('POLARIZATION_TYPE')
    tb.close()
    return ', '.join(polarization[:,0])

def get_directions(msfile):
    directions = collections.OrderedDict()
    msmd.open(msfile)
    field_names = msmd.namesforfields()
    msmd.done()
    for field_id, field_name in enumerate(field_names):
        vhead = vishead(msfile, mode = 'list', listitems = 'ptcs')
        ra_float = vhead['ptcs'][0]['r'+str(field_id+1)][0][0][0]*180./np.pi
        de_float = vhead['ptcs'][0]['r'+str(field_id+1)][1][0][0]*180./np.pi
        directions[field_name] = me.direction('J2000', '{}deg'.format(ra_float), '{}deg'.format(de_float))
    return directions


def get_distances(msfile, directions=''):
    if directions == '':
        directions = get_directions(msfile)
    msmd.open(msfile)
    field_names = msmd.namesforfields()
    msmd.done()
    separations = collections.OrderedDict()
    # Write all separations in a txt file
    sep_file = open(info_dir+'source_separations.txt', 'wb')
    for i in range(len(field_names)):
        for j in range(i+1, len(field_names)):
            f1 = field_names[i]
            f2 = field_names[j]
            separations[f1+'-'+f2] = me.separation(directions[f1],directions[f2])['value']
            sep_file.write('{0:10} {1:10} {2:7.2f}deg\n'.format(f1, f2, separations[f1+'-'+f2]))
    return separations

def get_integration_time(msfile, usemsmd=True):
    if usemsmd:
        msmd.open(msfile)
        first_scan_number = msmd.scannumbers()[0]
        int_time = msmd.exposuretime(first_scan_number)['value']
        msmd.done()
    else: # For CASA versions <5
        ms.open(msfile)
        axis_info = ms.getdata(['axis_info'],ifraxis=True)
        t_mjd   = axis_info['axis_info']['time_axis']['MJDseconds']
        ms.close()
        int_time = stats.mode(np.diff(t_mjd))[0][0]
    return int_time


def get_msinfo(eMCP, msfile, doprint=False):
    inputs = eMCP['inputs']
    logger.info('Found MS file: {0}'.format(msfile))
    logger.info('Reading ms file information for MS: {0}'.format(msfile))
    msinfo = collections.OrderedDict()
    msinfo['msfile'] = msfile
    msinfo['msfilename'] = os.path.splitext(msfile)[0].split('/')[-1]
    msinfo['project'] = get_project(msfile)
    msinfo['run'] = inputs['inbase']
    msinfo['sources'] = user_sources(inputs)
    msinfo['sources']['mssources'] = find_mssources(msfile)
    msinfo['sources']['source_intent'] = find_source_intent(msinfo)
    msinfo['antennas'] = get_antennas(msfile)
    msinfo['band'] = check_band(msfile)
    msinfo['baselines'] = get_baselines(msfile)
    msinfo['num_spw'] = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    t_ini, t_end = get_obstime(msfile)
    freq_ini, freq_end, chan_res, nchan = get_obsfreq(msfile)
    msinfo['t_ini'] = t_ini
    msinfo['t_end'] = t_end
    msinfo['freq_ini'] = freq_ini
    msinfo['freq_end'] = freq_end
    msinfo['int_time'] = get_integration_time(msfile)
    msinfo['chan_res'] = chan_res
    msinfo['nchan'] = nchan
    msinfo['innerchan'] = '{0:.0f}~{1:.0f}'.format(0.1*(nchan-nchan/512.), 0.9*(nchan-nchan/512.))
    msinfo['polarizations'] = get_polarization(msfile)
    msinfo['refant'] = define_refant(eMCP, msfile)
    msinfo['directions'] = get_directions(msfile)
    msinfo['separations'] = get_distances(msfile, directions=msinfo['directions'])
    # If eMCP['Lo_dropout_scans'] already there, it may have been recomputed in
    # flag_apriori step.
    try:
        msinfo['Lo_dropout_scans'] = eMCP['msinfo']['Lo_dropout_scans']
    except:
        msinfo['Lo_dropout_scans'] = eMCP['defaults']['flag_apriori']['Lo_dropout']
    logger.info('> Sources ({0}): {1}'.format(len(msinfo['sources']['mssources'].split(',')),
                                                 msinfo['sources']['mssources']))
    logger.info('> Number of spw: {0}'.format(msinfo['num_spw']))
    logger.info('> Channels per spw: {0}'.format(msinfo['nchan']))
    logger.info('> Itegration time {0:3.1f}s'.format(msinfo['int_time']))
    #save_obj(msinfo, info_dir + msinfo['msfilename']+'.msinfo.pkl')
    if doprint:
        prt_dict(msinfo)
    eMCP['msinfo'] = msinfo
    save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    return eMCP, msinfo, msfile

def get_unique_field(caltable):
    # If more than one field was used to generate the caltable but
    # combine='field' was used, only one field is included in the table. This
    # function will find which was actually used.
    tb.open(caltable)
    field_id = np.unique(tb.getcol('FIELD_ID'))[0]
    tb.close()
    tb.open(caltable+'/FIELD')
    unique_field = tb.getcol('NAME')[field_id]
    tb.close()
    return unique_field

def backup_table(caltable):
    try:
        shutil.rmtree(caltable+'_backup_missing')
    except OSError:
        pass
    shutil.copytree(caltable, caltable+'_backup_missing')
    logger.info('Backup of table {0} to {1}'.format(caltable, caltable+'_backup_missing'))

def remove_missing_scans(caltable, scans2flag):
    # Backup original table
    backup_table(caltable)
    tb.open(caltable+'/ANTENNA')
    anten = tb.getcol('NAME')
    tb.close()
    anten_Lo = np.argwhere('Lo')[0][0]
    tb.open(caltable, nomodify=False)
    scan_number = tb.getcol('SCAN_NUMBER')
    antenna1 = tb.getcol('ANTENNA1') == anten_Lo
    missing_rows = np.array([str(r) in scans2flag.replace(' ','').split(',') for r in scan_number])
    index_missing_rows = np.where(antenna1*missing_rows)[0]
    tb.removerows(index_missing_rows)
    logger.info('Removing Lo solutions for dropout scans from {0}: {1}'.format(caltable, scans2flag))
    tb.close()

def run_listobs(msfile):
    outfile = info_dir + msfile + '.listobs.txt'
    listobs(vis=msfile, listfile=outfile, overwrite=True)
    logger.info('Listobs file in: {0}'.format(outfile))

def decide_hanning(import_eM, msfile):
    run_hanning = import_eM['run_hanning']
    if run_hanning == 'auto':
        #logger.info('run_hanning = "auto" means check if data are L band.')
        band = check_band(msfile)
        if band == 'L':
            logger.info('L band dataset. Hanning smoothing will be executed.')
            apply_hanning = True
        else:
            logger.info('Dataset is not L band. Hanning smoothing not needed.')
            apply_hanning = False
    elif run_hanning > 1:
        apply_hanning = True
    elif run_hanning == 0:
        apply_hanning = False
    else:
        logger.warning('Not valid run_hanning parameter (use "auto", 0 or 1)')
        exit_pipeline()
    return apply_hanning


def mixed_mode(msfile):
    tb.open(msfile + '/SPECTRAL_WINDOW')
    bw_spw = np.array(tb.getcol('TOTAL_BANDWIDTH'))
    tb.close()
    is_mixed_mode = len(np.unique(bw_spw)) != 1
    if is_mixed_mode:
        cont_spw = ','.join(np.where(bw_spw==np.max(np.unique(bw_spw)))[0].astype('str'))
        line_spw =','.join(np.where(bw_spw!=np.max(np.unique(bw_spw)))[0].astype('str'))
        spw_separation = [cont_spw, line_spw] # Example ['0,1,2,3','4,5']
    else:
        spw_separation = ['']
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
    rmdir(eMCP['inputs']['inbase'] + '.ms')
    rmdir(eMCP['inputs']['inbase'] + '.mms')
    rmdir(eMCP['inputs']['inbase'] + '_sp.ms')
    rmdir(eMCP['inputs']['inbase'] + '_sp.mms')
    import_eM = eMCP['defaults']['import_eM']
    logger.info('Starting import eMERLIN fitsIDI procedure')
    fits_path = backslash_check(eMCP['inputs']['fits_path'])
    msfile_name = eMCP['inputs']['inbase']
    msg = 'first execution'
    t0 = datetime.datetime.utcnow()
    eMCP = add_step_time('start_pipeline', eMCP, msg, t0, doweblog=False)

    # importfitsIDI
    constobsid = import_eM['constobsid']
    scanreindexgap_s = import_eM['scanreindexgap_s']
    fitsfiles =[]
    for infile in os.listdir(fits_path):
        if infile.endswith('fits') or infile.endswith('FITS'):
            fitsfiles = fitsfiles + [fits_path+infile]
            logger.info('FITS file found to be imported: {0}'.format(infile))
    logger.info('Start importfitsIDI')
    t0 = datetime.datetime.utcnow()
    msfile0 = msfile_name+'_imported.ms'
    rmdir(msfile0)
    importfitsidi(fitsidifile=fitsfiles, vis=msfile0,
                  constobsid=constobsid, scanreindexgap_s=scanreindexgap_s)
    logger.info('Created file: {}'.format(msfile0))
    logger.info('End importfitsIDI')
    msg = 'constobsid={0}, scanreindexgap_s={1}'.format(constobsid,
                                                       scanreindexgap_s)
    eMCP['msfile'] = eMCP['inputs']['inbase']+'.ms'
    msfile = msfile0
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    eMCP = add_step_time('importfitsIDI', eMCP, msg, t0, doweblog=True)

    # mstransform
    datacolumn = 'data'
    if import_eM['antenna'] == '':
        antenna = '*&*'
    else:
        antenna = '*&*;'+import_eM['antenna']
        logger.info('Splitting antennas: {}'.format(antenna))
    timeaverage = import_eM['timeaverage']
    timebin = import_eM['timebin']
    chanaverage = import_eM['chanaverage']
    chanbin = import_eM['chanbin']
    usewtspectrum = import_eM['usewtspectrum']
    do_hanning = decide_hanning(import_eM, msfile0)
    do_ms2mms = import_eM['ms2mms']
    is_mixed_mode, spw_separation = mixed_mode(msfile0)
    eMCP['is_mixed_mode'] = is_mixed_mode
    ext_ms = {False:'.ms',True:'.mms'}
    msfile1 = msfile_name+'_transformed' + ext_ms[do_ms2mms]
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
    rmdir(msfile1)
    mstransform(vis = msfile0,
                outputvis = msfile1,
                keepflags = True,
                datacolumn = datacolumn,
                antenna = antenna,
                spw = spw_separation[0],
                timeaverage = timeaverage, timebin=timebin,
                chanaverage = chanaverage, chanbin=chanbin,
                usewtspectrum = usewtspectrum,
                hanning = do_hanning,
                createmms = do_ms2mms
                )
    logger.info('Transformed: {0} into {1}'.format(msfile0, msfile1))
    if do_hanning:
        ms.writehistory(message='Hanning smoothing applied',msname=msfile1)
    if is_mixed_mode: # No hanning and no channel average
        logger.info('Mixed mode data detected')
        msfile1_sp = msfile_name+'_transformed_sp' + ext_ms[do_ms2mms]
        logger.info('Running mstransform on spectral line data')
        rmdir(msfile1_sp)
        mstransform(vis = msfile0,
                outputvis = msfile1_sp,
                keepflags = True,
                datacolumn = datacolumn,
                antenna = antenna,
                spw = spw_separation[1],
                timeaverage = timeaverage, timebin=timebin,
                usewtspectrum = usewtspectrum,
                createmms = do_ms2mms
                )
        logger.info('Transformed: {0} into {1}'.format(msfile0, msfile1_sp))
    if os.path.isdir(msfile1) == True:
        rmdir(msfile0)
    else:
        logger.critical('Problem generating {}. Stopping ' \
                        'pipeline'.format(msfile1))
        exit_pipeline(eMCP)
    logger.info('End mstransform')
    msg = 'hanning={0}, createmms={1}'.format(do_hanning,
                                              do_ms2mms)
    if timeaverage:
        msg += ', timebin={0}'.format(timebin)
    if chanaverage:
        msg += ', chanbin={0}'.format(chanbin)
    if import_eM['antenna'] != '':
        msg += ', antenna="{}"'.format(antenna)
    msfile = msfile1
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    eMCP = add_step_time('mstransform', eMCP, msg, t0, doweblog=True)

    # FIXVIS
    msfile = eMCP['inputs']['inbase'] + ext_ms[do_ms2mms]
    logger.info('Start FIXVIS')
    t0 = datetime.datetime.utcnow()
    fixvis(vis = msfile1,
           outputvis = msfile,
           reuse = False)
    logger.info('Fixed {0} into {1}'.format(msfile1, msfile))
    run_listobs(msfile)
    if os.path.isdir(msfile) == True:
        rmdir(msfile1)
    else:
        logger.critical('Problem generating {}. Stopping ' \
                        'pipeline'.format(msfile))
        exit_pipeline(eMCP)
    if is_mixed_mode:
        msfile_sp = eMCP['inputs']['inbase'] + '_sp' + ext_ms[do_ms2mms]
        logger.info('Running FIXVIS on spectral line data')
        fixvis(vis = msfile1_sp,
           outputvis = msfile_sp,
           reuse = False)
        logger.info('Fixed {0} into {1}'.format(msfile1_sp, msfile_sp))
        run_listobs(msfile_sp)
        if os.path.isdir(msfile_sp) == True:
            rmdir(msfile1_sp)
        else:
            logger.critical('Problem generating FIXVIS ms. Stopping pipeline')
            exit_pipeline(eMCP)
    logger.info('End FIXVIS')
    msg = ''
    eMCP, msinfo, msfile = get_msinfo(eMCP, msfile)
    eMCP = add_step_time('fixvis', eMCP, msg, t0, doweblog=True)
    eMCP['msfile'] = msfile
    flag_statistics(eMCP, step='import')
    return eMCP


def run_aoflagger_fields(eMCP):
    """This version of the autoflagger iterates through the ms within the mms structure selecting individual fields. It uses pre-defined strategies. The np.unique in the loop below is needed for single source files. The mms thinks there are many filds (one per mms). I think it is a bug from virtualconcat."""
    fields = eMCP['defaults']['aoflagger']['fields']
    separate_bands = eMCP['defaults']['aoflagger']['separate_bands']
    pipeline_path = eMCP['pipeline_path']
    msfile = eMCP['msinfo']['msfile']

    logger.info('Start run_aoflagger_fields')
    t0 = datetime.datetime.utcnow()
    if separate_bands == True:
        logger.info('Bands will be processed separately.')
    elif separate_bands == False:
        logger.info('Bands will be processed all together.')
    else:
        logger.warning('separate_bands can only be True or False')
        exit_pipeline()
    # Check if aoflagger is available:
    aoflagger_available = check_command('aoflagger')
    if not aoflagger_available:
        logger.critical('aoflagger requested but not available.')
        logger.warning('Exiting pipeline.')
        exit_pipeline()
    # Check that version is at least 2.9+
    old_aoflagger = check_aoflagger_version()
    if old_aoflagger:
        logger.critical('aoflagger version <2.9 does not work correctly.')
        logger.warning('Exiting pipeline.')
        exit_pipeline()
    vis_fields = vishead(msfile,mode='list',listitems='field')['field'][0]
    fields_num = {f:i for i,f in enumerate(vis_fields)}
    if fields == 'all':
        fields = vis_fields
    else:
        fields = np.atleast_1d(fields)
    for field in np.unique(fields):
        # First, check if user has a new strategy for this field in the local folder.
        # If not, check if user has produced a new strategy for this field in the pipeline folder (for typical sources, etc).
        # If not, check for default strategies for this field
        # If nothing is found, just use the default strategy
        if os.path.isfile('./aoflagger_strategies/user/{0}.rfis'.format(field))==True:
            aostrategy = './aoflagger_strategies/user/{0}.rfis'.format(field)
        elif os.path.isfile(pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field))==True:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field)
        else:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format('default_faint')
        logger.info('Running AOFLagger for field {0} ({1}) using strategy {2}'.format(field,fields_num[field], aostrategy))
        if separate_bands:
            num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
            for b in range(num_spw):
                logger.info('Processing source {0}, band {1}'.format(field, b))
                flagcommand = 'time aoflagger -fields {2} -bands {3} -strategy {0} {1}'.format(aostrategy, msfile, fields_num[field], b)
                #os.system(flagcommand+' | tee -a pre-cal_flag_stats.txt')
                os.system(flagcommand)
            logger.info('Last AOFlagger command: {}'.format(flagcommand))
        else:
            logger.info('Processing source {0}, all bands'.format(field))
            flagcommand = 'time aoflagger -fields {2} -strategy {0} {1}'.format(aostrategy, msfile, fields_num[field])
            os.system(flagcommand+' | tee -a pre-cal_flag_stats.txt')
            logger.info('Last AOFlagger command: {}'.format(flagcommand))
        ms.writehistory(message='eMER_CASA_Pipeline: AOFlag field {0} with strategy {1}:'.format(field, aostrategy),msname=msfile)
    #flag_applied(flags, 'flagdata0_aoflagger')
    flag_statistics(eMCP, step='aoflagger')
    logger.info('End run_aoflagger_fields')
    msg = ''
    eMCP = add_step_time('aoflagger', eMCP, msg, t0)
    return eMCP

def check_command(command):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([command], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def check_aoflagger_version():
    logger.info('Checking AOflagger version')
    from subprocess import Popen, PIPE
    try:
        process = Popen(['aoflagger'], stdout=PIPE)
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

#def ms2mms(eMCP):
#    mode = eMCP['defaults']['ms2mms']['mode']
#    msfile = eMCP['msfile']
#    logger.info('Start ms2mms')
#    if mode == 'parallel':
#        partition(vis=msfile,outputvis=msfile[:-3]+'.mms',createmms=True,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
#"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
#,taql=None)
#        if os.path.isdir(msfile[:-3]+'.mms') == True:
#            rmdir(msfile)
#            rmdir(msfile+'.flagversions')
#        ms.writehistory(message='eMER_CASA_Pipeline: Converted MS to MMS for parallelisation',msname=msfile[:-3]+'.mms')
#
#    ## Need to use single if you need to aoflag the data later
#    if mode == 'single':
#        partition(vis=msfile,outputvis=msfile[:-3]+'.ms',createmms=False,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
#"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
#,taql=None)
#        if os.path.isdir(msfile[:-3]+'.ms') == True:
#           rmdir(msfile)
#           rmdir(msfile+'.flagversions')
#    run_listobs(msfile[:-3]+'.mms')
#    logger.info('End ms2mms')
#    msg = ''
#    eMCP = add_step_time('ms2mms', eMCP, msg)
#    return eMCP

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
            separation = float('{0:5.2f}'.format(separations[s1+'-'+s2]))
        except:
            try:
                separation = '{0:5.2f}'.format(separations[s2+'-'+s1])
            except:
                separation = 0.0
    if separation < 1.0:
        quacktime = 20.
    elif 1.0 <= separation < 2.0:
        quacktime = 25.
    elif 2.0 <= separation < 3.5:
        quacktime = 30.
    elif separation >= 3.5:
        quacktime = 35.
    else:
        quacktime = 0
    return quacktime

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
    msmd.open(msfile)
    nchan = len(msmd.chanwidths(0))
    msmd.done()
    # Subband edges
    channels_to_flag = '*:0~{0};{1}~{2}'.format(nchan/128-1, nchan-nchan/128, nchan-1)
    logger.info('MS has {} channels/spw'.format(nchan))
    logger.info('Flagging edge channels {0}'.format(channels_to_flag))
    flagdata(vis=msfile, mode='manual', spw=channels_to_flag)
    # Slewing (typical):
    do_quack = eMCP['defaults']['flag_apriori']['do_quack']
    if do_quack:
        ## Target and phase reference, 20 sec
        logger.info('Flagging first 20 sec of all sources.')
        flagdata(vis=msfile, mode='quack', quackinterval=20)
        # Main calibrators, 5 min
        bright_cal = join_lists([si for si in ['1331+305','1407+284','0319+415'] if si in
                      msinfo['sources']['mssources']])
        if bright_cal != '':
            logger.info('Flagging 5 min from bright calibrators')
            flagdata(vis=msfile, field=bright_cal, mode='quack', quackinterval=300)
        else:
            logger.warning('No main calibrators (1331+305, 1407+284, 0319+415) found in data set')
        for s1, s2 in zip(msinfo['sources']['targets'].split(','),
                          msinfo['sources']['phscals'].split(',')):
            sources_in_ms = msinfo['sources']['mssources'].split(',')
            missing_sources = [si for si in [s1,s2] if si not in sources_in_ms]
            if missing_sources == []:
                quacktime = find_quacktime(msinfo, s1, s2)
                if s1 != '':
                    logger.info('Flagging first {0} sec of target {1} and phasecal {2}'.format(quacktime, s1, s2))
                    flagdata(vis=msfile, field=','.join([s1,s2]), mode='quack', quackinterval=quacktime)
            else:
                logger.warning('Warning, source(s) {} not present in MS, will not flag this pair'.format(','.join(missing_sources)))
    else:
        logger.info('No quacking selected')
    # Search for Lo dropouts
    Lo_dropout_scans = eMCP['defaults']['flag_apriori']['Lo_dropout']
    if 'Lo' in antennas:
        if Lo_dropout_scans == 'none':
            eMCP['msinfo']['Lo_dropout_scans'] = ''
        elif Lo_dropout_scans == '':
            if msinfo['sources']['phscals'] != '':
                logger.info('Searching for Lo dropout scans')
                logger.info("To avoid this, set deafult Lo_dropout_scans = 'none' "
                            "when rerun flag_apriori")
                phscals = msinfo['sources']['phscals'].split(',')
                Lo_drop_list = find_Lo_drops(msfile,
                                             phscals, eMCP)
                eMCP['msinfo']['Lo_dropout_scans'] = ','.join(Lo_drop_list.astype('str'))
                logger.info('Flagging Lo dropout scans: '
                            '{0}'.format(eMCP['msinfo']['Lo_dropout_scans']))
                flagdata(vis=msfile, antenna='Lo', scan=eMCP['msinfo']['Lo_dropout_scans'])
        else:
            eMCP['msinfo']['Lo_dropout_scans'] = Lo_dropout_scans
    else:
        eMCP['msinfo']['Lo_dropout_scans'] = ''

    # Flag Lo-Mk2
    if 'Lo' in antennas and 'Mk2' in antennas:
        logger.info('Flagging Lo-Mk2 baseline')
        flagdata(vis=msfile, mode='manual', antenna='Lo*&Mk2*')
    flag_statistics(eMCP, step='apriori')
    msg = ''
    logger.info('End flagdata1_apriori')
    eMCP = add_step_time('flag_apriori', eMCP, msg, t0)
    return eMCP


def flagdata_manual(eMCP, run_name='flag_manual'):
    msfile = eMCP['msinfo']['msfile']
    if run_name == 'flag_manual':
        inpfile = './inputfg.flags'
    elif run_name == 'flag_manual_avg':
        inpfile = './inputfg_avg.flags'
    else:
        logger.warning('Wrong run_name specified')
        inpfile = ''
    logger.info('Start {}'.format(run_name))
    t0 = datetime.datetime.utcnow()
    if not os.path.isfile(inpfile) == True:
        logger.critical('Manual flagging step requested but cannot access file: {0}'.format(inpfile))
        logger.warning('Stopping pipeline at this step')
        exit_pipeline()
    logger.info('Applying manual flags from file: {0}'.format(inpfile))
    flagdata(vis=msfile, mode='list', inpfile=inpfile)
    flag_statistics(eMCP, step=run_name)
    logger.info('End {}'.format(run_name))
    msg = 'file={0}'.format(inpfile)
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
    logger.info("winsize={0}, timecutoff={1}, freqcutoff={2}, maxnpieces={3}".format(
        tfcrop['winsize'], tfcrop['timecutoff'], tfcrop['freqcutoff'],
        tfcrop['maxnpieces']))
    logger.info("usewindowstats='{0}', halfwin={1}, extendflags={2}".format(
        tfcrop['uwstats'], tfcrop['halfwin'], tfcrop['extendflags']))
    flagdata(vis=msfile,
             mode = 'tfcrop',
             field = eMCP['msinfo']['sources'][tfcrop['sources']],
             antenna = tfcrop['antenna'],
             scan = tfcrop['scan'],
             spw = tfcrop['spw'],
             correlation = tfcrop['correlation'],
             ntime = tfcrop['ntime'],
             combinescans = tfcrop['combinescans'],
             datacolumn = tfcrop['datacolumn'],
             winsize = tfcrop['winsize'],
             timecutoff= tfcrop['timecutoff'],
             freqcutoff = tfcrop['freqcutoff'],
             maxnpieces = tfcrop['maxnpieces'],
             usewindowstats = tfcrop['uwstats'],
             halfwin = tfcrop['halfwin'],
             extendflags = tfcrop['extendflags'],
             action = tfcrop['action'],
             display = tfcrop['display'],
             flagbackup = tfcrop['flagbackup'])
    if defaults == 'flag_target':
        msg = "mode={0}, maxnpieces={1}, timecutoff={2}, freqcutoff={3}".format(
                'tfcrop',
                tfcrop['maxnpieces'],
                tfcrop['timecutoff'],
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
    if defaults == 'flag_target':
        flag_statistics(eMCP, step='flag_target')
        msg = "mode={0}, ntime={1}, timedevscale={2}, freqdevscale={3}".format(
                'rflag',
                rflag['ntime'],
                rflag['timedevscale'],
                rflag['freqdevscale'])
        eMCP = add_step_time('flag_target', eMCP, msg, t0)
        logger.info('End flag_target')
    return eMCP


def define_refant(eMCP, msfile):
    if eMCP['inputs']['refant'] == 'compute':
        recompute = True
        logger.info('Forcing recompute of refant')
    elif eMCP['inputs']['refant'] == '':
        try:
            refant = eMCP['msinfo']['refant']
            logger.info('Refant already in eMCP["msinfo"], will not recompute')
            recompute = False
        except:
            logger.info('No refant specified. Will compute optimal')
            recompute = True
    else:
        recompute = False
        refant = eMCP['inputs']['refant']
    if recompute:
        logger.info('Recomputing best reference antenna')
        refant0 = eMCP['inputs']['refant']
        refant_user = refant0.replace(' ', '').split(',')
        antennas = get_antennas(msfile)
        refant_in_ms = (np.array([ri in antennas for ri in refant_user])).all()
        if not refant_in_ms:
            if refant0 != '' and eMCP['inputs']['refant'] != 'compute':
                logger.warning('Selected reference antenna(s) {0} not in MS! User selection will be ignored'.format(refant0))
            # Finding best antennas for refant
            logger.info('Estimating best reference antenna.')
            logger.info('To avoid this, set a reference antenna in the inputs file.')
            refant, refant_pref = find_refant(msfile, field=eMCP['inputs']['bpcal'])
        else:
            refant = ','.join(refant_user)
    logger.info('Refant in eMCP: {}'.format(refant))
    return refant

def find_refant(msfile, field, spw='3'):
    antennas = get_antennas(msfile)
    v = {}
    snr = {}
    for anten in antennas:
        try:
            v[anten] = visstat(msfile, field='1407+284', antenna='{}&*'.format(anten), axis='amplitude', spw=str(spw), correlation='RR,LL')
            snr[anten] = v[anten]['DATA_DESC_ID=3']['median']/v[anten]['DATA_DESC_ID=3']['stddev']
        except:
            v[anten] = None
            snr[anten] = 0.0
            logger.warning('No data found for bpcal for antenna: '
                           '{}'.format(anten))
    for i in np.argsort(snr.values())[::-1]:
        logger.info('{0:3}: {1:3.1f}'.format(snr.keys()[i], snr.values()[i]))
    pref_ant = [snr.keys()[i] for i in np.argsort(snr.values())][::-1]
    snr_sorted = [snr.values()[i] for i in np.argsort(snr.values())][::-1]
    if 'Lo' in antennas:
        priorities = ['Pi','Da','Kn']
        secondary = ['Pi','Da','Kn', 'Cm', 'De']
    else:
        priorities = ['Mk2','Pi','Da']
        secondary = ['Mk2','Pi','Da','Kn', 'Cm', 'De']
    u, ind = np.unique([a for a in pref_ant[:3] if a in priorities] +\
                       [a for a  in pref_ant if a in secondary],
                       return_index=True)
    refant = ','.join(u[np.argsort(ind)])
    refant_pref = collections.OrderedDict()
    refant_pref['snr'] = snr_sorted
    refant_pref['sorted_antennas'] = pref_ant
    refant_pref['refant'] = refant
    logger.info('Preference antennas refant = {0}'.format(refant))
    return refant, refant_pref


def plot_caltable(msinfo, caltable, plot_file, xaxis='', yaxis='', title='',
                  ymin=-1, ymax=-1, coloraxis='spw', symbolsize=8):
    num_anten = len(msinfo['antennas'])
    gridcols = 1
    showgui=False
    for i, anten in enumerate(msinfo['antennas']):
        if i == len(msinfo['antennas'])-1:
            plotfile = plot_file
        else:
            plotfile = ''
        plotms(vis=caltable['table'], xaxis=xaxis, yaxis=yaxis, title='{0} {1}'.format(title, anten),
               gridrows=num_anten, gridcols=gridcols, rowindex=i, colindex=0, plotindex=i,
               #timerange='{}~{}'.format(msinfo['t_ini'].time(), msinfo['t_end'].time()),
               antenna = str(anten),
               xselfscale = True, xsharedaxis = True, coloraxis = coloraxis, plotrange=[-1,-1,ymin, ymax],
               plotfile = plotfile, expformat = 'png', customsymbol = True, symbolshape = 'circle', symbolsize=symbolsize,
               width=1000, height=240*num_anten, clearplots=False, overwrite=True, showgui=showgui)


def saveflagstatus(eMCP):
    msinfo = eMCP['msinfo']
    logger.info('Starting saveflagstatus')
    t0 = datetime.datetime.utcnow()
    logger.info('Saving current flagging status to versionname=\'initialize_flags\'')
    flagmanager(msinfo['msfile'], mode='save', versionname='initialize_flags',
             comment='Restore this version to restart calibration without the flags produced by the calibration',
             merge='replace')
    msg = 'versionname=initialize_flags'
    eMCP = add_step_time('save_flags', eMCP, msg, t0)
    return eMCP


def restoreflagstatus(eMCP):
    msinfo = eMCP['msinfo']
    t0 = datetime.datetime.utcnow()
    logger.info('Starting restoreflagstatus')
    logger.info('Restoring flagging status in versionname=\'initialize_flags\'')
    flagmanager(msinfo['msfile'], mode='restore', versionname='initialize_flags',
             merge='replace')
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
        exit_pipeline()
    # Check that all sources are in the MS:
    sources_not_in_msfile = [s for s in
                             sources['allsources'].split(',')
                             if s not in mssources.split(',')]
    if len(sources_not_in_msfile) > 0:
        fields = mssources
        logger.critical('Fields {} not present in MS but listed in ' \
                       'inputs file.'.format(','.join(sources_not_in_msfile)))
        logger.warning('Stopping pipeline at this step')
        exit_pipeline()

def check_table_exists(caltables, tablename):
    try:
        if os.path.isdir(caltables[tablename]['table']):
            it_exists = True
        else:
            it_exists = False
    except:
            it_exists = False
    if it_exists == False:
        logger.critical('Calibration table {} required but ' \
                        'not available'.format(tablename))
        logger.warning('Stopping pipeline at this step')
        exit_pipeline()


### Run CASA calibration functions

def run_split(eMCP):
    logger.info(line0)
    logger.info('Start average')
    t0 = datetime.datetime.utcnow()
    width = eMCP['defaults']['average']['width']
    msfile = eMCP['msinfo']['msfile']
    sources = eMCP['msinfo']['sources']
    timebin = '{}s'.format(eMCP['inputs']['average'])
    datacolumn = eMCP['defaults']['average']['datacolumn']
    scan = eMCP['defaults']['average']['scan']
    antenna = eMCP['defaults']['average']['antenna']
    timerange = eMCP['defaults']['average']['timerange']
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    fields = sources['allsources']
    name = '.'.join(msfile.split('.')[:-1])
    exte = ''.join(msfile.split('.')[-1])
    outputmsfile = name+'_avg.'+exte
    rmdir(outputmsfile)
    rmdir(outputmsfile+'.flagversions')
    logger.info('Input MS: {0}'.format(msfile))
    logger.info('Output MS: {0}'.format(outputmsfile))
    logger.info('width={0}, timebin={1}'.format(width, timebin))
    logger.info('Fields: {0}'.format(fields))
    logger.info('Data column: {0}'.format(datacolumn))
    split(vis=msfile, outputvis=outputmsfile, field=fields,
          timerange=timerange, scan=scan, antenna=antenna,
          timebin=timebin,  width=width,
          datacolumn=datacolumn, keepflags=True)
    run_listobs(outputmsfile)
    logger.info('End average')
    msg = 'width={0}, timebin={1}, datacolumn={2}'.format(width, timebin,
                                                          datacolumn)
    eMCP = add_step_time('average', eMCP, msg, t0)
    return eMCP


def run_initialize_models(eMCP):
    logger.info(line0)
    logger.info('Start init_models')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    msfile = eMCP['msinfo']['msfile']
    init_models = eMCP['defaults']['init_models']
    models_path = eMCP['pipeline_path']+init_models['calibrator_models']
    fluxcal = eMCP['msinfo']['sources']['fluxcal']
    logger.info('Deleting model of all sources')
    delmod(vis=msfile, otf=True, scr=True) #scr to delete MODEL column
    # Check dataset frequency:
    band = check_band(msfile)
    if band == 'C':
        model_3C286 = models_path+'3C286_C.clean.model.tt0'
        logger.info('Dataset is band C. Using C band model of 3C286')
    elif band == 'L':
        model_3C286 = models_path+'1331+305.clean.model.tt0'
        logger.info('Dataset is band L. Using L band model of 3C286')
    else:
        logger.warning('No model available!')
        model_3C286 = ''
    logger.info('Initializing model for 3C286')
    logger.info('Model {0}'.format('./'+'/'.join(model_3C286.split('/')[-2:])))
    if fluxcal != '1331+305':
        logger.warning('Using a model for 3C286 (1331+305) but your flux calibrator source is: {0}. Model may be wrong for that source'.format(fluxcal))
    setjy(vis=msfile, field=fluxcal,
          model=model_3C286, scalebychan=True, usescratch=True)
    logger.info('End init_models')
    msg = ''
    eMCP = add_step_time('init_models', eMCP, msg, t0)
    return eMCP


def initialize_cal_dict(inputs, eMCP):
    # All the calibration steps will be saved in the dictionary caltables.pkl
    # located in the calib directory. If it does not exist a new one is created.
    any_calsteps = ['bandpass', 'initial_gaincal','fluxscale','bandpass_sp','gain_amp_sp','applycal_all']
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
        logger.info('New caltables dictionary created. Saved to: {0}'.format(calib_dir+'caltables.pkl'))
        caltables['Lo_dropout_scans'] = eMCP['msinfo']['Lo_dropout_scans']
        caltables['refant'] = msinfo['refant']
        caltables['refantmode'] = eMCP['defaults']['global']['refantmode']
        save_obj(caltables, calib_dir+'caltables.pkl')
        return caltables


def run_gaincal(msfile, caltables, caltable_name):
    logger.info(line0)
    rmdir(caltables[caltable_name]['table'])
    logger.info('Running gaincal to generate: {0}'.format(caltables[caltable_name]['name']))
    logger.info('Field(s) = {0}, gaintype = {1}, calmode = {2}'.format(
                caltables[caltable_name]['field'],
                caltables[caltable_name]['gaintype'],
                caltables[caltable_name]['calmode']))
    logger.info('solint = {0}, spw = {1},  combine = {2}'.format(
                caltables[caltable_name]['solint'],
                caltables[caltable_name]['spw'],
                caltables[caltable_name]['combine']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['gainfield'] if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(caltables[caltable_name]['table']))
    # Run CASA task gaincal
    gaincal(vis=msfile,
            caltable  = caltables[caltable_name]['table'],
            field     = caltables[caltable_name]['field'],
            gaintype  = caltables[caltable_name]['gaintype'],
            calmode   = caltables[caltable_name]['calmode'],
            solint    = caltables[caltable_name]['solint'],
            combine   = caltables[caltable_name]['combine'],
            spw       = caltables[caltable_name]['spw'],
            refant    = caltables['refant'],
            refantmode= caltables['refantmode'],
            gaintable = gaintable,
            gainfield = gainfield,
            interp    = interp,
            spwmap    = spwmap,
            minblperant= caltables[caltable_name]['minblperant'],
            minsnr    = caltables[caltable_name]['minsnr'])
#    if caltables[caltable_name]['calmode'] == 'p' and  caltables[caltable_name]['gaintype'] == 'G':
#        refant = caltables['refant'].split(',')[0]
#        logger.info('Running rerefant to refant'
#                    ' {}'.format(refant))
#        rerefant(vis=msfile, tablein=caltables[caltable_name]['table'],
#                 refantmode='strict', refant =refant)
    logger.info('caltable {0} in {1}'.format(caltables[caltable_name]['name'],
                                              caltables[caltable_name]['table']))

def run_bandpass(msfile, caltables, caltable_name, minblperant=3, minsnr=2):
    logger.info(line0)
    rmdir(caltables[caltable_name]['table'])
    logger.info('Running bandpass to generate: {0}'.format(caltables[caltable_name]['name']))
    logger.info('Field(s) {0}, solint = {1}, spw = {2}, combine = {3}, solnorm = {4}'.format(
                 caltables[caltable_name]['field'],
                 caltables[caltable_name]['solint'],
                 caltables[caltable_name]['spw'],
                 caltables[caltable_name]['combine'],
                 caltables[caltable_name]['solnorm']))
    #logger.info('uvrange = {0}'.format(caltables[caltable_name]['uvrange']))
    # Previous calibration
    previous_cal = caltables[caltable_name]['previous_cal']
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['gainfield'] if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(caltables[caltable_name]['table']))
    # Run CASA task bandpass
    bandpass(vis=msfile,
             caltable  = caltables[caltable_name]['table'],
             field     = caltables[caltable_name]['field'],
             solint    = caltables[caltable_name]['solint'],
             combine   = caltables[caltable_name]['combine'],
             spw       = caltables[caltable_name]['spw'],
             solnorm   = caltables[caltable_name]['solnorm'],
             uvrange   = caltables[caltable_name]['uvrange'],
             fillgaps  = caltables[caltable_name]['fillgaps'],
             refant    = caltables['refant'],
             gaintable = gaintable,
             gainfield = gainfield,
             interp    = interp,
             spwmap    = spwmap,
             minblperant=minblperant,
             minsnr=minsnr)
    logger.info('caltable {0} in {1}'.format(caltables[caltable_name]['name'],
                                             caltables[caltable_name]['table']))


def smooth_caltable(msfile, tablein, plotdir, caltable='', field='', smoothtype='median', smoothtime=120.):
    logger.info('Smoothing table: {0}, field {1}, smoothtype {2}, smoothtime {3}'.format(tablein, field, smoothtype, smoothtime))
    basename = os.path.basename(tablein)
    caltableplot_phs = plotdir + basename +'_phs.png'
    caltableplot_amp = plotdir + basename +'_amp.png'
    if caltable=='':
        os.system('mv {0} {1}'.format(caltableplot_phs, plotdir + basename +'_phs_pre_smooth.png'))
        os.system('mv {0} {1}'.format(caltableplot_amp, plotdir + basename +'_amp_pre_smooth.png'))
    smoothcal(vis=msfile, tablein=tablein, caltable=tablein+'smooth', field='', smoothtype='median', smoothtime=60*20.)
    logger.info('Pre-smoothing table saved to: {0}'.format(tablein+'_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein, tablein+'_pre_smooth'))
    os.system('mv {0} {1}'.format(tablein+'smooth', tablein))
    #sub_plot = int('{}21'.format(int(len(msinfo['antennas'])/2.+0.5)))
    #plotcal(caltable=tablein, xaxis='time', yaxis='phase', subplot=sub_plot, iteration='antenna', showgui=False, figfile=caltableplot_phs, fontsize = 8, plotrange=[-1,-1,-180,180])
    #plotcal(caltable=tablein, xaxis='time', yaxis='amp', subplot=sub_plot, iteration='antenna', showgui=False, figfile=caltableplot_amp, fontsize = 8, plotrange=[-1,-1,-1,-1])
    return


def run_applycal(eMCP, caltables, step, dotarget=False):
    logger.info(line0)
    logger.info('Start applycal')
    logger.info('Applying calibration up to step: {}'.format(step))
    previous_cal = eMCP['defaults'][step]['apply_calibrators']
    previous_cal_targets= eMCP['defaults'][step]['apply_targets']
    msfile = eMCP['msinfo']['msfile']
    sources = eMCP['msinfo']['sources']
    # 1 correct non-target sources:
    logger.info('Applying calibration to calibrator sources')
    logger.info('Fields: {0}'.format(sources['calsources']))
    # Check if tables exist:
    for table_i in previous_cal:
        check_table_exists(caltables, table_i)
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['gainfield'] if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))<2 else '' for p in previous_cal]

    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    applycal(vis=msfile,
             field = sources['calsources'],
             gaintable = gaintable,
             gainfield = gainfield,
             interp    = interp,
             spwmap    = spwmap)

    # 2 correct targets
    if previous_cal_targets != []:
        previous_cal_targets= eMCP['defaults'][step]['apply_targets']
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
                gaintable = [caltables[p]['table'] for p in previous_cal_targets]
                interp    = [caltables[p]['interp'] for p in previous_cal_targets]
                spwmap    = [caltables[p]['spwmap'] for p in previous_cal_targets]
                gainfield = [caltables[p]['gainfield'] if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))<2 else phscal for p in previous_cal_targets]
                logger.info('Field: {0}. Phase calibrator: {1}'.format(s, phscal))
                logger.info('Previous calibration applied: {0}'.format(str(previous_cal_targets)))
                logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
                logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
                logger.info('Previous calibration interp: {0}'.format(str(interp)))
                applycal(vis=msfile,
                         field = s,
                         gaintable = gaintable,
                         gainfield = gainfield,
                         interp    = interp,
                         spwmap    = spwmap)
            else:
                logger.warning('Source {} is not phase-referenced'.format(s))
    else:
        logger.info('Not applying calibration to target sources at this stage')

    logger.info('End applycal')



### Calibration steps

def initial_bp_cal(eMCP, caltables):
    logger.info('Start initial_bpcal')
    t0 = datetime.datetime.utcnow()
    # Pass 1
    logger.info('Starting pass 1 of initial_bpcal')
    eMCP, caltables = run_bpcal(eMCP, caltables, doplots=False)

    # Apply solutions
    run_applycal(eMCP, caltables, step = 'bandpass')

    # Flagging
    eMCP = flagdata_tfcrop(eMCP, defaults='bandpass')

    # Pass 2
    logger.info('Starting pass 2 of initial_bpcal')
    eMCP, caltables = run_bpcal(eMCP, caltables)

    # Apply calibration if requested:
    if eMCP['inputs']['bandpass'] == 2:
        run_applycal(eMCP, caltables, step='bandpass')

    flag_statistics(eMCP, step='initial_bpcal')
    save_obj(caltables, caltables['calib_dir']+'caltables.pkl')
    bp = eMCP['defaults']['bandpass']
    msg = 'field={0}, combine={1}, solint={2}'.format(
                    eMCP['msinfo']['sources']['bpcal'],
                    bp['bp_combine'],
                    bp['bp_solint'])
    eMCP = add_step_time('bandpass', eMCP, msg, t0)
    logger.info('End initial_bpcal')
    return eMCP, caltables


def run_bpcal(eMCP, caltables, doplots=True):
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    bp = eMCP['defaults']['bandpass']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # 0 Delay calibration of bpcal
    caltable_name = bp['delay_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = bp['delay_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = bp['delay_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['delay_spw'])
    caltables[caltable_name]['combine'] = bp['delay_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,bp['delay_combine'])
    caltables[caltable_name]['interp'] = bp['delay_interp']
    caltables[caltable_name]['minblperant'] = bp['delay_minblperant']
    caltables[caltable_name]['minsnr'] = bp['delay_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Delay calibration {0}: {1}'.format(caltable_name, caltable))
    # Plots
    if doplots:
        caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_1.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Delay',
                      xaxis='time', yaxis='delay', ymin=-1, ymax=-1, coloraxis='corr', symbolsize=8)
        logger.info('Delay plot in: {0}'.format(caltableplot))

    # 1 Phase calibration
    caltable_name = bp['phase_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = bp['phase_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = bp['phase_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['phase_spw'])
    caltables[caltable_name]['combine'] = bp['phase_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,bp['phase_combine'])
    caltables[caltable_name]['interp'] = bp['phase_interp']
    caltables[caltable_name]['minblperant'] = bp['phase_minblperant']
    caltables[caltable_name]['minsnr'] = bp['phase_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Bandpass0 phase calibration {0}: {1}'.format(caltable_name,caltable))
    # Plots
    if doplots:
        caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Phase',
                      xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=3)
        logger.info('BP0_p phase plot: {0}'.format(caltableplot))
    logger.info('Apply phase calibration flags to bandpass, applymode=flagonly')
    applycal(vis=msfile, gaintable=caltables[caltable_name]['table'],
             field=msinfo['sources']['bpcal'],
             applymode='flagonly')

    # 2 Amplitude calibration
    caltable_name = bp['ap_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = bp['ap_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['solint'] = bp['ap_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['ap_spw'])
    caltables[caltable_name]['combine'] = bp['ap_combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,bp['ap_combine'])
    caltables[caltable_name]['interp'] = bp['ap_interp']
    caltables[caltable_name]['minblperant'] = bp['ap_minblperant']
    caltables[caltable_name]['minsnr'] = bp['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Bandpass0 amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    if doplots:
        caltableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot_phs, title='Phase',
                      xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=5)
        logger.info('BP0_ap phase plot: {0}'.format(caltableplot_phs))
        caltableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot_amp, title='Amp',
                      xaxis='time', yaxis='amp', ymin=-1, ymax=-1, coloraxis='spw', symbolsize=5)
        logger.info('BP0_ap amp plot: {0}'.format(caltableplot_amp))

    # 3 Bandpass calibration
    caltable_name = bp['bp_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = bp['bp_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['solint'] = bp['bp_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp['bp_spw'])
    caltables[caltable_name]['combine'] = bp['bp_combine']
    caltables[caltable_name]['uvrange'] = bp['bp_uvrange']
    caltables[caltable_name]['fillgaps'] = bp['bp_fillgaps']
    caltables[caltable_name]['solnorm'] = bp['bp_solnorm']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,bp['bp_combine'])
    caltables[caltable_name]['interp'] = bp['bp_interp']
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name)
    caltables[caltable_name]['gainfield'] = get_unique_field(caltables[caltable_name]['table'])
#    logger.info('Bandpass0 {0}: {1}'.format(caltable_name,bptable))
    # Plots
    if doplots:
        bptableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
        bptableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
        plot_caltable(msinfo, caltables[caltable_name], bptableplot_phs, title='Bandpass phase',
                      xaxis='freq', yaxis='phase', ymin=-180, ymax=180, coloraxis='corr', symbolsize=5)
        logger.info('BP0 phase plot: {0}'.format(bptableplot_phs))
        plot_caltable(msinfo, caltables[caltable_name], bptableplot_amp, title='Bandpass amp',
                      xaxis='freq', yaxis='amp', ymin=-1, ymax=-1, coloraxis='corr', symbolsize=5)
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

    # Initial gain calibration #
    eMCP, caltables = gain_p_ap(eMCP, caltables, doplots=False)
    run_applycal(eMCP, caltables, step = 'initial_gaincal')

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

    # Initial gain calibration #
    eMCP, caltables = gain_p_ap(eMCP, caltables)

    # Apply calibration if requested:
    if eMCP['inputs']['initial_gaincal'] == 2:
        run_applycal(eMCP, caltables, step='initial_gaincal')
    flag_statistics(eMCP, step='initial_gaincal')
    save_obj(caltables, caltables['calib_dir']+'caltables.pkl')
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


def solve_delays(eMCP, caltables, doplots=True):
    logger.info('Starting solve_delays')
    #t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    delay = eMCP['defaults']['initial_gaincal']['delay']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    caltable_name = delay['tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = delay['prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['calsources']
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = delay['solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, delay['spw'])
    caltables[caltable_name]['combine'] = delay['combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['calsources']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,delay['combine'])
    caltables[caltable_name]['interp'] = delay['interp']
    caltables[caltable_name]['minblperant'] = delay['minblperant']
    caltables[caltable_name]['minsnr'] = delay['minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Delay calibration {0}: {1}'.format(caltable_name, caltable))
    # Plots
    if doplots:
        # 1 No range
        caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_1.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Delay',
                      xaxis='time', yaxis='delay', ymin=-1, ymax=-1, coloraxis='field', symbolsize=8)
        caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_2.png'
#        plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Delay',
#                      xaxis='time', yaxis='delay', ymin=-20, ymax=20, coloraxis='field', symbolsize=8)
        logger.info('Delay plot: {0}'.format(caltableplot))
        logger.info('End solve_delays')
    return eMCP, caltables



def run_fringefit(msfile, caltables, caltable_name, minblperant=3, minsnr=2, smodel=[]):
    previous_cal = caltables[caltable_name]['previous_cal']
    rmdir(caltables[caltable_name]['table'])
    logger.info('Running fringefit to generate: {0}'.format(caltables[caltable_name]['name']))
    logger.info('Field(s) = {0}, zerorates = {1}'.format(
                caltables[caltable_name]['field'],
                caltables[caltable_name]['zerorates']))
    logger.info('solint = {0}, spw = {1},  combine = {2}'.format(
                caltables[caltable_name]['solint'],
                caltables[caltable_name]['spw'],
                caltables[caltable_name]['combine']))
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['gainfield'] if len(np.atleast_1d(caltables[p]['gainfield'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating cal table: {0}'.format(caltables[caltable_name]['table']))
    # Run CASA task fringefit
    fringefit(vis=msfile,
            caltable  = caltables[caltable_name]['table'],
            field     = caltables[caltable_name]['field'],
            solint    = caltables[caltable_name]['solint'],
            combine   = caltables[caltable_name]['combine'],
            spw       = caltables[caltable_name]['spw'],
            refant    = caltables['refant'],
            antenna   = '*&*',
            gaintable = gaintable,
            gainfield = gainfield,
            interp    = interp,
            spwmap    = spwmap,
            zerorates = caltables[caltable_name]['zerorates'])
    logger.info('caltable {0} in {1}'.format(caltables[caltable_name]['name'],
                                              caltables[caltable_name]['table']))

def delay_fringefit(eMCP, caltables):
    # Currently does not combine spws because that is not correctly implemented
    # in the pre-release version of task fringefit
    logger.info('Starting delay_fringefit')
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    delay = eMCP['defaults']['initial_gaincal']['delay']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    caltable_name = delay['tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = delay['prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['calsources']
    caltables[caltable_name]['solint'] = delay['solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, delay['spw'])
    caltables[caltable_name]['combine'] = delay['combine']
    caltables[caltable_name]['gainfield'] = msinfo['sources']['calsources']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,delay['combine'])
    caltables[caltable_name]['interp'] = delay['interp']
    caltables[caltable_name]['zerorates'] = delay['zerorates']
    caltables[caltable_name]['minblperant'] = delay['minblperant']
    caltables[caltable_name]['minsnr'] = delay['minsnr']
    caltable = caltables[caltable_name]['table']

    # Calibration
    run_fringefit(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Fringe calibration {0}: {1}'.format(caltable_name, caltable))
    # Plots
    # 1 Phases
    caltableplot =  caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Phase',
                  xaxis='time', yaxis='phase', ymin=-180, ymax=-180, coloraxis='field')
    logger.info('Fringe calibration (phases) plot in: {0}'.format(caltableplot))
    # 2 Delays
    caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_dela.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Phase',
                  xaxis='time', yaxis='delay', ymin=-1, ymax=-1, coloraxis='field', symbolsize=8)
    logger.info('Fringe calibration (delays) plot in: {0}'.format(caltableplot))
    caltableplot =  caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_dela2.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Phase',
                  xaxis='time', yaxis='delay', ymin=-20, ymax=-20, coloraxis='field', symbolsize=8)
    logger.info('Fringe calibration (delays range) plot in: {0}'.format(caltableplot))
    ## 3 Rates (they are zeroed)
    #caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_rate.png'
    #plotcal(caltable=caltable,xaxis='time',yaxis='rate',subplot=sub_plot,iteration='antenna',
    #        showgui=False,figfile=caltableplot, fontsize = 8, plotrange =
    #        [-1,-1,-1,-1])
    #logger.info('Fringe calibration (rates) plot in: {0}'.format(caltableplot))
    logger.info('Delay calibration plot: {0}'.format(caltableplot))
    logger.info('End delay_fringefit')
    return eMCP, caltables

def make_spwmap(caltables, combine):
    if 'spw' in combine:
        spwmap = [0]*caltables['num_spw']
    else:
        spwmap = []
    return spwmap

def make_spw(msinfo, spw_list):
    spws = spw_list[0]
    chans = spw_list[1]
    if chans == 'innerchan':
        chans = ':'+msinfo['innerchan']
    elif chans == '':
        chans = ''
    else:
        chans = ':'+chans
    return '{0}{1}'.format(spws, chans)


def gain_p_ap(eMCP, caltables, doplots=True):
    logger.info('Running gain_p_ap')
    #t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    gain_p_ap = eMCP['defaults']['initial_gaincal']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']

    # 1 Phase calibration
    caltable_name = gain_p_ap['p_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = gain_p_ap['p_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = msinfo['sources']['calsources']
    caltables[caltable_name]['solint'] = gain_p_ap['p_solint']
    caltables[caltable_name]['combine'] = gain_p_ap['p_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_p_ap['p_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['calsources']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,gain_p_ap['p_combine'])
    caltables[caltable_name]['interp'] = gain_p_ap['p_interp']
    caltables[caltable_name]['minblperant'] = gain_p_ap['p_minblperant']
    caltables[caltable_name]['minsnr'] = gain_p_ap['p_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('G0 phase {0}: {1}'.format(caltable_name,caltable))
    # Plots
    if doplots:
        caltableplot = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot, title='Phase',
                      xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=3)
        logger.info('G0 phase plot: {0}'.format(caltableplot))
    logger.info('Apply phase calibration flags to calibrators, applymode=flagonly')
    applycal(vis=msfile, gaintable=caltables[caltable_name]['table'],
             field=msinfo['sources']['calsources'],
             applymode='flagonly')

    # 2 Amplitude calibration
    caltable_name = gain_p_ap['ap_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = gain_p_ap['ap_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = msinfo['sources']['calsources']
    caltables[caltable_name]['solint'] = gain_p_ap['ap_solint']
    caltables[caltable_name]['combine'] = gain_p_ap['ap_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, gain_p_ap['ap_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['calsources']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,gain_p_ap['ap_combine'])
    caltables[caltable_name]['interp'] = gain_p_ap['ap_interp']
    caltables[caltable_name]['minblperant'] = gain_p_ap['ap_minblperant']
    caltables[caltable_name]['minsnr'] = gain_p_ap['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Gain amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    if doplots:
        caltableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot_phs, title='Phase',
                      xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=3)
        logger.info('G1_ap phase plot: {0}'.format(caltableplot_phs))
        caltableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
        plot_caltable(msinfo, caltables[caltable_name], caltableplot_amp, title='Amp',
                      xaxis='time', yaxis='amp', ymin=-1, ymax=-1, coloraxis='spw', symbolsize=5)
        logger.info('G1_ap amp plot: {0}'.format(caltableplot_amp))

    logger.info('Finished gain_p_ap')
    return eMCP, caltables



def find_anten_fluxscale(antennas):
    if 'Lo' in antennas:
        anten_for_flux = ['!Lo']
    else:
        anten_for_flux = ['']
    # Previous logic, antenna selection
    # This function tries to remove Lo and De from the fluxscale determination.
    # But only if there are enough antennas to have at least 4 of them.
    # I found unstable solutions if too few antennas are used
    #if len(antennas) > 4:
    #     anten_for_flux = [x for x in antennas if x != "Lo"]
    #     if len(anten_for_flux) > 4:
    #         anten_for_flux = [x for x in anten_for_flux if x != "De"]
    #else:
    #    anten_for_flux = antennas
    return anten_for_flux


def read_source_model(model, field, msinfo, eMcalflux):
    logger.info('Found tt0 model for source {0}: {1}'.format(field, model[0]))
    logger.info('Found tt1 model for source {0}: {1}'.format(field, model[1]))
    model_path, modelfilename = os.path.split(model[0])
    model_name = os.path.splitext(modelfilename)[0]
    scaled_model_base = '{0}/{1}_{2}.'.format(model_path,msinfo['msfilename'],model_name)
    scaled_model = [scaled_model_base+ext for ext in ['tt0','tt1']]
    rmdir(scaled_model[0])
    rmdir(scaled_model[1])
    flux_in_model = imstat(model[0])
    logger.info('Flux in model tt0 (sum): {0:5.3g}, rms: {1:5.3g}, mean: {2:5.3g}'.format(
        flux_in_model['sum'][0],
        flux_in_model['rms'][0],
        flux_in_model['mean'][0]))
    logger.info('Scaling model {0} to eMcalflux {1} Jy to create {2}'.format(
                            model[0],
                            eMcalflux,
                            scaled_model[0]))
    logger.info('Scaling model {0} to eMcalflux {1} Jy to create {2}'.format(
                            model[1],
                            eMcalflux,
                            scaled_model[1]))
    logger.info('The multiplying factor is: {0:5.3f}'.format(eMcalflux/flux_in_model['sum'][0]))
    immath(imagename= model[0],
        mode='evalexpr',
        expr = 'IM0*{}'.format(eMcalflux/flux_in_model['sum'][0]),
        outfile = scaled_model[0])
    immath(imagename= model[1],
        mode='evalexpr',
        expr = 'IM0*{}'.format(eMcalflux/flux_in_model['sum'][0]),
        outfile = scaled_model[1])
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
    cals_to_scale = sources['cals_no_fluxcal']
    fluxcal = sources['fluxcal']
    caltable_name = flux['tablename']
    # Check if table allcal_ap.G1 exists:
    check_table_exists(caltables, ampcal_table)
    caltables[caltable_name] = copy.copy(caltables[ampcal_table])
    caltables[caltable_name]['table']=caltables[ampcal_table]['table']+'_fluxscaled'
    fluxes_txt = info_dir+ caltables[caltable_name]['name']+'_fluxes.txt'
    logger.info('Flux density scale from: {0}'.format(fluxcal))
    logger.info('Transfered to: {0}'.format(cals_to_scale))
    logger.info('Input caltable: {0}'.format(caltables[ampcal_table]['table']))
    logger.info('Antennas used to scale: {0}'.format(anten_for_flux))
    calfluxes = fluxscale(vis=msfile, reference=fluxcal,
                          transfer=cals_to_scale,
                          antenna = ','.join(anten_for_flux),
                          caltable  = caltables[ampcal_table]['table'],
                          fluxtable = caltables[caltable_name]['table'],
                          listfile  = fluxes_txt)
    logger.info('Modified caltable: {0}'.format(caltables[caltable_name]['table']))
    logger.info('Spectrum information: {0}'.format(caltables[caltable_name]['table']+'_fluxes.txt'))
    save_obj(calfluxes, info_dir + 'calfluxes.pkl')
    if calfluxes == None:
        logger.critical('Something went wrong with fluxscale')
        logger.warning('This probably means that necessary data are missing:')
        logger.warning('Required sources: {0}'.format(cals_to_scale))
        logger.warning('Required antennas: {0}'.format(','.join(anten_for_flux)))
        exit_pipeline(msinfo)
    # Compute correction to scale the flux density of 3C286 according to
    # resolution provided by the shortest available baseline of e-MERLIN
    eMfactor = calc_eMfactor(msfile, field=fluxcal)
    # Include a note in the fluxes.txt file warning that the values in that
    # file should be corrected by eMfactor
    with open(fluxes_txt, 'a') as file:
        file.write('# WARNING: All flux densities in this file need to be multiplied by eMfactor={0:6.4f} to match the corrections that have been applied to the data.'.format(eMfactor))
    # Get fitted flux density and spectral index, correctly scaled for e-MERLIN
    eMcalfluxes = collections.OrderedDict()
    for k in calfluxes.keys():
        if len(calfluxes[k]) > 4:
            try:
                a=[]
                a.append(calfluxes[k]['fitFluxd']*eMfactor)
                a.append(calfluxes[k]['spidx'][1])
                a.append(calfluxes[k]['fitRefFreq'])
                eMcalfluxes[calfluxes[k]['fieldName']]=a
                logger.info('Spectrum for {0:>9s}: Flux density = {1:6.3f} +/-{2:6.3f}, spidx ={3:5.2f}+/-{4:5.2f}'.format(calfluxes[k]['fieldName'],
                    calfluxes[k]['fitFluxd']*eMfactor, calfluxes[k]['fitFluxdErr']*eMfactor,
                    calfluxes[k]['spidx'][1], calfluxes[k]['spidxerr'][1]))
            except:
                pass
    # Scale and fill model column:
    for field in eMcalfluxes.keys():
        # Check if there are image models (tt0, tt1) for each field:
        model = ['./source_models/{0}.model.tt0'.format(field),
                 './source_models/{0}.model.tt1'.format(field)]
        if os.path.exists(model[0]) and os.path.exists(model[1]):
            scaled_model = read_source_model(model, field, msinfo,
                                             eMcalfluxes[field][0])
            logger.info('New model for this observation: {0}, {1}'.format(scaled_model[0], scaled_model[1]))
            ft(vis=msfile, field=field, model=scaled_model, nterms=2, usescratch=True)
            logger.info('Model for {0} included in MODEL column in {1}'.format(field, msfile))
        else:
            # If there is no model, just assume point source:
            logger.info('Filling model column for point-like source: {0}'.format(field))
            logger.info('Model: flux={0:6.3g}Jy, spix={1:6.3g}, reffreq={2:6.3g}Hz'.format(
                                                eMcalfluxes[field][0],
                                                eMcalfluxes[field][1],
                                                eMcalfluxes[field][2]))
            setjy(vis = msfile,
                  field = field,
                  standard = 'manual',
                  fluxdensity = eMcalfluxes[field][0],
                  spix = eMcalfluxes[field][1],
                  reffreq = str(eMcalfluxes[field][2])+'Hz',
                  usescratch = True)
    logger.info('End eM_fluxscale')
    # Apply calibration if requested:
    if eMCP['inputs']['fluxscale'] == 2:
        #logger.warning('Applying short-solint tables to target.')
        run_applycal(eMCP, caltables, step = 'fluxscale')
    save_obj(caltables, caltables['calib_dir']+'caltables.pkl')
    msg = ''
    eMCP = add_step_time('fluxscale', eMCP, msg, t0)
    return eMCP, caltables



def bandpass_sp(eMCP, caltables):
    logger.info('Start bandpass_sp')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    bp_sp = eMCP['defaults']['bandpass_sp']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # Bandpass calibration
    caltable_name = bp_sp['bp_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = bp_sp['bp_prev_cal']
    caltables[caltable_name]['field'] = msinfo['sources']['bpcal']
    caltables[caltable_name]['solint'] = bp_sp['bp_solint']
    caltables[caltable_name]['spw'] = make_spw(msinfo, bp_sp['bp_spw'])
    caltables[caltable_name]['combine'] = bp_sp['bp_combine']
    caltables[caltable_name]['uvrange'] = bp_sp['bp_uvrange']
    caltables[caltable_name]['fillgaps'] = bp_sp['bp_fillgaps']
    caltables[caltable_name]['solnorm'] = bp_sp['bp_solnorm']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,bp_sp['bp_combine'])
    caltables[caltable_name]['interp'] = bp_sp['bp_interp']
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name)
    caltables[caltable_name]['gainfield'] = get_unique_field(caltables[caltable_name]['table'])
    logger.info('Bandpass_sp BP {0}: {1}'.format(caltable_name,bptable))
    # Plots
    bptableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
    bptableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plot_caltable(msinfo, caltables[caltable_name], bptableplot_phs, title='Bandpass phase',
                  xaxis='freq', yaxis='phase', ymin=-180, ymax=180, coloraxis='corr', symbolsize=5)
    logger.info('Bandpass_sp BP phase plot: {0}'.format(bptableplot_phs))
    plot_caltable(msinfo, caltables[caltable_name], bptableplot_amp, title='Bandpass amp',
                  xaxis='freq', yaxis='amp', ymin=-1, ymax=-1, coloraxis='corr', symbolsize=5)
    logger.info('Bandpass_sp BP amplitude plot: {0}'.format(bptableplot_amp))
    logger.info('End bandpass_sp')
    # Apply calibration if requested:
    if eMCP['inputs']['bandpass_sp'] == 2:
        run_applycal(eMCP, caltables, step = 'bandpass_sp')
    save_obj(caltables, caltables['calib_dir']+'caltables.pkl')
    msg = 'field={0}, combine={1}, solint={2}'.format(
                    msinfo['sources']['bpcal'],
                    bp_sp['bp_combine'],
                    bp_sp['bp_solint'])
    eMCP = add_step_time('bandpass_sp', eMCP, msg, t0)
    return eMCP, caltables


def gain_amp_sp(eMCP, caltables):
    logger.info('Start gain_amp_sp')
    t0 = datetime.datetime.utcnow()
    # Check if all sources are in the MS: 
    check_sources_in_ms(eMCP)
    amp_sp = eMCP['defaults']['gain_amp_sp']
    msfile = eMCP['msinfo']['msfile']
    msinfo = eMCP['msinfo']
    # 1 Amplitude calibration
    caltable_name = amp_sp['ap_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = amp_sp['ap_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = msinfo['sources']['calsources']
    caltables[caltable_name]['solint'] = amp_sp['ap_solint']
    caltables[caltable_name]['combine'] = amp_sp['ap_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, amp_sp['ap_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['calsources']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,amp_sp['ap_combine'])
    caltables[caltable_name]['interp'] = amp_sp['ap_interp']
    caltables[caltable_name]['minblperant'] = amp_sp['ap_minblperant']
    caltables[caltable_name]['minsnr'] = amp_sp['ap_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('Gain amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot_phs, title='Phase',
                  xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=5)
    logger.info('{0} phase plot: {1}'.format(caltable_name,caltableplot_phs))
    caltableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot_amp, title='Amp',
                  xaxis='time', yaxis='amp', ymin=-1, ymax=-1, coloraxis='spw', symbolsize=5)
    logger.info('{0} amp plot: {1}'.format(caltable_name,caltableplot_amp))

    # 2 Scan-averaged phase calibration
    caltable_name = amp_sp['p_scan_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = amp_sp['p_scan_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['field'] = msinfo['sources']['phscals']
    caltables[caltable_name]['solint'] = amp_sp['p_scan_solint']
    caltables[caltable_name]['combine'] = amp_sp['p_scan_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, amp_sp['p_scan_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['phscals']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,amp_sp['p_scan_combine'])
    caltables[caltable_name]['interp'] = amp_sp['p_scan_interp']
    caltables[caltable_name]['minblperant'] = amp_sp['p_scan_minblperant']
    caltables[caltable_name]['minsnr'] = amp_sp['p_scan_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('{0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot_phs, title='Phase',
                  xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=5)

    # 3 Amplitude calibration on phasecal: scan-averaged amplitude solutions
    caltable_name = amp_sp['ap_scan_tablename']
    caltables[caltable_name] = collections.OrderedDict()
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['previous_cal'] = amp_sp['ap_scan_prev_cal']
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['field'] = msinfo['sources']['phscals']
    caltables[caltable_name]['solint'] = amp_sp['ap_scan_solint']
    caltables[caltable_name]['combine'] = amp_sp['ap_scan_combine']
    caltables[caltable_name]['spw'] = make_spw(msinfo, amp_sp['ap_scan_spw'])
    caltables[caltable_name]['gainfield'] = msinfo['sources']['phscals']
    caltables[caltable_name]['spwmap'] = make_spwmap(caltables,amp_sp['ap_scan_combine'])
    caltables[caltable_name]['interp'] = amp_sp['ap_scan_interp']
    caltables[caltable_name]['minblperant'] = amp_sp['ap_scan_minblperant']
    caltables[caltable_name]['minsnr'] = amp_sp['ap_scan_minsnr']
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name)
    if caltables['Lo_dropout_scans'] != '':
        remove_missing_scans(caltable, caltables['Lo_dropout_scans'])
#    logger.info('{0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot_phs, title='Phase',
                  xaxis='time', yaxis='phase', ymin=-180, ymax=180, coloraxis='spw', symbolsize=5)
    logger.info('{0} phase plot: {1}'.format(caltable_name,caltableplot_phs))
    caltableplot_amp = caltables['plots_dir']+'caltables/'+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plot_caltable(msinfo, caltables[caltable_name], caltableplot_amp, title='Amp',
                  xaxis='time', yaxis='amp', ymin=-1, ymax=-1, coloraxis='spw', symbolsize=5)
    logger.info('{0} amp plot: {1}'.format(caltable_name,caltableplot_amp))
    logger.info('End gain_amp_sp')
    # Apply calibration if requested:
    if eMCP['inputs']['gain_amp_sp'] == 2:
        run_applycal(eMCP, caltables, step = 'gain_amp_sp')
    save_obj(caltables, caltables['calib_dir']+'caltables.pkl')
    msg = 'ap_solint={0}'.format(amp_sp['ap_solint'])
    eMCP = add_step_time('gain_amp_sp', eMCP, msg, t0)
    return eMCP, caltables

def applycal_all(eMCP, caltables):
    #logger.info('Start applycal_all')
    t0 = datetime.datetime.utcnow()
    run_applycal(eMCP, caltables, step='applycal_all', dotarget=True)
    #logger.info('End applycal_all')
    flag_statistics(eMCP, step='applycal_all')
    msg = ''
    eMCP = add_step_time('applycal_all', eMCP, msg, t0)
    return eMCP


def compile_statistics(msfile, tablename=''):
    logger.info('Start compile_stats')
    # Num of spw and baselines
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    baselines = get_baselines(msfile)
    # Date and time of observation
    ms.open(msfile)
    axis_info = ms.getdata(['axis_info'],ifraxis=True)
    ms.close()
    vis_field = vishead(msfile,mode='list',listitems='field')['field'][0][0] # Only one expected
    t_mjd, t = get_dates(axis_info)
    freq_ini = np.min(axis_info['axis_info']['freq_axis']['chan_freq'])/1e9
    freq_end = np.max(axis_info['axis_info']['freq_axis']['chan_freq'])/1e9
    chan_res = np.mean(axis_info['axis_info']['freq_axis']['resolution'])/1e9
    band = check_band(msfile)
    data_stats = []
    for bsl in baselines:
        print('Processing baseline: {}'.format(bsl))
        for spw in range(num_spw):
            for corr in ['RR', 'LL']:
                d = visstat(msfile, antenna=bsl, axis='amplitude',
                            spw=str(spw)+':100~400', correlation=corr)
                data_stats.append([t[0],
                                   t[-1],
                                   np.mean(t_mjd),
                                   band,
                                   freq_ini,
                                   freq_end,
                                   chan_res,
                                   vis_field,
                                   bsl,
                                   spw,
                                   corr,
                                   d['DATA']['mean'],
                                   d['DATA']['npts'],
                                   d['DATA']['median'],
                                   d['DATA']['stddev']])
                print t[0], vis_field, bsl, spw, corr, d['DATA']['median'], d['DATA']['stddev']
    data_stats = np.asarray(data_stats)
    ini_date = str(t[0]).replace('-', '').replace(':','').replace(' ','_').split('.')[0]
    outname = '{0}_{1}_{2}.npy'.format(ini_date, vis_field, band)
    np.save('data_'+outname, data_stats)
    logger.info('Data statistics saved to: {0}'.format(outname))
    if tablename != '':
       compile_delays(tablename, outname)
    logger.info('End compile_stats')


def compile_delays(tablename, outname):
    tb.open(tablename+'/ANTENNA')
    antennas = tb.getcol('NAME')
    tb.close()
    tb.open(tablename)
    a = tb.getcol('ANTENNA1')
    times  = tb.getcol('TIME')
    delays = tb.getcol('FPARAM')
    tb.close()
    delay_stats = []
    for i in range(len(times)):
        delay_stats.append([antennas[a[i]], 'RR', delays[0,0,i]])
        delay_stats.append([antennas[a[i]], 'LL', delays[1,0,i]])
    delay_stats = np.asarray(delay_stats)
    np.save('delay_'+outname, delay_stats)
    logger.info('Delay statistics saved to: {0}'.format(outname))


#def monitoring(msfile, msinfo, caltables, previous_cal):
#    # This is intented to run on a single-source file for daily monitoring on
#    # unaveraged data
#    logger.info('Starting monitoring')
#    band = check_band(msfile)
#    if band == 'L':
#        hanning(inputvis=msfile,deloriginal=True)
#
#    flagdata1_apriori(msfile=msfile, msinfo=msinfo, do_quack=True)
#
#    flagdata_tfcrop_bright(msfile=msfile, sources=msinfo['sources'])
#    caltables = solve_delays(msfile=msfile, caltables=caltables,
#                   previous_cal=[], calsources=msinfo['sources']['calsources'],
#                             solint='60s')
#    run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
#                    previous_cal=['delay.K1'],
#                    previous_cal_targets=['delay.K1'])
#    compile_statistics(msfile, tablename=caltables['delay.K1']['table'])
#    return caltables


def calc_eMfactor(msfile, field='1331+305'):
    logger.info('Computing eMfactor')
    if field != '1331+305':
        logger.warning('Scaling flux assuming 3C286 is the flux calibrator. Your flux calibrator is: {}. Scaling is probably wrong.'.format(field))
    tb.open(msfile+'/FIELD')
    names = tb.getcol('NAME')
    field_id = np.argwhere(names == field)[0][0]
    tb.close()

    tb.open(msfile+'/ANTENNA')
    anten = tb.getcol('NAME')
    tb.close()

    tb.open(msfile+'/ANTENNA')
    anten = tb.getcol('NAME')
    tb.close()
    try:
        Lo_id = np.argwhere(anten=='Lo')[0][0]
    except:
        Lo_id = -1

    tb.open(msfile)
    uvw = tb.getcol('UVW')
    a1 = tb.getcol('ANTENNA1')
    a2 = tb.getcol('ANTENNA2')
    field = tb.getcol('FIELD_ID')
    tb.close()

    uvdist = np.sqrt(uvw[0]**2+uvw[1]**2)
    mask = (uvdist==0) + (field!=field_id) + (anten[a1]=='Lo')
    uvdist_nonzero = np.ma.array(uvdist, mask=mask)

    # To exclude completely flagged data:
    # I comment this out because for some data sets tb.getcol returns a flatten
    # array without shape
    #flags = tb.getcol('FLAG')
    #allflag = np.sum(flags, axis=(0,1))
    #flags_entries = flags.shape[0]*flags.shape[1]
    #mask = (uvdist==0) | (allflag==flags_entries) | field!=field_id

    n = np.argmin(uvdist_nonzero)

    tb.open(msfile+'/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    tb.close()

    shortest_baseline = uvdist_nonzero[n] # Shortest baseline in m
    center_freq = np.mean(chan_freq)/1e6 # Center frequency in MHz
    dfluxpy_output = dfluxpy(center_freq, shortest_baseline)
    eMfactor = dfluxpy_output[1]/dfluxpy_output[0]
    logger.info('Shortest projected baseline: {0} [{1}] - {2} [{3}] {4:10.2f}m'.format(
        anten[a1[n]], a1[n], anten[a2[n]], a2[n], uvdist_nonzero[n]))
    logger.info('Central frequency of the MS: {0} MHz'.format(center_freq))
    logger.info('eMfactor: {0:6.4f}'.format(eMfactor))
    return eMfactor


def dfluxpy(freq,baseline):
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

    # Reworked to use the 1999 VLA flux formula, and a 2nd formula to give a continuous estimate of the resolved fraction, by Ian Stewart, JBO, 8 Aug 2007.
    # Minor changes by amsr, 8 Aug 2007

    lowest_freq = 300.0;
    highest_freq = 50000.0;
    if (freq < lowest_freq or freq > highest_freq):
        print "Frequency must be between $lowest_freq and $highest_freq MHz. \n"

    # Perley & Butler 2012 values
    A = 1.2515
    B = -0.4605
    C = -0.1715
    D = 0.0336

    log10f = (math.log(freq)/2.3025851) - 3.0; # Why the -3? Because freq has to be GHz for the formula to work.
    log_flux = A + B*log10f + C*log10f*log10f + D*log10f*log10f*log10f
    vlaflux = math.pow(10.0,log_flux)

    # The VLA flux must now be corrected to account for the higher resolving power of merlin. The formula used was obtained with the help of Peter Thomasson. If we assume that 3C286 is represented by a gaussian of angular size theta_s, and represent the resolving power of the given baseline as a function of frequency f and antenna separation L by theta_b(f,L), then the reduction in central flux density A(0) due to the finite theta_s is given by
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

    ref_bl_length = 11236.79 # MK-TA separation in metres.
    ref_freq = 5000.0
    ref_rho = 0.04
    thisbl = "this baseline (Mk-Ta)"

    bl_length = baseline

    frac = (freq / ref_freq) * (bl_length / ref_bl_length)
    rho = frac * frac * ref_rho
    merlinflux = vlaflux / (1.0 + rho)

    # Another useful quantity is the resolved percentage:
    resolved_percent = 100.0 * rho / (1.0 + rho)
    caution_res_pc = 10.0

    return vlaflux, merlinflux, resolved_percent, caution_res_pc, thisbl


def plot_image(eMCP, imagename, ext='.tt0', dozoom=False):
    imstat_residual = imstat(imagename+'.residual'+ext)
    imstat_image = imstat(imagename+'.image'+ext)
    noise = imstat_residual['rms'][0]
    peak = imstat_image['max'][0]
    scaling =  np.min([0, -np.log(1.0*peak/noise)+4])
    #scaling = np.min([0.0, -int(np.log(1.0*peak/noise)+1)])
    imgpar = eMCP['defaults']['first_images']
    level0 = imgpar['level0']
    levels = level0*np.sqrt(3**np.arange(20))
    center = 512
    zoom_range = 150
    for extension in ['image'+ext]:
        filename = '{0}.{1}'.format(imagename, extension)
        logger.info('Creating png for image: {0}'.format(filename))
        logger.info('Peak: {0:5.2e} mJy, noise: {1:5.2e} mJy, ' \
                    'scaling: {2:3.1f}'.format(peak*1000.,
                                               noise*1000.,
                                               scaling))
        imview(raster={'file':filename,
                       'scaling': float(scaling),
                       'colorwedge':True},
               contour = {'file':filename,
                   'levels':list(levels),
                   'base':0,
                   'unit':float(noise)*2.},
               out = filename+'.png')
        if dozoom:
            imview(raster={'file':filename,
                           'scaling': float(scaling),
                           'colorwedge':True},
                   contour = {'file':filename,
                              'levels':list(levels),
                              'base':0,
                              'unit':float(noise)*2.},
                   zoom = {'blc':[center-zoom_range,center-zoom_range],'trc':[center+zoom_range, center+zoom_range]},
                   out = filename+'_zoom.png')
    return peak, noise

def plot_image_add(imagename, ext='.tt0', dozoom=False):
    filename = imagename
    center = 512
    zoom_range = 150
    imview(raster={'file':filename+'.residual'+ext,
                   'colorwedge':True},
           contour = {'file':filename+'.mask',
               'levels':[1]},
           out = filename+'.residual'+ext+'.png')
    if dozoom:
        imview(raster={'file':filename+'.residual'+ext,
               'colorwedge':True},
        contour = {'file':filename+'.mask',
                   'levels':[1]},
        zoom = {'blc':[center-zoom_range,center-zoom_range],'trc':[center+zoom_range, center+zoom_range]},
        out = filename+'.residual'+ext+'_zoom.png')
    return

def single_tclean(eMCP, s, num):
    msinfo = eMCP['msinfo']
    logger.info('Producing tclean images for {0}, field: {1}'.format(msinfo['msfile'], s))
    makedir(images_dir+'{}'.format(s))
    imagename = images_dir+'{0}/{1}_{0}_img{2:02d}'.format(s, msinfo['msfilename'], num)
    prev_images = glob.glob(imagename+'*')
    for prev_image in prev_images:
        rmdir(prev_image)
    cellsize = {'C':'0.008arcsec', 'L':'0.02arcsec'}
    imgpar = eMCP['defaults']['first_images']
    imsize = imgpar['imsize']
    cell = cellsize[msinfo['band']]
    niter = imgpar['niter']
    gain = imgpar['gain']
    deconvolver = imgpar['deconvolver']
    nterms = imgpar['nterms']
    weighting = imgpar['weighting']
    robust = imgpar['robust']
    usemask = 'auto-multithresh'
    nsigma = imgpar['nsigma']
    sidelobethreshold = imgpar['sidelobethreshold']
    noisethreshold = imgpar['noisethreshold']
    lownoisethreshold = imgpar['lownoisethreshold']
    minbeamfrac = imgpar['minbeamfrac']
    growiterations = imgpar['growiterations']
    logger.info('imsize = {0}, cell = {1}, niter = {2}'.format(
                imsize, cell, niter))
    logger.info('weighting = {0}, robust = {1}'.format(
                weighting, robust))
    logger.info('usemask = {0}, growiterations = {1}'.format(usemask,
                                                             growiterations))
    tclean(vis=msinfo['msfile'], field=s, datacolumn='corrected',
           imagename=imagename,
           imsize=imsize, cell=cell, deconvolver=deconvolver,
           gain=gain,
           nterms=nterms,
           weighting=weighting,
           robust=robust, niter=niter, usemask=usemask,
           nsigma=nsigma,
           sidelobethreshold=sidelobethreshold,
           noisethreshold=noisethreshold,
           lownoisethreshold=lownoisethreshold,
           minbeamfrac=minbeamfrac,
           growiterations=growiterations,
           savemodel='none', parallel=False)
    if nterms > 1:
        ext = '.tt0'
    else:
        ext = ''
    peak, noise = plot_image(eMCP, imagename, ext=ext, dozoom=True)
    plot_image_add(imagename, ext=ext, dozoom=True)
    eMCP['img_stats'][s] = [peak, noise]
    return eMCP


def run_first_images(eMCP):
    msinfo = eMCP['msinfo']
    logger.info(line0)
    logger.info('Start run_first_images')
    t0 = datetime.datetime.utcnow()
    if eMCP['defaults']['first_images']['run_statwt']:
        logger.info('Running statwt with default parameters')
        statwt(vis=msinfo['msfile'])
    num = 0
    eMCP['img_stats'] = collections.OrderedDict()
    for s in msinfo['sources']['targets_phscals'].split(','):
        eMCP = single_tclean(eMCP, s, num)
    logger.info('End run_first_images')
    msg = ''
    eMCP = add_step_time('first_images', eMCP, msg, t0)
    return eMCP

def shift_field_position(msfile, shift):
    field = shift['field']
    new_pos = shift['new_position']
    position_name = shift['new_field_name']
    logger.info('Field {0} will be shifted to {1} on {2}'.format(field, position_name, new_pos))
    msfile_split = '{0}_{1}'.format(msfile, position_name)
    mssources = vishead(msfile,mode='list',listitems='field')['field'][0]
    if field not in mssources:
        logger.critical('Requested field to shift: {} not in MS! Closing '.format(field))
        exit_pipeline(msinfo)
    rmdir(msfile_split)
    # Split
    logger.info('Splitting field: {}'.format(field))
    mstransform(msfile, outputvis=msfile_split, field=field, datacolumn='data')
    #FIXVIS
    logger.info('Changing phase center to: {}'.format(new_pos))
    fixvis(vis=msfile_split, field=field, outputvis='', phasecenter=new_pos, datacolumn='data')
    # Change field name
    tb.open(msfile_split+'/FIELD',nomodify=False)
    st=tb.selectrows(0)
    st.putcol('NAME', '{0}'.format(position_name))
    st.done()
    tb.close()
    # Concatenate again
    logger.info('Concatenating {0} into {1}'.format(msfile_split, msfile))
    concat(vis= msfile_split, concatvis= msfile)
    logger.info('Updating listobs for MS: {}'.format(msfile))
    run_listobs(msfile)
    rmdir(msfile_split)
    logger.warning('New field: {} in MS. Make sure you include it in the inputs file'.format(position_name))

def read_shifts_file(shifts_file):
    shifts_list = []
    with open(shifts_file, 'rb') as shifts_txt:
        lines = shifts_txt.readlines()
        for line in lines:
            if line.strip() != '':
                shift = {'field': line.split(',')[0].strip(),
                         'new_field_name': line.split(',')[1].strip(),
                         'new_position': line.split(',')[2].strip()}
                shifts_list.append(shift)
    return shifts_list


def shift_all_positions(eMCP):
    msfile = eMCP['msinfo']['msfile']
    logger.info('Start shift_all_pos')
    t0 = datetime.datetime.utcnow()
    shifts_file = './shift_phasecenter.txt'
    try:
        shifts_list = read_shifts_file(shifts_file)
        logger.info('Reading shifts from {0}'.format(shifts_file))
    except:
        logger.critical('Unable to open {0} with new position information'.format(shifts_file))
        exit_pipeline(msinfo)
    listobs_file = info_dir + msfile
    mvdir(listobs_file + '.listobs.txt', listobs_file+'preshift_listobs.txt')
    logger.info('Found {0} shifts to apply. {0} new fields will be added'.format(len(shifts_list)))
    for shift in shifts_list:
        shift_field_position(msfile, shift)
    run_listobs(msfile)
    logger.info('Listobs file in: {0}'.format(msfile+'.listobs.txt'))
    logger.info('End shift_all_pos')
    msg = 'file={0}'.format(shifts_file)
    eMCP = add_step_time('shift_field_pos', eMCP, msg, t0)
    return eMCP


def eMCP_info_start_steps():
    default_value = [0, 0, '']
    steps = collections.OrderedDict()
    steps['start_pipeline'] = default_value
    steps['importfitsIDI'] = default_value
    steps['mstransform'] = default_value
    steps['fixvis'] = default_value
    steps['aoflagger'] = default_value
    steps['flag_apriori'] = default_value
    steps['flag_manual'] = default_value
    steps['shift_field_pos'] = default_value
    steps['average'] = default_value
    steps['plot_data'] = default_value
    steps['save_flags'] = default_value
    steps['restore_flags'] = default_value
    steps['flag_manual_avg'] = default_value
    steps['init_models'] = default_value
    steps['bandpass'] = default_value
    steps['initial_gaincal'] = default_value
##    steps['flag_tfcropBP'] = default_value
##    steps['delay'] = default_value
##    steps['gain_p_ap'] = default_value
    steps['fluxscale'] = default_value
    steps['bandpass_sp'] = default_value
    steps['gain_amp_sp'] = default_value
    steps['applycal_all'] = default_value
    steps['flag_target'] = default_value
    steps['plot_corrected'] = default_value
    steps['first_images'] = default_value
    return steps


#########  Search Lo dropout scans  #########

def find_fields_scans(msfile):
    msmd.open(msfile)
    scans = msmd.scannumbers()
    dict_scans = msmd.fieldsforscans(scans, True, asmap=True, obsid=0, arrayid=0)
    msmd.done()
    field_for_scan = np.array([dict_scans[str(scan)][0] for scan in scans])
    return scans, field_for_scan

def find_Lo_amp_spw(msfile, phscal, phscal_scans, spw, eMCP):
    amp_mean = np.ones_like(phscal_scans) * np.nan
    amp_std = np.ones_like(phscal_scans) * np.nan
    Lo_defaults = eMCP['defaults']['flag_apriori']
    results_tmp = visstat(msfile, scan='',
                          field = phscal,
                          antenna='Lo&*',
                          spw = spw + ':'+eMCP['msinfo']['innerchan'],
                          datacolumn = Lo_defaults['Lo_datacolumn'],
                          useflags = Lo_defaults['Lo_useflags'],
                          correlation = 'RR,LL',
                          timeaverage = True,
                          timebin = '9999999s',
                          timespan = '')
    for key in results_tmp.keys():
        scan = key.split(',')[1].split('=')[-1]
        scan_idx = np.where(phscal_scans==int(scan))[0][0]
        amp_mean[scan_idx] = results_tmp[key]['median']
        amp_std[scan_idx] = results_tmp[key]['stddev']
    return  amp_mean, amp_std

def find_Lo_amp(msfile, phscal, phscal_scans, eMCP, spws):
    amp_means = np.zeros((len(spws), len(phscal_scans)))
    amp_stds = np.zeros((len(spws), len(phscal_scans)))
    logger.info('Analysing {0} scans for phscal: {1}'.format(len(phscal_scans),
                                                       phscal))
    for i, spw in enumerate(spws):
        logger.info('Processing spw: {}'.format(spw))
        amp_means_i, amp_std_i = find_Lo_amp_spw(msfile, phscal, phscal_scans,
                                                 spw, eMCP)
        amp_means[i] = amp_means_i
        amp_stds[i] = amp_std_i
    amp_mean = np.average(amp_means, weights=1.0/amp_stds, axis=0)
    amp_std = np.average(amp_stds, axis=0)
    return amp_mean, amp_std

def plot_Lo_drops(phscal_scans, amp_mean, lo_dropout_scans, phscal, eMCP):
    msinfo = eMCP['msinfo']
    drops = np.array([scan in lo_dropout_scans for scan in phscal_scans])
    fig = plt.figure(figsize=(30,8))
    ax1 = fig.add_subplot(111)

    scans, field_for_scan = find_scans(msfile)
    ax1.bar(scans-0.5, np.ones_like(scans)*np.max(amp_mean)*1.2, alpha=0.2,
            color='0.5', width=1),

    ax1.bar(phscal_scans-0.5, amp_mean, alpha=1.0,
            color='0.5', width=1,
            label='{0}'.format(phscal))
    ax1.bar(phscal_scans[drops]-0.5, amp_mean[drops], alpha=1.0,
            color='r', width=1,
            label='{0} Lo dropouts'.format(phscal))

    ax1.legend(loc=0)
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xlim(np.min(phscal_scans)-0.5, np.max(phscal_scans)+0.5)
    ax1.set_ylim(0, np.max(amp_mean)*1.2)
    ax1.set_xlabel('Scan number')
    ax1.set_ylabel('Mean spw Lo raw amplitude')

    plots_obs_dir = './weblog/plots/plots_flagstats/'
    plot_file_Lo = plots_obs_dir+'{0}_Lo_dropout_scans{1}.png'.format(msinfo['msfilename'],
                                                       phscal)
    fig.savefig(plot_file_Lo, bbox_inches='tight')

def calc_Lo_drops(amp_mean, phscal_scans, threshold=0.5):
    # Sort amplitude values:
    a = np.sort(amp_mean[~np.isnan(amp_mean)])
    # Iterate to create two groups. Compute sum of std of both groups
    # That sum will be minimum when the two groups are divided by the 
    # amplitude intersecting the two real distributions
    astd = np.array([0.0])
    for i in range(1,len(a)):
        astd=np.append(astd,np.std(a[:i])+np.std(a[i:]))

    if astd[1:].min() > threshold*astd[1:].max():
        logger.info('No evidence of bimodality. Lo_drop_scans will be empty')
        lo_dropout_scans = []
    else:
        astd[0]=astd.max()
        idxmin = np.argmin(astd)
        separation_amp = a[idxmin]
        drops = amp_mean<separation_amp
        lo_dropout_scans = phscal_scans[drops]
    return lo_dropout_scans

def find_Lo_drops(msfile, phscals, eMCP):
    spws = eMCP['defaults']['flag_apriori']['Lo_spws']
    min_scans = eMCP['defaults']['flag_apriori']['Lo_min_scans']
    logger.info('Searching for possible Lo dropout scans on phscal scans')
    lo_dropout_scans = np.array([], dtype='int')
    scans, field_for_scan = find_fields_scans(msfile)
    for phscal in phscals:
        logger.info('Now searching for phasecal: {}'.format(phscal))
        phscal_scans = scans[field_for_scan == phscal]
        amp_mean, amp_std = find_Lo_amp(msfile, phscal, phscal_scans, eMCP,
                                        spws)
        threshold = eMCP['defaults']['flag_apriori']['Lo_threshold']
        lo_dropout_scans_i = calc_Lo_drops(amp_mean, phscal_scans, threshold)
        emplt.plot_Lo_drops(phscal_scans, scans, amp_mean,
                            lo_dropout_scans_i, phscal, eMCP)
        logger.info('Potential dropout scans: '
                    '{0}'.format(','.join(lo_dropout_scans_i.astype('str'))))
        if len(lo_dropout_scans_i) >= min_scans:
            lo_dropout_scans = np.hstack([lo_dropout_scans, lo_dropout_scans_i])
        else:
            logger.info('Less than {} dropout scans, not considered '
                        'persistent drops'.format(min_scans))
            lo_dropout_scans = []
    return lo_dropout_scans


def flag_statistics(eMCP, step):
    msinfo = eMCP['msinfo']
    logger.info(line0)
    logger.info('Start flagstatistics')
    plots_obs_dir = './weblog/plots/plots_flagstats/'
    makedir(plots_obs_dir)
    logger.info('Running flagdata on {0}'.format(step))
    logger.info('mode="summary", action="calculate", antenna="*&*"'.format(msinfo['msfile']))
    flag_stats = flagdata(vis=msinfo['msfile'], mode='summary', action='calculate', display='none', antenna='*&*')
    outfile = weblog_dir + 'plots/plots_flagstats/flagstats_{}.pkl'.format(step)
    save_obj(flag_stats, outfile)
    logger.info('flagstats file saved to: {}'.format(outfile))
    logger.info('Flag statistics ready. Now plotting.')
    emplt.plot_flagstatistics(flag_stats, msinfo, step)
    logger.info('End flagstatistics')
