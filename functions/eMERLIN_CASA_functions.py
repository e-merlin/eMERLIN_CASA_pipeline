#!/usr/local/python
import os
import numpy as np
import pickle
from Tkinter import *
import tkMessageBox
import sys, shutil
import copy
import getopt
import datetime
from eMERLIN_CASA_GUI import GUI_pipeline

# CASA imports
from taskinit import *
from tasks import *
from recipes.setOrder import setToCasaOrder

import logging
logger = logging.getLogger('logger')

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
            inputs = headless(a) ## read input file
            inputs['quit'] = 0 ##needed to add to be compatible with GUI
            logger.info('Inputs from file: {}'.format(a))
            logger.info('Inputs used: {}'.format(inputs))
        elif o in ('-g','--gui'):
            inputs = GUI_pipeline(pipeline_path).confirm_parameters() ## read input file
            logger.info('inputs from GUI: {}'.format(inputs))
        elif o in ('-h','--help'):
            logger.debug('help will be written soon')
            sys.exit()
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
    INPUTFILE = open(inputfile, "r")
    control = {}
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
                    control[param] = int(valuelist[0])
                else:
                    control[param] = str(valuelist[0])
            else:
                control[param] = ','.join(valuelist)
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



# Functions to save and load dictionaries
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


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
    ms.open(msfile)
    freq = ms.getdata(['axis_info'])['axis_info']['freq_axis']['chan_freq'][0][0]/1e9
    ms.close()
    band = ''
    if (freq > 1.2) and (freq < 1.7):
        band = 'L'
    if (freq > 4) and (freq < 8):
        band = 'C'
    if (freq > 22) and (freq < 24):
        band = 'K'
    return band

def get_baselines(msfile):
    ms.open(msfile)
    baselines0 = ms.getdata(['axis_info'],ifraxis=True)['axis_info']['ifr_axis']['ifr_name']
    ms.close()
    baselines = []
    for bsl_name in baselines0:
        ant =  bsl_name.split('-')
        bsl = bsl_name.replace('-', '&')
        if ant[0] != ant[1]:
            baselines.append(bsl)
    return baselines

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
    sources = {}
    sources['targets'] = inputs['targets']
    sources['phscals'] = inputs['phscals']
    sources['fluxcal'] = inputs['fluxcal']
    sources['bpcal']   = inputs['bpcal']
    sources['ptcal']   = inputs['ptcal']
    sources['calsources'] = join_lists([sources['phscals'], sources['fluxcal'],sources['bpcal'], sources['ptcal']])
    sources['maincal'] = join_lists([sources['fluxcal'],sources['bpcal'], sources['ptcal']])
    sources['allsources'] = join_lists([sources['calsources'], sources['targets']])
    sources['no_fluxcal'] = sources['allsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
    sources['cals_no_fluxcal'] = sources['calsources'].replace(sources['fluxcal'], '').replace(',,',',').strip(',')
    sources['targets_phscals'] = join_lists([sources['targets'],sources['phscals']])
    logger.info('Targets:   {0}'.format(sources['targets']))
    logger.info('Phasecals: {0}'.format(sources['phscals']))
    logger.info('Fluxcal:   {0}'.format(sources['fluxcal']))
    logger.info('Bandpass:  {0}'.format(sources['bpcal']))
    logger.info('Pointcal:  {0}'.format(sources['ptcal']))
    return sources

def get_antennas(msfile):
    # Antenna list 
    ms.open(msfile)
    d = ms.getdata(['axis_info'],ifraxis=True)
    ms.close()
    antennas = np.unique('-'.join(d['axis_info']['ifr_axis']['ifr_name']).split('-'))
    logger.info('Antennas in MS {0}: {1}'.format(msfile, antennas))
    return antennas

def prt_dict(d):
    for key in d.keys():
        if type(d[key]) == dict:
            prt_dict(d[key])
        else:
            print('{0:20s} {1}'.format(key, d[key]))

def get_timefreq(msfile):
    # Date and time of observation
    ms.open(msfile)
    axis_info = ms.getdata(['axis_info'], ifraxis=True)
    ms.close()
    t_mjd, t = get_dates(axis_info)
    freq_ini = np.min(axis_info['axis_info']['freq_axis']['chan_freq'])/1e9
    freq_end = np.max(axis_info['axis_info']['freq_axis']['chan_freq'])/1e9
    chan_res = np.mean(axis_info['axis_info']['freq_axis']['resolution'])/1e9
    nchan = len(axis_info['axis_info']['freq_axis']['chan_freq'][:,0])
    return t[0], t[-1], freq_ini, freq_end, chan_res, nchan

def ms_sources(msfile):
    mssources = ','.join(vishead(msfile,mode='list',listitems='field')['field'][0])
    logger.info('Sources in MS {0}: {1}'.format(msfile, mssources))
    return mssources

def get_msinfo(msfile, inputs, doprint=False):
    logger.info('Reading ms file information for MS: {0}'.format(msfile))
    msinfo = {}
    msinfo['msfile'] = msfile
    msinfo['sources'] = user_sources(inputs)
    msinfo['mssources'] = ms_sources(msfile)
    msinfo['antennas'] = get_antennas(msfile)
    msinfo['band'] = check_band(msfile)
    msinfo['baselines'] = get_baselines(msfile)
    msinfo['num_spw'] = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    t_ini, t_end, freq_ini, freq_end, chan_res, nchan = get_timefreq(msfile)
    msinfo['t_ini'] = t_ini
    msinfo['t_end'] = t_end
    msinfo['freq_ini'] = freq_ini
    msinfo['freq_end'] = freq_end
    msinfo['chan_res'] = chan_res
    msinfo['nchan'] = nchan
    save_obj(msinfo, msfile)
    logger.info('Saving information of MS {0} in: {1}'.format(msfile, msfile+'.pkl'))
    if doprint:
        prt_dict(msinfo)
    return msinfo


def run_importfitsIDI(data_dir,vis, setorder=False):
    logger.info('Starting importfitsIDI procedure')
    rmdir(vis)
    fitsfiles =[]
    for file in os.listdir(data_dir):
        if file.endswith('fits') or file.endswith('FITS'):
            fitsfiles = fitsfiles + [data_dir+file]
            logger.info('FITS file found to be imported: {0}'.format(file))
    logger.info('Start importfitsIDI')
    importfitsidi(fitsidifile=fitsfiles, vis=vis, constobsid=True, scanreindexgap_s=15.0)
    ms.writehistory(message='eMER_CASA_Pipeline: Import fitsidi to ms, complete',msname=vis)
    if setorder:
        logger.info('Setting MS order with setToCasaOrder')
        mvdir(vis, vis+'_noorder')
        setToCasaOrder(inputMS=vis+'_noorder', outputMS=vis)
        rmdir(vis+'_noorder')
    ms.writehistory(message='eMER_CASA_Pipeline: setToCasaOrder, complete',msname=vis)
    logger.info('End importfitsIDI')
    logger.info('Start UVFIX')
    fixvis(vis=vis,outputvis=vis+'.uvfix',reuse=False)
    logger.info('End UVFIX')
    rmdir(vis)
    mvdir(vis+'.uvfix', vis)
    logger.info('Start flagdata_autocorr')
    flagdata(vis=vis,mode='manual',autocorr=True)
    ms.writehistory(message='eMER_CASA_Pipeline: Fixed uv coordinates & remove autocorr',msname=vis)
    logger.info('End flagdata_autocorr')
    logger.debug('You have been transformed from an ugly UVFITS to beautiful MS')
    listobs(vis=vis, listfile=vis+'.listobs')
    logger.info('Listobs file in: {0}'.format(vis+'.listobs'))
    return

##Hanning smoothing and flag of autocorrelations, will delete original and rename
def hanning(inputvis,deloriginal):
    logger.info('Start hanning')
    if inputvis[-3:].lower() == '.ms':
        outputvis = inputvis[:-3]+'_hanning'+inputvis[-3:]
        rmdir(outputvis)
        hanningsmooth(vis=inputvis,outputvis=outputvis,datacolumn='data')
    elif inputvis[-3:].lower() == 'mms':
        outputvis = inputvis[:-4]+'_hanning'+inputvis[-4:]
        rmdir(outputvis)
        mstransform(vis=inputvis,outputvis=outputvis,hanning=True,datacolumn='data')
    if deloriginal==True:
        rmdir(inputvis)
        rmdir(inputvis+'.flagversions')
    else:
        mvdir(inputvis, inputvis+'_prehanning')
        mvdir(inputvis+'.flagversions', inputvis+'_prehanning.flagversions')
    mvdir(outputvis, inputvis)
    mvdir(outputvis+'.flagversions', inputvis+'.flagversions')
    ms.writehistory(message='eMER_CASA_Pipeline: Hanning smoothed data, complete',msname=inputvis)
    logger.info('End hanning')
    return

def run_rfigui(vis):
    logger.info('Start run_rfigui')
    """This function should output the new strategies to /aoflagger_strategies/user/<field>.rfis in
    either the local folder or the pipeline folder."""
    os.system('rfigui '+vis)
    logger.info('End run_rfigui')


def run_aoflagger_fields(vis, flags, fields='all', pipeline_path='./'):
    """This version of the autoflagger iterates through the ms within the mms structure selecting individual fields. It uses pre-defined strategies. The np.unique in the loop below is needed for single source files. The mms thinks there are many filds (one per mms). I think it is a bug from virtualconcat."""
    logger.info('Start run_aoflagger_fields')
    vis_fields = vishead(vis,mode='list',listitems='field')['field'][0]
    fields_num = {f:i for i,f in enumerate(vis_fields)}
    if fields == 'all':
        fields = vis_fields
    else:
        fields = np.atleast_1d(fields)
    old_aoflagger = check_aoflagger_version()
    for field in np.unique(fields):
        # First, check if user has a new strategy for this field in the local folder.
        # If not, check if user has produced a new strategy for this field in the pipeline folder (for typical sources, etc).
        # If not, check for default strategies for this field
        # If nothing is found, just use the default strategy
        if os.path.isfile('./aoflagger_strategies/user/{0}.rfis'.format(field))==True:
            aostrategy = './aoflagger_strategies/user/{0}.rfis'.format(field)
        elif os.path.isfile(pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field))==True:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field)
        elif os.path.isfile(pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field))==True:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format(field)
        else:
            aostrategy = pipeline_path+'aoflagger_strategies/default/{0}.rfis'.format('default_faint')
        logger.info('Running AOFLagger for field {0} ({1}) using strategy {2}'.format(field,fields_num[field], aostrategy))
        if old_aoflagger: # < 2.9
            flagcommand = 'time aoflagger -strategy {0} {1}'.format(aostrategy, vis+'/SUBMSS/*{0}.mms.*.ms'.format(field))
        else: # >= 2.9
            flagcommand = 'time aoflagger -fields {2} -strategy {0} {1}'.format(aostrategy, vis, fields_num[field])
        os.system(flagcommand+' | tee -a pre-cal_flag_stats.txt')
        ms.writehistory(message='eMER_CASA_Pipeline: AOFlag field {0} with strategy {1}:'.format(field, aostrategy),msname=vis)
    flag_applied(flags, 'flagdata0_aoflagger')
    logger.info('End run_aoflagger_fields')
    return flags

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

def ms2mms(vis,mode):
    logger.info('Start ms2mms')
    if mode == 'parallel':
        partition(vis=vis,outputvis=vis[:-3]+'.mms',createmms=True,separationaxis="baseline",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)
        if os.path.isdir(vis[:-3]+'.mms') == True:
            rmdir(vis)
            rmdir(vis+'.flagversions')
        ms.writehistory(message='eMER_CASA_Pipeline: Converted MS to MMS for parallelisation',msname=vis[:-3]+'.mms')

    ## Need to use single if you need to aoflag the data later
    if mode == 'single':
        partition(vis=vis,outputvis=vis[:-3]+'.ms',createmms=False,separationaxis="auto",numsubms="auto",flagbackup=True,datacolumn=
"all",field="",spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",disableparallel=None,ddistart=None
,taql=None)
        if os.path.isdir(vis[:-3]+'.ms') == True:
           rmdir(vis)
           rmdir(vis+'.flagversions')
    logger.info('End ms2mms')

def ms2mms_fields(msfile):
    logger.info('Start ms2mms_fields')
    output_mmsfile = msfile[:-3]+'.mms'
    fields = vishead(msfile, mode = 'list', listitems = 'field')['field'][0]
    mmsfiles = []
    for field in fields:
        logger.info('Running partition on field found: {}'.format(field))
        mmsfile = msfile[:-3]+'_'+field+'.mms'
        mmsfiles.append(mmsfile)
        partition(vis=msfile, outputvis=mmsfile, createmms=True, separationaxis="baseline", numsubms="auto", flagbackup=False, datacolumn="all", field= field, spw="", scan="", antenna="", correlation="", timerange="", intent="", array="", uvrange="", observation="", feed="", disableparallel=None, ddistart=None, taql=None)
    # Virtual concatenation. No data copied, just moved to SUBMMS directory
    if len(mmsfiles) == 1: # No need to concatenate because there is only one file
        os.system('mv {0} {1}'.format(mmsfiles[0], output_mmsfile))
    elif len(mmsfiles) > 1:
        virtualconcat(vis = mmsfiles, concatvis = output_mmsfile, copypointing=True)
    else:
        logger.critical('No MMS files found.')
    if os.path.isdir(msfile) == True:
        rmdir(msfile)
        rmdir(msfile+'.flagversions')
    logger.info('End ms2mms_fields')

def do_prediagnostics(vis,plot_dir):
    ##Pre diagnostics for measurement sets##
    ## Includes:
    ## - Antenna positions
    ## - Amplitude vs. Time
    ## - Amplitude vs. Frequency
    ## - Phase vs. Time
    ## - Phase vs. Frequency
    ## - Closures (if task is available)
    ## - Listobs summary
    logger.info('Start prediagnostics')
    if os.path.isdir(plot_dir) == False:
        os.system('mkdir '+plot_dir)
    if os.path.isdir('./'+plot_dir+'pre-calibration') == False:
        os.system('mkdir ./'+plot_dir+'pre-calibration')
    directory = plot_dir
    ## Get information from ms
    x = vishead(vis,mode='list',listitems='field')['field'][0]
    tb.open(vis+'/SPECTRAL_WINDOW')
    nChan = str(tb.getcol('NUM_CHAN')[0])
    time = str(10E6)
    for i in range(len(x)):
        ## - Amplitude vs. Time
        plotms(vis=vis,xaxis='time',yaxis='amplitude',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel=str(nChan), iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time.pdf',highres=True ,dpi=1200,expformat='pdf',exprange='all', showgui=False)
        os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time_* '+directory+'Pre-cal_amp_vs_time_'+vis+'_'+x[i]+'.pdf')
        os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_time_*')
        ## - Amplitude vs frequency
        plotms(vis=vis,xaxis='frequency',yaxis='amplitude',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel='1', avgtime=time, iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency.pdf',highres=True ,dpi=1200,expformat='pdf',exprange='all', showgui=False)
        os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency_* '+directory+'Pre-cal_amp_vs_frequency_'+vis+'_'+x[i]+'.pdf')
        os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_amp_vs_frequency_*')

        ## - Phase vs time
        plotms(vis=vis,xaxis='time',yaxis='phase',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgchannel=str(nChan), iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time.pdf', expformat='pdf',highres=True,dpi=1200,exprange='all', showgui=False)
        os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time_* '+directory+'Pre-cal_phase_vs_time_'+vis+'_'+x[i]+'.pdf')
        os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_time_*')

        ## - Phase vs frequency
        plotms(vis=vis,xaxis='frequency',yaxis='phase',xdatacolumn='data',ydatacolumn='data',\
field=x[i], antenna='*&*', averagedata=True, avgtime=time, iteraxis='baseline', plotfile=directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency.pdf', expformat='pdf',highres=True,dpi=1200,exprange='all', showgui=False)
        os.system('convert '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency_* '+directory+'Pre-cal_phase_vs_frequency_'+vis+'_'+x[i]+'.pdf')
        os.system('rm '+directory+'pre-cal_'+vis+'_'+x[i]+'_phase_vs_frequency_*')
    #vishead(vis=vis,listfile=directory+vis+'.listobs')
    #plotants(vis=vis,figfile=directory+vis+'.plotants.png')
    ## Amplitude vs Time:
    os.system('mv ./'+plot_dir+'/*pdf ./'+plot_dir+'/pre-calibration')
    logger.info('End prediagnostics')


def flagdata1_apriori(msfile, msinfo, flags, do_quack=True):
    logger.info('Start flagdata1_apriori')
    # Find number of channels in MS:
    nchan = msinfo['nchan']
    # Flag Lo-Mk2
    if 'Lo' in msinfo['antennas'] and 'Mk2' in msinfo['antennas']:
        logger.info('Flagging Lo-Mk2 baseline')
        flagdata(vis=msfile, mode='manual', field=msinfo['sources']['allsources'], antenna='Lo*&Mk2*')
    # Subband edges
    channels_to_flag = '*:0~{0};{1}~{2}'.format(nchan/128-1, nchan-nchan/128, nchan-1)
    logger.info('MS has {} channels'.format(nchan))
    logger.info('Flagging edge channels {0}'.format(channels_to_flag))
    flagdata(vis=msfile, mode='manual', field=msinfo['sources']['allsources'], spw=channels_to_flag)
    # Slewing (typical):
    # Main calibrators, 5 min
    logger.info('Flagging 5 min from bright calibrators')
    flagdata(vis=msfile, field=msinfo['sources']['maincal'], mode='quack', quackinterval=300)
    # Target and phase reference, 20 sec
    logger.info('Flagging first 20 sec of target and phasecal')
    flagdata(vis=msfile, field=msinfo['sources']['targets_phscals'], mode='quack', quackinterval=20)
    # We can add more (Lo is slower, etc).
    flag_applied(flags, 'flagdata1_apriori')
    logger.info('End flagdata1_apriori')
    return flags

def flagdata2_manual(msfile, inpfile, flags):
    logger.info('Start flagdata_manual')
    logger.info('Applying manual flags from file: {0}'.format(inpfile))
    flagdata(vis=msfile, mode='list', inpfile=inpfile)
    flag_applied(flags, 'flagdata_manual')
    logger.info('End flagdata_manual')
    return flags

def flagdata3_tfcropBP(msfile, msinfo, flags):
    logger.info('Start flagdata3_tfcropBP')
    logger.info("Running flagdata, mode = 'tfcrop'")
    logger.info("correlation='ABS_ALL'")
    logger.info("ntime='90min', combinescans=True, datacolumn='corrected'")
    logger.info("winsize=3, timecutoff=4.0, freqcutoff=3.0, maxnpieces=1")
    logger.info("usewindowstats='sum', halfwin=3, extendflags=True")
    flagdata(vis=msfile, mode='tfcrop', field=msinfo['sources']['allsources'],
             antenna='', scan='',spw='', correlation='ABS_ALL',
             ntime='90min', combinescans=True, datacolumn='corrected',
             winsize=3, timecutoff=4.0, freqcutoff=3.0, maxnpieces=1,
             usewindowstats='sum', halfwin=3, extendflags=True,
             action='apply', display='', flagbackup=True)
    flag_applied(flags, 'flagdata3_tfcropBP')
    logger.info('End flagdata3_tfcropBP')
    return flags

def flagdata_tfcrop_bright(msfile, sources, flags, datacolumn='DATA'):
    logger.info('Start flagdata_tfcrop_bright')
    logger.info("Running flagdata, mode = 'tfcrop'")
    logger.info("correlation='ABS_ALL'")
    logger.info("ntime='90min', combinescans=True, datacolumn='{}'".format(datacolumn))
    logger.info("winsize=3, timecutoff=4.0, freqcutoff=3.0, maxnpieces=5")
    logger.info("usewindowstats='sum', halfwin=3, extendflags=True")
    flagdata(vis=msfile, mode='tfcrop', field=sources['allsources'],
             antenna='', scan='',spw='', correlation='ABS_ALL',
             ntime='90min', combinescans=True, datacolumn=datacolumn,
             winsize=3, timecutoff=3.6, freqcutoff=3.6, maxnpieces=2,
             usewindowstats='sum', halfwin=3, extendflags=True,
             action='apply', display='', flagbackup=True)
    flagdata(vis=msfile, mode='extend', ntime='90min', extendpols=False,
             growtime=50, growfreq=0)
    flag_applied(flags, 'flagdata_tfcrop_bright')
    logger.info('End flagdata_tfcrop_bright')
    return flags

def flagdata4_rflag(msfile, msinfo, flags):
    timedevscale = 5
    freqdevscale = 5
    logger.info('Start flagdata4_rflag')
    logger.info("Running flagdata, mode = 'rflag'")
    logger.info("ntime='90min', combinescans=True, datacolumn='corrected'")
    logger.info("timedevscale = {0}, freqdevscale = {1}".format(timedevscale,
                                                                freqdevscale))
    flagdata(vis=msfile, mode='rflag', field=msinfo['sources']['allsources'],
             antenna='', scan='',spw='', correlation='',
             ntime='90min', combinescans=True, datacolumn='corrected',
             timedevscale=timedevscale, freqdevscale=freqdevscale,
             action='apply', display='', flagbackup=True)
    flag_applied(flags, 'flagdata3_rflag')
    logger.info('End flagdata3_rflag')
    return flags


def flag_applied(flags, new_flag):
    if new_flag in flags:
        logger.warning('Flags from {0} were already applied by the pipeline! They are applied more than once.'.format(new_flag))
    flags.append(new_flag)
    logger.info('New flags applied: {0}'.format(new_flag))
    logger.info('Current flags applied: {0}'.format(flags))
    save_obj(flags, './flags')
    return flags

def define_refant(msfile, msinfo, inputs):
    refant0 = inputs['refant']
    refant_user = refant0.replace(' ', '').split(',')
    refant_in_ms = (np.array([ri in msinfo['antennas'] for ri in refant_user])).all()
    if not refant_in_ms:
        if refant0 != '':
            logger.warning('Selected reference antenna(s) {0} not in MS! User selection will be ignored'.format(refant0))
        # Finding best antennas for refant
        refant, refant_pref = find_refant(msfile, field=msinfo['sources']['bpcal'],
                                             antennas='Mk2,Pi,Da,Kn', spws='2,3', scan='')
    else:
        refant = ','.join(refant_user)
    logger.info('Refant: {}'.format(refant))
    return refant


def find_refant(msfile, field, antennas='', spws='', scan=''):
    logger.info('Searching refant automatically')
    ms.open(msfile)
    d = ms.getdata2(['axis_info'],ifraxis=True)
    ms.close()
    if len(antennas)==0:
        antennas = np.unique('-'.join(d['axis_info']['ifr_axis']['ifr_name']).split('-'))
    else:
        antennas = np.array(antennas.strip(' ').split(','))
    nchan = len(d['axis_info']['freq_axis']['chan_freq'][:,0])
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    channels = ':{0}~{1}'.format(nchan - 3*nchan/4, nchan - 1*nchan/4)
    if spws == '':
        spws = range(num_spw)
    else:
        spws = spws.strip(' ').split(',')

    logger.info('Searching best antenna among {0}'.format(','.join(antennas)))
    logger.info('Considering best S/N for spw {0}, channels {1}, field {2}'.format(','.join(spws), channels, field))

    snr = np.zeros((2, len(antennas), len(spws)))
    for i, corr in enumerate(['RR', 'LL']):
        for j, anten in enumerate(antennas):
            print('Checking antenna {0}, corr {1}'.format(anten, corr))
            for k, spw in enumerate(spws):
                snrj = []
                for basel in antennas:
                    if basel != anten:
                        try:
                            d2 = visstat(msfile, antenna='{}&{}'.format(anten, basel),
                                         field=field, axis='amplitude', scan=scan,
                                         spw=str(spw)+channels, correlation=corr)
                            snrj.append(d2['DATA']['median']/d2['DATA']['stddev'])
                        except:
                            snrj.append(0.0)
                snr[i, j, k] = np.mean(snrj)

    # Sorting antennas by average S/N, averaged in polarization and spw.
    snr_avg = np.mean(snr, axis=(0,2))
    pref_ant = antennas[np.argsort(snr_avg)][::-1]
    snr_avg_sorted = snr_avg[np.argsort(snr_avg)][::-1]

    logger.info('Average S/N on baselines to each antenna:')
    for p, s in zip(pref_ant, snr_avg_sorted):
        logger.info('{0:3s}  {1:4.1f}'.format(p, s))

    # This will check if Mk2, Pi or Da are part of the best three antennas and give them
    # priority. Then, it will add any other except Lo.
    u, ind = np.unique([a for a  in pref_ant[:3] if a in ['Mk2','Pi','Da']] +\
    [a for a  in pref_ant if a in ['Mk2','Pi','Da','Kn', 'Cm', 'De']], return_index=True)

    refant = ','.join(u[np.argsort(ind)])
    refant_pref = {}
    refant_pref['snr'] = snr
    refant_pref['pref_ant'] = pref_ant
    refant_pref['snr_avg_sorted'] = snr_avg_sorted
    logger.info('Preference antennas refant = {0}'.format(refant))
    return refant, refant_pref


### Run CASA calibration functions

def run_split(msfile, msinfo, width, timebin, datacolumn='data'):
    logger.info('Start split')
    # Check that all sources are there:
    sources_not_in_msfile = [s for s in
                             msinfo['sources']['allsources'].split(',')
                             if s not in msinfo['mssources'].split(',')]
    if len(sources_not_in_msfile) > 0:
        fields = msinfo['mssources']
        logger.warning('Fields {} not present in MS but listed in inputs file.'.format(','.join(sources_not_in_msfile)))
        logger.warning('All fields will be included in the averaged MS.')
    else:
        fields = msinfo['sources']['allsources']
    name = '.'.join(msfile.split('.')[:-1])
    exte = ''.join(msfile.split('.')[-1])
    outputmsfile = name+'_avg.'+exte
    rmdir(outputmsfile)
    rmdir(outputmsfile+'.flagversions')
    rmdir(outputmsfile+'.listobs')
    logger.info('Input MS: {0}'.format(msfile))
    logger.info('Output MS: {0}'.format(outputmsfile))
    logger.info('width={0}, timebin={1}'.format(width, timebin,))
    logger.info('Fields: {0}'.format(fields))
    logger.info('Data column: {0}'.format(datacolumn))
    split(vis=msfile, outputvis=outputmsfile, field=fields, width=width,
          timebin=timebin, datacolumn=datacolumn, keepflags=True)
    listobs(vis=outputmsfile, listfile=outputmsfile+'.listobs',overwrite=True)
    logger.info('Listobs file in: {0}'.format(outputmsfile+'.listobs'))
    logger.info('End split')


def run_initialize_models(msfile, fluxcal, models_path, delmod_sources):
    logger.info('Start init_models')
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
    logger.info('Initializing 3C286 model using: {0}'.format(model_3C286))
    if fluxcal != '1331+305':
        logger.warning('Using a model for 3C286 (1331+305) but your flux calibrator source is: {0}. Model may be wrong for that source'.format(fluxcal))
    setjy(vis=msfile, field=fluxcal, standard='Perley-Butler 2013',
          model=model_3C286, scalebychan=True, usescratch=True)
    # Is usescratch needed? Probably not, but I see amp=1 in the model column
    # otherwise when running on an MMS file.
    logger.info('Deleting model for all other sources: '+delmod_sources)
    delmod(vis=msfile, field=delmod_sources)
    logger.info('End init_models')

def run_fringefit(msfile, caltables, caltable_name, previous_cal, minblperant=3, minsnr=2, smodel=[]):
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
    gainfield = [caltables[p]['field'] if len(np.atleast_1d(caltables[p]['field'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating calibration table: {0}'.format(caltables[caltable_name]['table']))
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


def run_gaincal(msfile, caltables, caltable_name, previous_cal, minblperant=3, minsnr=2):
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
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['field'] if len(np.atleast_1d(caltables[p]['field'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating calibration table: {0}'.format(caltables[caltable_name]['table']))
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
            gaintable = gaintable,
            gainfield = gainfield,
            interp    = interp,
            spwmap    = spwmap,
            minblperant=minblperant,
            minsnr=minsnr)
    logger.info('caltable {0} in {1}'.format(caltables[caltable_name]['name'],
                                              caltables[caltable_name]['table']))

def run_bandpass(msfile, caltables, caltable_name, previous_cal, minblperant=3, minsnr=2):
    rmdir(caltables[caltable_name]['table'])
    logger.info('Running bandpass to generate: {0}'.format(caltables[caltable_name]['name']))
    logger.info('Field(s) {0}, solint = {1}, spw = {2}, combine = {3}, solnorm = {4}'.format(
                 caltables[caltable_name]['field'],
                 caltables[caltable_name]['solint'],
                 caltables[caltable_name]['spw'],
                 caltables[caltable_name]['combine'],
                 caltables[caltable_name]['solnorm']))
    logger.info('uvrange = {0}'.format(caltables[caltable_name]['uvrange']))
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['field'] if len(np.atleast_1d(caltables[p]['field'].split(',')))<2 else '' for p in previous_cal]
    logger.info('Previous calibration applied: {0}'.format(str(previous_cal)))
    logger.info('Previous calibration gainfield: {0}'.format(str(gainfield)))
    logger.info('Previous calibration spwmap: {0}'.format(str(spwmap)))
    logger.info('Previous calibration interp: {0}'.format(str(interp)))
    logger.info('Generating calibration table: {0}'.format(caltables[caltable_name]['table']))
    # Run CASA task bandpass
    bandpass(vis=msfile,
             caltable  = caltables[caltable_name]['table'],
             field     = caltables[caltable_name]['field'],
             solint    = caltables[caltable_name]['solint'],
             combine   = caltables[caltable_name]['combine'],
             spw       = caltables[caltable_name]['spw'],
             solnorm   = caltables[caltable_name]['solnorm'],
             uvrange   = caltables[caltable_name]['uvrange'],
             fillgaps  = 16,
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
    plotcal(caltable=tablein, xaxis='time', yaxis='phase', subplot=321, iteration='antenna', showgui=False, figfile=caltableplot_phs, fontsize = 8, plotrange=[-1,-1,-180,180])
    plotcal(caltable=tablein, xaxis='time', yaxis='amp', subplot=321, iteration='antenna', showgui=False, figfile=caltableplot_amp, fontsize = 8, plotrange=[-1,-1,-1,-1])
    return


def run_applycal(msfile, caltables, sources, previous_cal, previous_cal_targets=''):
    logger.info('Start applycal')
    # 1 correct non-target sources:
    logger.info('Applying calibration to calibrator sources')
    logger.info('Fields: {0}'.format(sources['calsources']))
    # Previous calibration
    gaintable = [caltables[p]['table'] for p in previous_cal]
    interp    = [caltables[p]['interp'] for p in previous_cal]
    spwmap    = [caltables[p]['spwmap'] for p in previous_cal]
    gainfield = [caltables[p]['field'] if len(np.atleast_1d(caltables[p]['field'].split(',')))<2 else '' for p in previous_cal]
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
    logger.info('Applying calibration to target sources')
    logger.info('Target fields: {0}'.format(sources['targets']))
    if previous_cal_targets == '':
        previous_cal_targets = previous_cal
    for i, s in enumerate(sources['targets'].split(',')):
        if s != '':
            phscal = sources['phscals'].split(',')[i]
            # Previous calibration
            gaintable = [caltables[p]['table'] for p in previous_cal_targets]
            interp    = [caltables[p]['interp'] for p in previous_cal_targets]
            spwmap    = [caltables[p]['spwmap'] for p in previous_cal_targets]
            gainfield = [caltables[p]['field'] if len(np.atleast_1d(caltables[p]['field'].split(',')))<2 else phscal for p in previous_cal_targets]
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
    logger.info('End applycal')


### Calibration steps

def solve_delays(msfile, caltables, previous_cal, calsources, solint='300s'):
    logger.info('Start solve_delays')
    caltable_name = 'delay.K1'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = calsources
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = solint
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = [0]*caltables['num_spw']
    caltables[caltable_name]['combine'] = 'spw'
    caltables[caltable_name]['spw'] = ''

    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal)
    logger.info('Delay calibration {0}: {1}'.format(caltable_name, caltable))
    # Plots
    # 1 No range
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_1.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange = [-1,-1,-1,-1])
    logger.info('Delay calibration plot in: {0}'.format(caltableplot))
    # 2 Only show typical delay values: from -20 to 20 nanosec
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_2.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange = [-1,-1,-20,20])
    logger.info('Delay calibration plot: {0}'.format(caltableplot))
    logger.info('End solve_delays')
    return caltables

def delay_fringefit(msfile, caltables, previous_cal, calsources):
    # Currently does not combine spws because that is not correctly implemented
    # in the pre-release version of task fringefit
    logger.info('Start delay_fringefit')
    caltable_name = 'delay.K1'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = calsources
    caltables[caltable_name]['solint'] = '300s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltables[caltable_name]['zerorates'] = True

    caltable = caltables[caltable_name]['table']
    # Calibration
    run_fringefit(msfile, caltables, caltable_name, previous_cal)
    logger.info('Fringe calibration {0}: {1}'.format(caltable_name, caltable))
    # Plots
    # 1 Phases
    caltableplot =  caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange =
            [-1,-1,-180,-180])
    logger.info('Fringe calibration (phases) plot in: {0}'.format(caltableplot))
    # 2 Delays
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_dela.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange =
            [-1,-1,-1,-1])
    logger.info('Fringe calibration (delays) plot in: {0}'.format(caltableplot))
    caltableplot =  caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_dela2.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange =
            [-1,-1,-20,-20])
    logger.info('Fringe calibration (delays range) plot in: {0}'.format(caltableplot))
    ## 3 Rates (they are zeroed)
    #caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_rate.png'
    #plotcal(caltable=caltable,xaxis='time',yaxis='rate',subplot=321,iteration='antenna',
    #        showgui=False,figfile=caltableplot, fontsize = 8, plotrange =
    #        [-1,-1,-1,-1])
    #logger.info('Fringe calibration (rates) plot in: {0}'.format(caltableplot))

    # 2 Only show typical delay values: from -20 to 20 nanosec
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_2.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange = [-1,-1,-20,20])
    logger.info('Delay calibration plot: {0}'.format(caltableplot))
    logger.info('End delay_fringefit')
    return caltables


def initial_bp_cal(msfile, caltables, previous_cal, bpcal):
    logger.info('Start initial_bp_cal')

    # 0 Delay calibration of bpcal
    caltable_name = 'bpcal_d.K0'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = bpcal
    caltables[caltable_name]['gaintype'] = 'K'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = '180s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = [0]*caltables['num_spw']
    caltables[caltable_name]['combine'] = 'spw'
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal)
    logger.info('Delay calibration of bpcal {0}: {1}'.format(caltable_name, caltable))
    # Plots
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_1.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='delay',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot, fontsize = 8, plotrange = [-1,-1,-1,-1])
    logger.info('Delay calibration plot in: {0}'.format(caltableplot))

    # 1 Phase calibration
    caltable_name = 'bpcal_p.G0'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = bpcal
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = '8s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    previous_cal_p = previous_cal + ['bpcal_d.K0']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal_p)
    logger.info('Bandpass0 phase calibration {0}: {1}'.format(caltable_name,caltable))
    # Plots
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot,fontsize=8, plotrange=[-1,-1,-180,180])
    logger.info('Bandpass0 phase calibration plot: {0}'.format(caltableplot))

    # 2 Amplitude calibration
    caltable_name = 'bpcal_ap.G1'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = bpcal
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['solint'] = '32s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    previous_cal_ap = previous_cal_p + ['bpcal_p.G0']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal_ap)
    logger.info('Bandpass0 amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_phs,fontsize=8, plotrange=[-1,-1,-180,180])
    logger.info('Bandpass0 amplitude calibration plot: {0}'.format(caltableplot_phs))
    caltableplot_amp = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='amp',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_amp,fontsize=8, plotrange=[-1,-1,-1,-1])
    logger.info('Bandpass0 amplitude calibration plot: {0}'.format(caltableplot_amp))

    # 3 Bandpass calibration
    caltable_name = 'bpcal.B0'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = bpcal
    caltables[caltable_name]['solint'] = 'inf'
    caltables[caltable_name]['interp'] = 'nearest,linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltables[caltable_name]['uvrange'] = ''
    caltables[caltable_name]['solnorm'] = True
    bptable = caltables[caltable_name]['table']
    previous_cal_ap_bp = previous_cal_ap + ['bpcal_ap.G1']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name, previous_cal_ap_bp)
    logger.info('Bandpass0 BP {0}: {1}'.format(caltable_name,bptable))
    # Plots
    bptableplot_phs = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    bptableplot_amp = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plotcal(caltable=bptable, xaxis='freq', yaxis='phase',
            subplot=321,iteration='antenna', showgui=False,
            figfile=bptableplot_phs, fontsize = 8, plotrange = [-1,-1,-180,180])
    logger.info('Bandpass0 BP phase plot: {0}'.format(caltableplot_phs))
    plotcal(caltable=bptable, xaxis='freq', yaxis='amp',  subplot=321,
            iteration='antenna', showgui=False, figfile=bptableplot_amp,
            fontsize = 8, plotrange = [-1,-1,-1,-1])
    logger.info('Bandpass0 BP phase plot: {0}'.format(caltableplot_amp))
    logger.info('End initial_bp_cal')
    return caltables


def initial_gaincal(msfile, caltables, previous_cal, calsources, phscals):
    logger.info('Start initial_gaincal')

    # 1 Phase calibration
    caltable_name = 'allcal_p.G0'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = calsources
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = '8s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal)
    logger.info('Gain phase calibration {0}: {1}'.format(caltable_name,caltable))
    # Plots
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',showgui=False,figfile=caltableplot, fontsize = 8, plotrange = [-1,-1,-180, 180])
    logger.info('Gain phase calibration plot: {0}'.format(caltableplot))

    # 2 Amplitude calibration
    caltable_name = 'allcal_ap.G1'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = calsources
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['solint'] = '32s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    previous_cal_ap = previous_cal + ['allcal_p.G0']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal_ap)
    logger.info('Gain amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_phs,fontsize=8,plotrange=[-1,-1,-180,180])
    logger.info('Bandpass0 phase calibration plot: {0}'.format(caltableplot))
    caltableplot_amp = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='amp',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_amp,fontsize=8,plotrange=[-1,-1,-1,-1])
    logger.info('Bandpass0 phase calibration plot: {0}'.format(caltableplot))

    # 3 Phase calibration on phasecal: scan-averaged phase solutions
    caltable_name = 'phscal_p_scan.G2'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = phscals
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'p'
    caltables[caltable_name]['solint'] = 'inf'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal)
    logger.info('Gain phase calibration {0}: {1}'.format(caltable_name,caltable))
    # Plots
    caltableplot = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot,fontsize=8,plotrange=[-1,-1,-180,180])
    logger.info('Gain phase calibration plot: {0}'.format(caltableplot))
    logger.info('End initial_gaincal')
    return caltables


def find_anten_fluxscale(antennas):
    # This function tries to remove Lo and De from the fluxscale determination.
    # But only if there are enough antennas to have at least 4 of them.
    # I found unstable solutions if too few antennas are used
    if len(antennas) > 4:
         anten_for_flux = [x for x in antennas if x != "Lo"]
         if len(anten_for_flux) > 4:
             anten_for_flux = [x for x in anten_for_flux if x != "De"]
    else:
        anten_for_flux = antennas
    return anten_for_flux


def eM_fluxscale(msfile, caltables, ampcal_table, sources, antennas):
    logger.info('Start eM_fluxscale')
    anten_for_flux = find_anten_fluxscale(antennas)
    cals_to_scale = sources['cals_no_fluxcal']
    fluxcal = sources['fluxcal']
    caltable_name = 'allcal_ap.G1_fluxscaled'
    caltables[caltable_name] = copy.copy(caltables[ampcal_table])
    caltables[caltable_name]['table']=caltables[ampcal_table]['table']+'_fluxscaled'
    logger.info('Flux density scale from: {0}'.format(fluxcal))
    logger.info('Transfered to: {0}'.format(cals_to_scale))
    logger.info('Input caltable: {0}'.format(caltables[ampcal_table]['table']))
    logger.info('Antennas used to scale: {0}'.format(anten_for_flux))
    calfluxes = fluxscale(vis=msfile, reference=fluxcal,
                          transfer=cals_to_scale,
                          antenna = ','.join(anten_for_flux),
                          caltable  = caltables[ampcal_table]['table'],
                          fluxtable = caltables[caltable_name]['table'],
                          listfile  = caltables[caltable_name]['table']+'_fluxes.txt')
    logger.info('Modified caltable: {0}'.format(caltables[caltable_name]['table']))
    logger.info('Spectrum information: {0}'.format(caltables[caltable_name]['table']+'_fluxes.txt'))
    # Compute correction to scale the flux density of 3C286 according to
    # resolution provided by the shortest available baseline of e-MERLIN
    eMfactor = calc_eMfactor(msfile, field=fluxcal)
    # Include a note in the fluxes.txt file warning that the values in that
    # file should be corrected by eMfactor
    with open(caltables[caltable_name]['table']+'_fluxes.txt', 'a') as file:
        file.write('# WARNING: All flux densities in this file need to be multiplied by eMfactor={0:6.4f} to match the corrections that have been applied to the data.'.format(eMfactor))
    # Get fitted flux density and spectral index, correctly scaled for e-MERLIN
    eMcalfluxes = {}
    for k in calfluxes.keys():
        if len(calfluxes[k]) > 4:
            try:
                a=[]
                a.append(calfluxes[k]['fitFluxd']*eMfactor)
                a.append(calfluxes[k]['spidx'][1])
                a.append(calfluxes[k]['fitRefFreq'])
                eMcalfluxes[calfluxes[k]['fieldName']]=a
                logger.info('Spectrum for {0:>9s}: Flux density ={1:6.3f} +/-{2:6.3f}, spidx ={3:5.2f}+/-{4:5.2f}'.format(calfluxes[k]['fieldName'],
                    calfluxes[k]['fitFluxd']*eMfactor, calfluxes[k]['fitFluxdErr']*eMfactor,
                    calfluxes[k]['spidx'][1], calfluxes[k]['spidxerr'][1]))
            except:
                pass

    for f in eMcalfluxes.keys():    # Phase calibrator and bandpass calibrator
        setjy(vis = msfile,
              field = f,
              standard = 'manual',
              fluxdensity = eMcalfluxes[f][0],
              spix = eMcalfluxes[f][1],
              reffreq = str(eMcalfluxes[f][2])+'Hz')
    logger.info('End eM_fluxscale')
    return caltables

def bandpass_sp(msfile, caltables, previous_cal, bpcal):
    logger.info('Start bandpass_sp')
    # Bandpass calibration
    caltable_name = 'bpcal_sp.B1'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = bpcal
    caltables[caltable_name]['solint'] = 'inf'
    caltables[caltable_name]['interp'] = 'nearest,linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    #caltables[caltable_name]['uvrange'] = '>15km'
    caltables[caltable_name]['uvrange'] = ''
    caltables[caltable_name]['solnorm'] = False
    bptable = caltables[caltable_name]['table']
    # Calibration
    run_bandpass(msfile, caltables, caltable_name, previous_cal)
    logger.info('Bandpass1 BP {0}: {1}'.format(caltable_name,bptable))
    # Plots
    bptableplot_phs = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    bptableplot_amp = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plotcal(caltable=bptable, xaxis='freq', yaxis='phase',
            subplot=321,iteration='antenna', showgui=False,
            figfile=bptableplot_phs, fontsize = 8, plotrange = [-1,-1,-180,180])
    logger.info('Bandpass1 BP phase plot: {0}'.format(bptableplot_phs))
    plotcal(caltable=bptable, xaxis='freq', yaxis='amp',  subplot=321,
            iteration='antenna', showgui=False, figfile=bptableplot_amp,
            fontsize = 8, plotrange = [-1,-1,-1,-1])
    logger.info('Bandpass1 BP phase plot: {0}'.format(bptableplot_amp))
    logger.info('End bandpass_sp')
    return caltables


def sp_amp_gaincal(msfile, caltables, previous_cal, calsources):
    logger.info('Start gaincal_amp_sp')

    # 1 Amplitude calibration
    caltable_name = 'allcal_ap.G3'
    caltables[caltable_name] = {}
    caltables[caltable_name]['name'] = caltable_name
    caltables[caltable_name]['table'] = caltables['calib_dir']+caltables['inbase']+'_'+caltable_name
    caltables[caltable_name]['field'] = calsources
    caltables[caltable_name]['gaintype'] = 'G'
    caltables[caltable_name]['calmode'] = 'ap'
    caltables[caltable_name]['solint'] = '32s'
    caltables[caltable_name]['interp'] = 'linear'
    caltables[caltable_name]['spwmap'] = []
    caltables[caltable_name]['combine'] = ''
    caltables[caltable_name]['spw'] = ''
    caltable = caltables[caltable_name]['table']
    # Calibration
    run_gaincal(msfile, caltables, caltable_name, previous_cal)
    logger.info('Gain amplitude calibration {0}: {1}'.format(caltable_name,caltable))
#    smooth_caltable(msfile=msfile, plotdir=plotdir, tablein=caltable2, caltable='', field='', smoothtype='median', smoothtime=60*20.)
    # Plots
    caltableplot_phs = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_phs.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='phase',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_phs,fontsize=8,plotrange=[-1,-1,-180,180])
    logger.info('Amplitude gain calibration plot: {0}'.format(caltableplot_phs))
    caltableplot_amp = caltables['plots_dir']+caltables['inbase']+'_'+caltable_name+'_amp.png'
    plotcal(caltable=caltable,xaxis='time',yaxis='amp',subplot=321,iteration='antenna',
            showgui=False,figfile=caltableplot_amp,fontsize=8,plotrange=[-1,-1,-1,-1])
    logger.info('Amplitude gain calibration plot: {0}'.format(caltableplot_amp))

    logger.info('End gaincal_amp_sp')
    return caltables

def compile_statistics(msfile, tablename=''):
    logger.info('Start compile_stats')
    # Num of spw and baselines
    num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
    baselines = get_baselines(msfile)
    # Date and time of observation
    ms.open(msfile)
    axis_info = ms.getdata2(['axis_info'],ifraxis=True)
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


def monitoring(msfile, msinfo, flags, caltables, previous_cal):
    # This is intented to run on a single-source file for daily monitoring on
    # unaveraged data
    logger.info('Starting monitoring')
    band = check_band(msfile)
    if band == 'L':
        hanning(inputvis=msfile,deloriginal=True)

    flags = flagdata1_apriori(msfile=msfile, msinfo=msinfo, flags=flags, do_quack=True)

    flags = flagdata_tfcrop_bright(msfile=msfile, sources=msinfo['sources'], flags=flags)
    caltables = solve_delays(msfile=msfile, caltables=caltables,
                   previous_cal=[], calsources=msinfo['sources']['calsources'],
                             solint='60s')
    run_applycal(msfile=msfile, caltables=caltables, sources=msinfo['sources'],
                    previous_cal=['delay.K1'],
                    previous_cal_targets=['delay.K1'])
    compile_statistics(msfile, tablename=caltables['delay.K1']['table'])
    return flags, caltables


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

    tb.open(msfile)
    uvw = tb.getcol('UVW')
    a1 = tb.getcol('ANTENNA1')
    a2 = tb.getcol('ANTENNA2')
    field = tb.getcol('FIELD_ID')
    flags = tb.getcol('FLAG')
    tb.close()

    uvdist = np.sqrt(uvw[0]**2+uvw[1]**2)
    allflag = np.sum(flags, axis=(0,1))
    flags_entries = flags.shape[0]*flags.shape[1]

    mask = (uvdist==0) | (allflag==flags_entries) | field!=field_id
    uvdist_nonzero = np.ma.array(uvdist, mask=mask)

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

