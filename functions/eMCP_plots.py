#!/usr/local/python
import os
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
plt.ioff()

from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.ticker import ScalarFormatter
import datetime
import shutil
import glob

from functions import eMCP_weblog as emwlog
from functions import eMCP_utils as emutils
from functions import eMCP_functions as em


#import casatasks
#from casaplotms import plotms
#from casatools import table
#from casatools import msmetadata
#from casatools import ms as myms
#from casatools import measures
#tb = table()
#msmd = msmetadata()
#ms = myms()
#me = measures()

## CASA imports
#from taskinit import *
#from tasks import *
#from casac import casac
#msmd = casac.msmetadata()

#plt.ioff()

import logging
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

line0 = '-'*15


def add_step_time(step, eMCP, msg, t0, doweblog=True):
    t1 = datetime.datetime.utcnow()
    timestamp = t1.strftime('%Y-%m-%d %H:%M:%S')
    delta_t_min = (t1-t0).total_seconds()/60.
    eMCP['steps'][step] = [timestamp, delta_t_min, msg]
    emutils.save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    os.system('cp casa_eMCP.log {}casa_eMCP.log.txt'.format(info_dir))
    if doweblog:
        emwlog.start_weblog(eMCP)
    return eMCP

#'#def get_scans(msfile, field):
#'#    ms.open(msfile)
#'#    ms.msselect({'field':field})
#'#    scans = ms.getdata('scan_number')['scan_number']
#'#    num_scans = len(np.unique(scans))
#'#    ms.close()
#'#    return num_scans, scans

#'#def get_freqs(msfile, allfreqs=False):
#'#    if allfreqs:
#'#        channels = emutils.read_keyword(msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW')
#'#        freqs=np.sort(channels.flatten())[::int(len(channels.T)/4)]
#'#    else:
#'#        freqs = channels.mean(axis=0)
#'#    return freqs

#'#def get_freqs(msfile, allfreqs=False):
#'#    ms.open(msfile)
#'#    axis_info = ms.getdata(['axis_info'],ifraxis=True)
#'#    ms.close()
#'#    if allfreqs:
#'#        channels = axis_info['axis_info']['freq_axis']['chan_freq']
#'#        freqs=np.sort(axis_info['axis_info']['freq_axis']['chan_freq'].flatten())[::int(len(channels)/4)]
#'#    else:
#'#        freqs = axis_info['axis_info']['freq_axis']['chan_freq'].mean(axis=0)
#'#    return freqs

#'#def count_active_baselines(msfile):
#'#    msmd.open(msfile)
#'#    a = np.array(msmd.antennanames())
#'#    b = np.array(msmd.baselines())
#'#    msmd.done()
#'#    active_baselines = 0
#'#    for i, ant1 in enumerate(a):
#'#        for ant2 in a[i+1:]:
#'#            j = np.argwhere(ant2==a)[0][0]
#'#            if b[i,j]:
#'#                active_baselines += 1
#'#    return active_baselines

def simple_plot_name(plot_file, i):
    try:
        actual_name = glob.glob('{0}{1}_*.png'.format(plot_file, i))[0]
        shutil.move(actual_name, '{0}{1}.png'.format(plot_file, i))
    except:
        pass

def count_active_antennas(msfile):
    msmd.open(msfile)
    a = np.array(msmd.antennanames())
    b = np.array(msmd.baselines())
    msmd.done()
    active_antennas = []
    for i, ant1 in enumerate(a):
        for ant2 in a[i+1:]:
            j = np.argwhere(ant2==a)[0][0]
            if b[i,j]:
                active_antennas.append(ant1)
                active_antennas.append(ant2)
    nice_order = ['Lo', 'Mk2', 'Pi', 'Da', 'Kn', 'De', 'Cm']
    active_antennas = [a for a in nice_order if a in active_antennas]
    return active_antennas

def plot_caltable(msinfo, caltable, plot_file, xaxis='', yaxis='', title='',
                  ymin=-1, ymax=-1, coloraxis='spw', symbolsize=8):
    gridcols = 1
    showgui=False
    tab = caltable['table']
    if xaxis == 'time':
        tb.open(tab)
        time_mjd = tb.getcol('TIME')
        tb.close()
        x_min, x_max = np.min(time_mjd), np.max(time_mjd)
    elif xaxis == 'freq':
        tb.open(tab+'/SPECTRAL_WINDOW')
        f = tb.getcol('CHAN_FREQ').flatten()/1e9
        tb.close()
        x_min0, x_max0 = np.min(f), np.max(f)
        x_span = x_max0 - x_min0
        x_min = x_min0 - x_span*0.1
        x_max = x_max0 + x_span*0.1
    else:
        x_min, x_max = -1, -1

    active_antennas = count_active_antennas(msinfo['msfile'])
    num_anten = len(active_antennas)
    for i, anten in enumerate(active_antennas):
        if i == len(active_antennas)-1:
            plotfile = plot_file
        else:
            plotfile = ''
        plotms(vis=caltable['table'], xaxis=xaxis, yaxis=yaxis, title='{0} {1}'.format(title, anten),
               gridrows=num_anten, gridcols=gridcols, rowindex=i, colindex=0, plotindex=i,
               #timerange='{}~{}'.format(msinfo['t_ini'].time(), msinfo['t_end'].time()),
               antenna = str(anten),
               xselfscale = True, xsharedaxis = True, coloraxis = coloraxis,
               plotrange=[x_min, x_max, ymin, ymax],
               plotfile = plotfile, expformat = 'png', customsymbol = True,
               symbolshape = 'circle', symbolsize=symbolsize,
               width=1000, height=240*num_anten, clearplots=False, overwrite=True,
               showgui=showgui)


def single_4plot(msinfo, field, datacolumn, plots_data_dir):
    logger.info('Visibility plots for field: {0}, datacolumn: {1}'.format(field, datacolumn))
    msfile = msinfo['msfile']
    nchan = str(msinfo['nchan'])
#    tb.open(msfile+'/SPECTRAL_WINDOW')
#    nchan = str(tb.getcol('NUM_CHAN')[0])
#    tb.close()
    plot_file = os.path.join(plots_data_dir, f"{msinfo['msfilename']}_4plot_{field}_{datacolumn}")
    num_baselines = len(msinfo['baselines'])
#'#    num_baselines = count_active_baselines(msfile)
    gridrows = num_baselines
    gridcols = 1
    showgui=False
    avgtime = '300'
    baseline = '*&*'
    # Find min and max times per field
    x_min_time, x_max_time = emutils.find_source_timerange(msfile, field)
#'#    x_min_time = times.min()
#'#    x_max_time = times.max()
#'#    msmd.open(msfile)
#'#    all_times = msmd.timesforscans(msmd.scansforfield(field))
#'#    x_min_time = np.min(all_times)
#'#    x_max_time = np.max(all_times)
#'#    msmd.close()
    w, h = 1200, num_baselines*200
    # 0
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'time', 'yaxis': 'amp',
            'title':  f"Amp vs Time {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':baseline, 'field':field, 'iteraxis':'baseline',
            'averagedata':True, 'avgchannel':nchan, 'avgtime':'4',
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotrange':[x_min_time,x_max_time,0,-1],
            'plotfile': plot_file+'0.png', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()

    simple_plot_name(plot_file, 0)

    # 1
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'time', 'yaxis': 'phase',
            'title':  f"Phase vs Time {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':baseline, 'field':field, 'iteraxis':'baseline',
            'averagedata':True, 'avgchannel':nchan, 'avgtime':'4',
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotrange':[x_min_time,x_max_time,-180,180],
            'plotfile': plot_file+'1.png', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()
    simple_plot_name(plot_file, 1)

    # 2
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'freq', 'yaxis': 'amp',
            'title':  f"Amp vs Freq {field} (color=corr)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':baseline, 'field':field, 'iteraxis':'baseline',
            'averagedata':True, 'avgtime':avgtime, 'avgchannel':'4',
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis':'corr',
            'plotrange':[-1,-1,0,-1],
            'plotfile': plot_file+'2.png', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()
    simple_plot_name(plot_file, 2)

#'#    plotms(vis=msfile, xaxis='freq', yaxis='amp', title='Amp vs Frequency {0} (color=corr)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
#'#    antenna=baseline, field=field, iteraxis='baseline',
#'#    averagedata = True, avgtime=avgtime, avgchannel='4',
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,0,-1],
#'#    plotfile = plot_file+'2.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)
#'#    simple_plot_name(plot_file, 2)

    # 3
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'freq', 'yaxis': 'phase',
            'title':  f"Phase vs Freq {field} (color=corr)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':baseline, 'field':field, 'iteraxis':'baseline',
            'averagedata':True, 'avgtime':avgtime, 'avgchannel':'4',
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'corr',
            'plotrange':[-1,-1,-180,180],
            'plotfile': plot_file+'3.png', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()
#'#    plotms(vis=msfile, xaxis='freq', yaxis='phase', title='Phase vs Frequency {0} (color=corr)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
#'#    antenna=baseline, field=field, iteraxis='baseline',
#'#    averagedata = True, avgtime=avgtime, avgchannel='4',
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,-180,180],
#'#    plotfile = plot_file+'3.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)
#    logger.info('Finished {0}:'.format(field))
    simple_plot_name(plot_file, 3)
    logger.info('{0}{{0-4}}.png'.format(plot_file))


def make_4plots(eMCP, datacolumn='data'):
    logger.info(line0)
    msinfo = eMCP['msinfo']
    msfile = eMCP['msinfo']['msfile']
    logger.info('Start plot_{}'.format(datacolumn))
    t0 = datetime.datetime.utcnow()
    if datacolumn == 'data':
        plots_data_dir = './weblog/plots/plots_data/'
    elif datacolumn == 'corrected':
        plots_data_dir = './weblog/plots/plots_corrected/'
    else:
        plots_data_dir = './weblog/plots/'
    emutils.makedir(plots_data_dir)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    logger.info('Producing plots for: {}'.format(','.join(allsources)))
    for field in allsources:
        if field in mssources:
            single_4plot(msinfo, field, datacolumn, plots_data_dir)
        else:
            logger.warning('Cannot plot {0}. Source not in ms.'.format(field))
##    num_proc = eMCP['defaults']['plot_data']['num_proc']
##    pool = multiprocessing.Pool(num_proc)
##    input_args = [(msinfo, field, datacolumn, plots_data_dir) for field in fields]
##    p = pool.map(single_4plot, input_args)
##    pool.close()
##    pool.join()
    logger.info('Visibility plots finished')
    if datacolumn == 'corrected':
        make_uvplt(eMCP)
    logger.info('End plot_{}'.format(datacolumn))
    msg = ''
    eMCP = add_step_time('plot_'+datacolumn, eMCP, msg, t0)
    return eMCP

def single_uvplt(msinfo, field, plots_data_dir):
    logger.info('uvplt for field: {}'.format(field))
    plot_file =  plots_data_dir+'{0}_uvplt_{1}.png'.format(msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn='corrected'
    avgtime = '16'
    showgui = False
    gridrows = 1
    gridcols = 2
    # Amp
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'UVwave', 'yaxis': 'amp',
            'title':  f"Amplitude vs UVWave {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':'*&*', 'field':field,
            'averagedata':True, 'avgchannel':str(nchan), 'avgtime':avgtime,
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotfile': '', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()
#'#    plotms(vis=msfile, xaxis='UVwave', yaxis='amp', title='Amplitude vs UVWave {0} (color=spw)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
#'#    antenna='*&*', field=field,
#'#    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw',
#'#    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    symbolsize=4, clearplots=True, overwrite=True, showgui=showgui)

    # Phase
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'UVwave', 'yaxis': 'phase',
            'title':  f"Phase vs UVWave {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':1, 'plotindex':1,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':'*&*', 'field':field,
            'averagedata':True, 'avgchannel':str(nchan), 'avgtime':avgtime,
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotrange':[-1,-1,-180,180],
            'plotfile': plot_file, 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':1200, 'height':573, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()

#'#    plotms(vis=msfile, xaxis='UVwave', yaxis='phase', title='Phase vs UVWave {0} (color=spw)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=1,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
#'#    antenna='*&*', field=field,
#'#    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
#'#    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    width=1200, height=573, symbolsize=4, clearplots=False, overwrite=True, showgui=showgui)

def single_uvplt_model(msinfo, field, plots_data_dir):
    logger.info('uvplt (model) for field: {}'.format(field))
    plot_file =  plots_data_dir+'{0}_uvpltmodel_{1}.png'.format(msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn='model'
    avgtime = '600'
    showgui = False
    gridrows = 1
    gridcols = 2
    # Amp
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'UVwave', 'yaxis': 'amp',
            'title':  f"Amplitude vs UVWave {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':0, 'plotindex':0,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':'*&*', 'field':field,
            'averagedata':True, 'avgchannel':str(int(nchan/16)), 'avgtime':avgtime,
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotfile': '', 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':w, 'height':h, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()
#'#
#'#    plotms(vis=msfile, xaxis='UVwave', yaxis='amp', title='Model Amplitude vs UVWave {0} (color=spw)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
#'#    antenna='*&*', field=field,
#'#    averagedata = True, avgtime=avgtime, avgchannel = str(int(nchan/16)),
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw',
#'#    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    symbolsize=4, clearplots=True, overwrite=True, showgui=showgui)
#'#
    # Phase
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'UVwave', 'yaxis': 'phase',
            'title':  f"Phase vs UVWave {field} (color=spw)",
            'gridrows':gridrows, 'gridcols':gridcols, 
            'rowindex':0, 'colindex':1, 'plotindex':1,
            'xdatacolumn':datacolumn, 'ydatacolumn':datacolumn,
            'correlation': 'RR,LL',
            'antenna':'*&*', 'field':field,
            'averagedata':True, 'avgchannel':str(nchan), 'avgtime':avgtime,
            'xselfscale':True, 'xsharedaxis': True, 'coloraxis': 'spw',
            'plotrange':[-1,-1,-180,180],
            'plotfile': plot_file, 'expformat': 'png', 'customsymbol': True, 'symbolshape': 'circle',
            'width':1200, 'height':573, 'symbolsize':4,'clearplots':False, 'overwrite':True, 'showgui':showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()

#'#
#'#    plotms(vis=msfile, xaxis='UVwave', yaxis='phase', title='Model Phase vs UVWave {0} (color=spw)'.format(field),
#'#    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=1,
#'#    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
#'#    antenna='*&*', field=field,
#'#    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
#'#    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
#'#    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#    width=1200, height=573, symbolsize=4, clearplots=False, overwrite=True, showgui=showgui)


def make_uvplt(eMCP):
    msinfo = eMCP['msinfo']
    num_proc = eMCP['defaults']['plot_data']['num_proc']
    plots_data_dir = './weblog/plots/plots_uvplt/'
    emutils.makedir(plots_data_dir)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    logger.info('Producing uvplot for: {}'.format(','.join(allsources)))
    # UVplot all sources
    for field in allsources:
        if field in mssources:
            single_uvplt(msinfo, field, plots_data_dir)
        else:
            logger.warning('Cannot plot {0}. Source not in ms.'.format(field))
##    pool = multiprocessing.Pool(num_proc)
##    input_args = [(msinfo, field, plots_data_dir) for field in fields]
##    pool.map(single_uvplt, input_args)
##    pool.close()
##    pool.join()
    # UVplot model, calibrators
    calsources = msinfo['sources']['calsources'].split(',')
    for field in calsources:
        if field in mssources:
            single_uvplt_model(msinfo, field, plots_data_dir)
        else:
            logger.warning('Cannot plot {0}. Source not in ms.'.format(field))
##    pool = multiprocessing.Pool(num_proc)
##    input_args = [(msinfo, field, plots_data_dir) for field in calsources]
##    pool.map(single_uvplt_model, input_args)
##    pool.close()
##    pool.join()
    logger.info('uvplts finished')

#def single_uvcov(field, u, v, freqs, plot_file):
#    c = 299792458.
#    # Color scale:
#    norm = matplotlib.colors.Normalize(vmin=np.min(freqs/1e9),
#                                       vmax=np.max(freqs/1e9))
#    cmap = matplotlib.cm.get_cmap('Spectral')
#    fig = plt.figure()
#    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#    ax1 = fig.add_axes([0.8, 0.1, 0.02, 0.8])
#    for i, freqi in enumerate(freqs):
#        col= cmap(norm(freqi/1e9)) #color
#        fc = freqi/c/1e6
#        ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
#                 color=col, mec=col, label='{0:4.1f}GHz'.format(freqi/1e9))
#        ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
#                 color=col, mec=col)
#    #lgnd = ax.legend(numpoints=1, markerscale=6, frameon=False, ncol=2, prop={'size':8})
#    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap,
#                                norm=norm, orientation='vertical')
#    ax1.yaxis.set_label_position("right")
#    ax1.set_ylabel('Frequency [GHz]')
#    ax.set_xlabel('U [Mlambda]')
#    ax.set_ylabel('V [Mlambda]')
#    ax.set_aspect('equal')
#    ax.set_title(field)
#    main_lim = np.max(np.abs(ax.get_ylim() + ax.get_xlim()))
#    ax.set_xlim(-main_lim, +main_lim)
#    ax.set_ylim(-main_lim, +main_lim)
#    fig.savefig(plot_file, dpi=120, bbox_inches='tight')
#
#def read_uvw(msfile, field):
#    ms.open(msfile)
#    staql={'field':field, 'spw':'0'}
#    ms.msselect(staql)
#    uv = ms.getdata(['u', 'v', 'flag'])
#    ms.close()
#    flags = np.mean(uv['flag'], axis=(0,1)) < 0.9
#    u = uv['u'][flags]
#    v = uv['v'][flags]
#    return u, v

def make_uvcov(msfile, msinfo):
    plots_obs_dir = './weblog/plots/plots_observation/'
    emutils.makedir(plots_obs_dir)
    # CASA 5.4 has a bug, it selects the uv limits of the first spw
    # I create this manual limit as a compromise
    max_freq = emutils.read_keyword(msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW').max()
#'#    msmd.open(msfile)
    c = 299792458.
#'#    max_freq = np.max(np.array([msmd.chanfreqs(spw) for spw in
#'#                                msmd.datadescids()]))
    max_uvdist = 217000.0/c*max_freq  # For 217 km baseline
    #freqs = get_freqs(msfile, allfreqs=True)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    logger.info('Plotting uvcov for:')
    for f in allsources:
        if f in mssources:
            plot_file = os.path.join(plots_obs_dir,'{0}_uvcov_{1}.png'.format(msinfo['msfilename'],f))
            logger.info('{0}'.format(f))
            avgtime = '32'
            nchan = msinfo['nchan']
            commands = {}
            commands['plotms'] = {
                    'vis' :  msfile,
                    'xaxis': 'Uwave',
                    'yaxis': 'Vwave',
                    'field': f,
                    'title': f,
                    'correlation': 'RR',
                    'spw': '',
                    'coloraxis': 'spw',
                    'width': 900,
                    'symbolsize': 1,
                    'plotrange': [-max_uvdist,+max_uvdist,-max_uvdist,+max_uvdist],
                    'averagedata':  True,
                    'avgtime': avgtime,
                    'avgchannel': str(int(nchan/8)),
                    'plotfile':  plot_file,
                    'expformat':  'png',
                    'customsymbol':  True,
                    'symbolshape':  'circle',
                    'overwrite': True,
                    'showlegend': True,
                    'showgui': False}
            em.run_casa_command(commands, 'plotms')
            em.find_casa_problems()
#'#            plotms(vis=msfile, xaxis='Uwave', yaxis='Vwave', field=f, title=f,
#'#            correlation = 'RR', spw='', coloraxis = 'spw',
#'#            width=900, height=900, symbolsize=1,
#'#            plotrange=[-max_uvdist,+max_uvdist,-max_uvdist,+max_uvdist],
#'#            averagedata = True, avgtime=avgtime, avgchannel = str(int(nchan/8)),
#'#            plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
#'#            overwrite=True, showlegend=False, showgui=False)
#            u, v = read_uvw(msfile, f)
#            single_uvcov(f, u, v, freqs, plot_file)
        else:
            logger.info('Cannot plot uvcov for {0}. Source not in ms.'.format(f))


def make_elevation(msfile, msinfo):
    plots_obs_dir = './weblog/plots/plots_observation/'
    emutils.makedir(plots_obs_dir)
    plot_file = plots_obs_dir+'{0}_elevation.png'.format(msinfo['msfilename'])
    logger.info('Plotting elevation to:')
    logger.info('{}'.format(plot_file))
    avgtime = '16'
    showgui = False
    commands = {}
    commands['plotms'] = {
            'vis' :  msfile,
            'xaxis': 'time',
            'yaxis': 'elevation',
            'correlation': 'RR',
            'spw': '',
            'coloraxis': 'field',
            'width': 900,
            'symbolsize': 5,
            'plotrange': [-1,-1,0,90],
            'averagedata':  True,
            'avgtime': avgtime,
            'plotfile':  plot_file,
            'expformat':  'png',
            'customsymbol':  True,
            'symbolshape':  'circle',
            'overwrite': True,
            'showlegend': True,
            'showgui': showgui}
    em.run_casa_command(commands, 'plotms')
    em.find_casa_problems()


### Flag statistics
def fperc(x):
    return 1.0*x['flagged']/x['total']

def sort_list(item, flagged, list_order):
    order = {a:i for i, a in enumerate(list_order)}
    item_sorted, flagged_sorted = np.asarray(sorted(zip(item,flagged),key=lambda d:order[d[0]])).T
    return item_sorted, np.asfarray(flagged_sorted)

def read_scan_summary(datain):
    ms.open(datain)
    scan_summary = ms.getscansummary()
    ms.close()
    return scan_summary

def count_flags(flag_stats, label, list_order=[]):
    item = []
    flagged = []
    for s in flag_stats[label].keys():
        flagged.append(fperc(flag_stats[label][s]))
        try:
            item.append(int(s))
        except:
            item.append(s)
    if len(list_order) == 0:
        order = np.array(item).argsort()
        item_sorted = np.array(item)[order]
        flagged_sorted = np.array(flagged)[order]
    else:
        item_sorted,flagged_sorted = sort_list(item, flagged, list_order=list_order)
    return item_sorted,flagged_sorted


def plot_flagstatistics(flag_stats, msinfo, step):
    # Different colors for each field
#'#    scan_summary = read_scan_summary(msinfo['msfile'])
#'#    scan_fieldID_dict =  {si:scan_summary[str(si)]['0']['FieldId'] for si in scan_summary.keys()}
#'#    #vis_fields = vishead(msinfo['msfile'],mode='list',listitems='field')['field'][0]
#'#    msmd.open(msinfo['msfile'])
#'#    vis_fields = np.array(msmd.fieldnames())
#'#    msmd.done()
    msfile = msinfo['msfile']
    scan_number = emutils.read_keyword(msfile, 'SCAN_NUMBER')
    field_id = emutils.read_keyword(msfile, 'FIELD_ID')
    scan_fieldID_dict = {}
    for scan in np.unique(scan_number):
        scan_fieldID_dict[str(scan)] = np.unique(field_id[np.where(scan_number == scan)[0]])[0]
    vis_fields = emutils.read_keyword(msfile, 'NAME', 'FIELD')


    # Compute % statistics
    i_scan, f_scan = count_flags(flag_stats, 'scan')
    i_field, f_field = count_flags(flag_stats, 'field', list_order=vis_fields)
    i_corr, f_corr = count_flags(flag_stats, 'correlation', list_order=['RR','LL','RL','LR'])
    i_spw, f_spw = count_flags(flag_stats, 'spw')
    i_ant, f_ant = count_flags(flag_stats, 'antenna', list_order=msinfo['antennas'])

    fig = plt.figure(figsize=(25,4))
    plt.subplots_adjust(wspace=0.01)
    ax1 = fig.add_subplot(1,5,(1,2))
    ax2 = fig.add_subplot(153)
#    ax3 = fig.add_subplot(153)
    ax4 = fig.add_subplot(154, sharey=ax2)
    ax5 = fig.add_subplot(155, sharey=ax2)

    scan_fieldID = np.array([scan_fieldID_dict[str(si)] for si in i_scan])
    for i, fi in enumerate(i_field):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), zorder=10)
        field_value = f_field[np.argwhere(i_field==fi)[0][0]]
        ax2.bar(i, field_value, alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), align='center', zorder=10)
        ax2.text(i-0.1, 0.9*field_value, "{0:2.0f}".format(field_value*100.),
                 color='k', va='center', zorder=12)

#    ax3.bar(range(len(i_corr)), f_corr, alpha=0.5, color='k', width=1, align='center')
    ax4.bar(range(len(i_spw)), f_spw, alpha=1.0, color='0.5', width=1,
            align='center', zorder=10)
    ax5.bar(range(len(i_ant)), f_ant, alpha=1.0, color='0.5', width=1,
            align='center', zorder=10)

    ax2.axes.set_xticks(range(len(i_field)))
#    ax3.axes.set_xticks(range(len(i_corr)))
    ax5.axes.set_xticks(range(len(i_ant)))
    ax2.axes.set_xticks(range(len(i_field)))
    ax2.axes.set_xticklabels([])
#    ax3.axes.set_xticks(range(len(i_corr)))
#    ax3.axes.set_xticklabels(i_corr)
    ax4.axes.set_xticks(range(len(i_spw)))
    ax4.axes.set_xticklabels(range(len(i_spw)))
    ax5.axes.set_xticks(range(len(i_ant)))
    ax5.axes.set_xticklabels(i_ant)
    ax2.axes.set_yticklabels([])
    ax4.axes.set_yticklabels([])
    ax5.axes.set_yticklabels([])

    [ax2.annotate('{0} ({1})'.format(si, i), (i+0.1, 0.95), va='top', ha='right', rotation=90, zorder=100) for i, si in enumerate(i_field)]
    #[ax3.annotate(si, (i+0.5, f_corr[i])) for i, si in enumerate(i_corr)]
    #[ax5.annotate(si, (i+0.3, f_ant[i])) for i, si in enumerate(i_ant)]

    for i, v in enumerate(f_spw):
        ax4.text(i-0.1, 0.9*v,"{0:2.0f}".format(v*100.), color='k',
                 va='center', zorder=12)
    for i, v in enumerate(f_ant):
        ax5.text(i-0.1, 0.9*v,"{0:2.0f}".format(v*100.), color='k',
                 va='center', zorder=12)

    ax1.set_title('Scan')
    ax2.set_title('Field')
#    ax3.set_title('Correlation')
    ax4.set_title('spw')
    ax5.set_title('Antenna')

    ax1.set_xlabel('Scan')
    ax4.set_xlabel('spw')
    ax1.xaxis.set_major_locator(MultipleLocator(10))
    ax1.set_ylabel('Flagged fraction')
    #ax2.set_ylabel('Flagged fraction')
    #ax4.set_ylabel('Flagged fraction')

    ax1.grid(axis='y', zorder = -1000, ls='-', color='0.6')
    ax2.grid(axis='y', zorder = -1000, ls='-', color='0.6')
    ax4.grid(axis='y', zorder = -1000, ls='-', color='0.6')
    ax5.grid(axis='y', zorder = -1000, ls='-', color='0.6')

    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)
#    ax3.set_ylim(0,1)
    ax4.set_ylim(0,1)
    ax5.set_ylim(0,1)

#'#    ax1.set_xlim(np.min(i_scan)-0.5, np.max(i_scan)+0.5)
    ax1.set_xlim(np.min(i_scan)-1.0, np.max(i_scan)+0.)
    ax2.set_xlim(-0.5, len(i_field)-0.5)
#    ax3.set_xlim(-0.5, len(i_corr)-0.5)
    ax4.set_xlim(-0.5, len(i_spw)-0.5)
    ax5.set_xlim(-0.5, len(i_ant)-0.5)

    #ax1.legend(loc=0)
    #ax2.legend(loc=0)

    plots_obs_dir = './weblog/plots/plots_flagstats/'
    plot_file1 = plots_obs_dir+'{0}_flagstats_{1}.png'.format(msinfo['msfilename'], step)
    #logger.info('Plot flagstats: {0}'.format(plot_file1))
    fig.savefig(plot_file1, bbox_inches='tight')

    # Plot only scans:
    fig = plt.figure(figsize=(50,8))
    ax1 = fig.add_subplot(111)

    for i, fi in enumerate(np.unique(i_field)):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), zorder=10)

    try:
        ax1.legend(loc=0)
    except:
        pass
    ax1.grid(axis='y', zorder = -1000, ls='-', color='0.6')

    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xlim(np.min(i_scan)-0.5, np.max(i_scan)+0.5)
    ax1.set_ylim(0,1)
    ax1.set_xlabel('Scan number')
    ax1.set_ylabel('Flagged fraction')

    plot_file2 = plots_obs_dir+'{0}_flagstats_scans_{1}.png'.format(msinfo['msfilename'],
                                                       step)
    #logger.info('Plot flagstats scans: {0}'.format(plot_file2))
    fig.savefig(plot_file2, bbox_inches='tight')


def plot_Lo_drops(phscal_scans, scans, amp_mean, lo_dropout_scans, phscal, eMCP):
    plots_obs_dir = './weblog/plots/plots_flagstats/'
    emutils.makedir(plots_obs_dir)
    msinfo = eMCP['msinfo']
    drops = np.array([scan in lo_dropout_scans for scan in phscal_scans])
    fig = plt.figure(figsize=(30,8))
    ax1 = fig.add_subplot(111)

    ax1.bar(scans-0.5, np.ones_like(scans)*np.max(amp_mean)*1.2, alpha=0.2,
            color='0.5', width=1),
    ax1.bar(phscal_scans-0.5, amp_mean, alpha=1.0,
            color='0.5', width=1,
            label='{0}'.format(phscal))
    if lo_dropout_scans != []:
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

def read_calfluxes(calfluxes, k, eMfactor):
    freq = calfluxes['freq']
    spws = calfluxes['spwID']
    fieldName = calfluxes[k]['fieldName']
    spindex = calfluxes[k]['spidx'][1]
    espindex = calfluxes[k]['spidxerr'][1]
    S0 = calfluxes[k]['fitFluxd']*eMfactor
    eS0 = calfluxes[k]['fitFluxdErr']*eMfactor
    freq0 = calfluxes[k]['fitRefFreq']
    flux = np.ones(len(spws))*np.nan
    eflux = np.ones(len(spws))*np.nan
    for i, spw in enumerate(spws):
        try:
            flux[i] = calfluxes[k][str(spw)]['fluxd'][0]
            eflux[i] = calfluxes[k][str(spw)]['fluxdErr'][0]
            if eflux[i] <= 0.0:
                flux[i] = np.nan
                eflux[i] = np.nan
                freq[i] = np.nan
        except:
            pass
    flux *= eMfactor
    eflux *= eMfactor
    return freq,spws,fieldName,spindex,espindex,S0,eS0,freq0,flux,eflux

def fluxscale_models(calfluxes, eMfactor, msinfo):
    factor_unit= 1e-9
    units = 'Jy'
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(111)

    freq_min, freq_max = 1e20,0.
    for k in calfluxes.keys():
        if type(calfluxes[k]) is dict:
            freq,spws,fieldName,spindex,espindex,S0,eS0,freq0,flux,eflux = read_calfluxes(calfluxes,k, eMfactor)
            freq_min = np.nanmin([freq_min, np.nanmin(freq)])
            freq_max = np.nanmax([freq_max, np.nanmax(freq)])
            freqspan = np.nanmax(freq) - np.nanmin(freq)
            freqlim = np.array([np.nanmin(freq) - 0.025*freqspan, np.nanmax(freq) + 0.1*freqspan])
            fluxfit = S0*(freqlim/freq0)**spindex
            ff_min  = (S0-eS0)*(freqlim/freq0)**(spindex-espindex)
            ff_max  = (S0+eS0)*(freqlim/freq0)**(spindex+espindex)
            ff_min2 = (S0-eS0)*(freqlim/freq0)**(spindex+espindex)
            ff_max2 = (S0+eS0)*(freqlim/freq0)**(spindex-espindex)
            label = '{0:>9s}: Flux density = {1:6.3f} +/-{2:6.3f}, '\
                    'spidx ={3:5.2f}+/-{4:5.2f}'.format(fieldName,
                        calfluxes[k]['fitFluxd']*eMfactor, calfluxes[k]['fitFluxdErr']*eMfactor,
                        calfluxes[k]['spidx'][1], calfluxes[k]['spidxerr'][1])
            p, = ax1.plot(freqlim*factor_unit, fluxfit, '-', linewidth = 2,  zorder = -5,
                     label=label)
            color1 = str(p.get_color())
            #ax1.errorbar(freq*factor_unit, flux, eflux, fmt = 'o', color =color1, mec = color1, zorder = 10)
            ax1.plot(freq*factor_unit, flux,  marker = 'o', ls='', color ='k',mec= 'k', zorder = 10)
            ax1.fill_between(freqlim*factor_unit, ff_min, ff_max,
                             facecolor=color1,color=color1, alpha = 1.0,linewidth = 0, zorder = -32)
            ax1.fill_between(freqlim*factor_unit, ff_max2,ff_min2,
                             facecolor=color1,color=color1, alpha = 1.0,linewidth = 0, zorder = -32)

    freq = calfluxes['freq']
    freqspan = freq_max - freq_min
    freqlim = np.array([freq_min - 0.025*freqspan, freq_max + 0.1*freqspan])
    ax1.set_xlabel("Frequency [GHz]")
    ax1.set_ylabel("Flux density [Jy]")
    ax1.set_xlim(freqlim[0]*factor_unit, freqlim[1]*factor_unit)
    #ax1.set_ylim(0.0, ax1.get_ylim()[1])
    ax1.xaxis.set_major_locator(MaxNLocator(5))
    ax1.xaxis.set_minor_locator(MaxNLocator(20))
    #ax1.yaxis.set_major_locator(MaxNLocator(6))
    #ax1.yaxis.set_minor_locator(MaxNLocator(6*5))
    leg = ax1.legend(loc=0, borderaxespad=1)
    ax1.grid(ls='-', alpha=0.3, zorder=-50)
    #ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.yaxis.set_major_formatter(ScalarFormatter())

    plots_obs_dir = calib_dir
    plot_file = plots_obs_dir+'{0}_fluxscale.png'.format(msinfo['msfilename'])
    fig.savefig(plot_file, bbox_inches='tight')
