#!/usr/local/python
import os
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
#import multiprocessing
from matplotlib.ticker import MultipleLocator
import datetime

plt.ioff()

import functions.weblog as emwlog

# CASA imports
from taskinit import *
from tasks import *
from casac import casac
msmd = casac.msmetadata()

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

line0 = '-'*10

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


def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass

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

def get_scans(msfile, field):
    ms.open(msfile)
    ms.msselect({'field':field})
    scans = ms.getdata('scan_number')['scan_number']
    num_scans = len(np.unique(scans))
    ms.close()
    return num_scans, scans

def get_freqs(msfile, allfreqs=False):
    ms.open(msfile)
    axis_info = ms.getdata(['axis_info'],ifraxis=True)
    ms.close()
    if allfreqs:
        channels = axis_info['axis_info']['freq_axis']['chan_freq']
        freqs=np.sort(axis_info['axis_info']['freq_axis']['chan_freq'].flatten())[::len(channels)/4]
    else:
        freqs = axis_info['axis_info']['freq_axis']['chan_freq'].mean(axis=0)
    return freqs

def count_active_baselines(msfile):
    msmd.open(msfile)
    a = np.array(msmd.antennanames())
    b = np.array(msmd.baselines())
    msmd.done()
    active_baselines = 0
    for i, ant1 in enumerate(a):
        for ant2 in a[i+1:]:
            j = np.argwhere(ant2==a)[0][0]
            if b[i,j]:
                active_baselines += 1
    return active_baselines


def single_4plot(msinfo, field, datacolumn, plots_data_dir):
    logger.info('Visibility plots for field: {0}, datacolumn: {1}'.format(field, datacolumn))
    msfile = msinfo['msfile']
    tb.open(msfile+'/SPECTRAL_WINDOW')
    nchan = str(tb.getcol('NUM_CHAN')[0])
    tb.close()
    plot_file = plots_data_dir+'{0}_4plot_{1}_{2}'.format(msinfo['msfilename'],
                                                          field, datacolumn)
    num_baselines = count_active_baselines(msfile)
    gridrows = num_baselines
    gridcols = 1
    showgui=False
    avgtime = '300'
    baseline = '*&*'
    w, h = 1200, num_baselines*200
    plotms(vis=msfile, xaxis='time', yaxis='amp', title='Amp vs Time {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgchannel = nchan, avgtime='4',
    xselfscale = True, xsharedaxis = True, coloraxis = 'spw', plotrange=[-1,-1,0,-1],
    plotfile = plot_file+'0.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='time', yaxis='phase', title='Phase vs Time {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgchannel = nchan, avgtime='4',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
    plotfile = plot_file+'1.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='freq', yaxis='amp', title='Amp vs Frequency {0} (color=corr)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgtime=avgtime, avgchannel='4',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,0,-1],
    plotfile = plot_file+'2.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='freq', yaxis='phase', title='Phase vs Frequency {0} (color=corr)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgtime=avgtime, avgchannel='4',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,-180,180],
    plotfile = plot_file+'3.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)
    logger.info('Finished {0}: {1}{{0-4}}.png'.format(field, plot_file))


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
    makedir(plots_data_dir)
    fields = msinfo['sources']['allsources'].split(',')
    logger.info('Producing visibility plots for: {}'.format(msinfo['sources']['allsources']))
    for field in fields:
        single_4plot(msinfo, field, datacolumn, plots_data_dir)
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
    plotms(vis=msfile, xaxis='UVwave', yaxis='amp', title='Amplitude vs UVWave {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
    antenna='*&*', field=field,
    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw',
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    symbolsize=4, clearplots=True, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='UVwave', yaxis='phase', title='Phase vs UVWave {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=1,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
    antenna='*&*', field=field,
    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=1200, height=573, symbolsize=4, clearplots=False, overwrite=True, showgui=showgui)

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
    plotms(vis=msfile, xaxis='UVwave', yaxis='amp', title='Model Amplitude vs UVWave {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
    antenna='*&*', field=field,
    averagedata = True, avgtime=avgtime, avgchannel = str(int(nchan/16)),
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw',
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    symbolsize=4, clearplots=True, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='UVwave', yaxis='phase', title='Model Phase vs UVWave {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=1,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
    antenna='*&*', field=field,
    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=1200, height=573, symbolsize=4, clearplots=False, overwrite=True, showgui=showgui)


def make_uvplt(eMCP):
    msinfo = eMCP['msinfo']
    num_proc = eMCP['defaults']['plot_data']['num_proc']
    plots_data_dir = './weblog/plots/plots_uvplt/'
    makedir(plots_data_dir)
    fields = msinfo['sources']['allsources'].split(',')
    logger.info('Producing uvplot for: {}'.format(msinfo['sources']['allsources']))
    # UVplot all sources
    for field in fields:
        single_uvplt(msinfo, field, plots_data_dir)
##    pool = multiprocessing.Pool(num_proc)
##    input_args = [(msinfo, field, plots_data_dir) for field in fields]
##    pool.map(single_uvplt, input_args)
##    pool.close()
##    pool.join()
    # UVplot model, calibrators
    calsources = msinfo['sources']['calsources'].split(',')
    for field in calsources:
        single_uvplt_model(msinfo, field, plots_data_dir)
##    pool = multiprocessing.Pool(num_proc)
##    input_args = [(msinfo, field, plots_data_dir) for field in calsources]
##    pool.map(single_uvplt_model, input_args)
##    pool.close()
##    pool.join()
    logger.info('uvplts finished')

def single_uvcov(field, u, v, freqs, plot_file):
    c = 299792458.
    # Color scale:
    norm = matplotlib.colors.Normalize(vmin=np.min(freqs/1e9),
                                       vmax=np.max(freqs/1e9))
    cmap = matplotlib.cm.get_cmap('Spectral')
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax1 = fig.add_axes([0.8, 0.1, 0.02, 0.8])
    for i, freqi in enumerate(freqs):
        col= cmap(norm(freqi/1e9)) #color
        fc = freqi/c/1e6
        ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
                 color=col, mec=col, label='{0:4.1f}GHz'.format(freqi/1e9))
        ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
                 color=col, mec=col)
    #lgnd = ax.legend(numpoints=1, markerscale=6, frameon=False, ncol=2, prop={'size':8})
    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm, orientation='vertical')
    ax1.yaxis.set_label_position("right")
    ax1.set_ylabel('Frequency [GHz]')
    ax.set_xlabel('V [Mlambda]')
    ax.set_ylabel('U [Mlambda]')
    ax.set_aspect('equal')
    ax.set_title(field)
    main_lim = np.max(np.abs(ax.get_ylim() + ax.get_xlim()))
    ax.set_xlim(-main_lim, +main_lim)
    ax.set_ylim(-main_lim, +main_lim)
    fig.savefig(plot_file, dpi=120, bbox_inches='tight')

def read_uvw(msfile, field):
    ms.open(msfile)
    staql={'field':field, 'spw':'0'}
    ms.msselect(staql)
    uv = ms.getdata(['u', 'v'])
    ms.close()
    u = uv['u']
    v = uv['v']
    return u, v

def make_uvcov(msfile, msinfo):
    logger.info('Plotting uv-coverage for all sources'.format())
    plots_obs_dir = './weblog/plots/plots_observation/'
    makedir(plots_obs_dir)
    freqs = get_freqs(msfile, allfreqs=True)
    fields_ms = msinfo['sources']['allsources'].split(',')
    logger.info('Plotting uvcov for:')
    for f in fields_ms:
        if f in fields_ms:
            plot_file = plots_obs_dir+'{0}_uvcov_{1}.png'.format(msinfo['msfilename'],f)
            logger.info('{0}'.format(f))
            u, v = read_uvw(msfile, f)
            single_uvcov(f, u, v, freqs, plot_file)
        else:
            logger.info('Cannot plot uvcov for {0}. Source not in ms.'.format(f))


def make_elevation(msfile, msinfo):
    plots_obs_dir = './weblog/plots/plots_observation/'
    makedir(plots_obs_dir)
    plot_file = plots_obs_dir+'{0}_elevation.png'.format(msinfo['msfilename'])
    logger.info('Plotting elevation to:')
    logger.info('{}'.format(plot_file))
    avgtime = '32'
    showgui = False
    plotms(vis=msfile, xaxis='time', yaxis='elevation', correlation = 'RR',
    spw='0', coloraxis = 'field', width=900, symbolsize=5, plotrange=[-1,-1,0,90],
    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,showlegend=True, showgui=showgui)



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
    scan_summary = read_scan_summary(msinfo['msfile'])
    scan_fieldID_dict =  {si:scan_summary[str(si)]['0']['FieldId'] for si in scan_summary.keys()}
    #vis_fields = vishead(msinfo['msfile'],mode='list',listitems='field')['field'][0]
    msmd.open(msinfo['msfile'])
    vis_fields = np.array(msmd.fieldnames())
    msmd.done()

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
    for i, fi in enumerate(vis_fields):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), zorder=10)
        field_value = f_field[np.argwhere(i_field==fi)[0][0]]
        ax2.bar(i, field_value, alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), align='center', zorder=10)
        ax2.text(i-0.1, 0.9*field_value, "{0:2.0f}".format(field_value*100.), color='k', va='center')

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

    [ax2.annotate('{0} ({1})'.format(si, i), (i+0.1, 0.95), va='top', ha='right', rotation=90) for i, si in enumerate(i_field)]
    #[ax3.annotate(si, (i+0.5, f_corr[i])) for i, si in enumerate(i_corr)]
    #[ax5.annotate(si, (i+0.3, f_ant[i])) for i, si in enumerate(i_ant)]

    for i, v in enumerate(f_spw):
        ax4.text(i-0.1, 0.9*v,"{0:2.0f}".format(v*100.), color='k', va='center')
    for i, v in enumerate(f_ant):
        ax5.text(i-0.1, 0.9*v,"{0:2.0f}".format(v*100.), color='k', va='center')

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

    ax1.set_xlim(np.min(i_scan)-0.5, np.max(i_scan)+0.5)
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

    for i, fi in enumerate(np.unique(vis_fields[scan_fieldID])):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), zorder=10)

    ax1.legend(loc=0)
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
    makedir(plots_obs_dir)
    msinfo = eMCP['msinfo']
    drops = np.array([scan in lo_dropout_scans for scan in phscal_scans])
    fig = plt.figure(figsize=(30,8))
    ax1 = fig.add_subplot(111)

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

