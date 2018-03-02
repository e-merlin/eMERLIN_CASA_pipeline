#!/usr/local/python
import os
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
from matplotlib.ticker import MultipleLocator

plt.ioff()

# CASA imports
from taskinit import *
from tasks import *

import logging
logger = logging.getLogger('logger')


# Functions to save and load dictionaries
def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
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


def single_4plot_old(msfile, field, baseline, datacolumn, plot_file):
    logger.info('Baseline: {0:6s}, Output: {1}'.format(baseline, plot_file))
    tb.open(msfile+'/SPECTRAL_WINDOW')
    nchan = str(tb.getcol('NUM_CHAN')[0])
    tb.close()
    gridrows = gridcols = 2
    showgui=False
    avgtime = '300'

    plotms(vis=msfile, xaxis='time', yaxis='amp', title='Amp vs Time {0} {1} (color=spw)'.format(field, baseline),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis = 'spw', plotrange=[-1,-1,0,-1],
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4)

    plotms(vis=msfile, xaxis='time', yaxis='phase', title='Phase vs Time {0} {1} (color=spw)'.format(field, baseline),
    gridrows=gridrows, gridcols=gridcols, rowindex=1, colindex=0, plotindex=1,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4, clearplots=False)

    plotms(vis=msfile, xaxis='freq', yaxis='amp', title='Amp vs Frequency {0} {1} (color=corr)'.format(field, baseline),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=2,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,0,-1],
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4, clearplots=False)

    plotms(vis=msfile, xaxis='freq', yaxis='phase', title='Phase vs Frequency {0} {1} (color=corr)'.format(field, baseline),
    gridrows=gridrows, gridcols=gridcols, rowindex=1, colindex=1, plotindex=3,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,-180,180],
    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=2600, height=1200, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)


def make_4plots_old(msfile, msinfo, datacolumn='data'):
    if datacolumn == 'data':
        plots_data_dir = './plots/plots_data/'
    elif datacolumn == 'corrected':
        plots_data_dir = './plots/plots_corrected/'
    else:
        plots_data_dir = './plots/'
    makedir(plots_data_dir)
    for f in msinfo['sources']['allsources'].split(','):
        logger.info('Generating plot for msfile: {0}, field: {1}, column: {2}'.format(
                                                    (msfile), f, datacolumn))
        for baseline in msinfo['baselines']:
            plot_file = plots_data_dir+'{0}_4plot_{1}_{2}_{3}.png'.format(
                                                         msinfo['msfilename'], f,
                                                         baseline.replace('&','-'),
                                                         datacolumn)
            single_4plot_old(msfile, f, baseline, datacolumn, plot_file)




def single_4plot((msinfo, field, datacolumn, plots_data_dir)):
    logger.info('Visibility plots for field: {0}, datacolumn: {1}'.format(field, datacolumn))
    msfile = msinfo['msfile']
    tb.open(msfile+'/SPECTRAL_WINDOW')
    nchan = str(tb.getcol('NUM_CHAN')[0])
    tb.close()
    plot_file = plots_data_dir+'{0}_4plot_{1}_{2}'.format(msinfo['msfilename'],
                                                          field, datacolumn)
    num_baselines = len(msinfo['baselines'])
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
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis = 'spw', plotrange=[-1,-1,0,-1],
    plotfile = plot_file+'0.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='time', yaxis='phase', title='Phase vs Time {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', plotrange=[-1,-1,-180,180],
    plotfile = plot_file+'1.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='freq', yaxis='amp', title='Amp vs Frequency {0} (color=corr)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,0,-1],
    plotfile = plot_file+'2.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)

    plotms(vis=msfile, xaxis='freq', yaxis='phase', title='Phase vs Frequency {0} (color=corr)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field, iteraxis='baseline',
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', plotrange=[-1,-1,-180,180],
    plotfile = plot_file+'3.png', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=w, height=h, symbolsize=4,clearplots=False, overwrite=True, showgui=showgui)
    logger.info('Finished field: {0}, output: {1}{{0-4}}.png'.format(field, plot_file))

#def make_4plots((msinfo, field, datacolumn, plots_data_dir)):
#    msfile = msinfo['msfile']
#    plot_file = plots_data_dir+'{0}_4plot_{1}_{2}'.format(msinfo['msfilename'],
#                                                          field, datacolumn)
#    if datacolumn == 'data':
#        plots_data_dir = './plots/plots_data/'
#    elif datacolumn == 'corrected':
#        plots_data_dir = './plots/plots_corrected/'
#    else:
#        plots_data_dir = './plots/'
#    makedir(plots_data_dir)
#    for f in msinfo['sources']['mssources'].split(','):
#        if f != '':
#            logger.info('Generating plot for msfile: {0}, field: {1}, column: {2}'.format(
#                                                    (msfile), field, datacolumn))
#            plot_file = plots_data_dir+'{0}_4plot_{1}_{2}'.format(
#                                                         msinfo['msfilename'],
#                                                         field, datacolumn)
#            num_baselines = len(msinfo['baselines'])
#            single_4plot(msfile, field, num_baselines, datacolumn, plot_file)

def make_4plots(msfile, msinfo, datacolumn='data'):
    logger.info('Start plot_{}'.format(datacolumn))
    if datacolumn == 'data':
        plots_data_dir = './plots/plots_data/'
    elif datacolumn == 'corrected':
        plots_data_dir = './plots/plots_corrected/'
    else:
        plots_data_dir = './plots/'
    makedir(plots_data_dir)
    fields = msinfo['sources']['allsources'].split(',')
    logger.info('Producing in parallel visibility plots for fields: {}'.format(msinfo['sources']['allsources']))
    num_proc = 1
    pool = multiprocessing.Pool(num_proc)
    input_args = [(msinfo, field, datacolumn, plots_data_dir) for field in fields]
    p = pool.map(single_4plot, input_args)
    pool.close()
    pool.join()
    logger.info('Visibility plots finished')
    logger.info('End plot_{}'.format(datacolumn))

def single_uvplt((msinfo, field, plots_data_dir)):
    logger.info('uvplt for field: {}'.format(field))
    plot_file =  plots_data_dir+'{0}_uvplt_{1}.png'.format(msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn='corrected'
    avgtime = ''
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

def single_uvplt_model((msinfo, field, plots_data_dir)):
    logger.info('uvplt (model) for field: {}'.format(field))
    plot_file =  plots_data_dir+'{0}_uvpltmodel_{1}.png'.format(msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn='model'
    avgtime = ''
    showgui = False
    gridrows = 1
    gridcols = 2
    plotms(vis=msfile, xaxis='UVwave', yaxis='amp', title='Model Amplitude vs UVWave {0} (color=spw)'.format(field),
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR,LL',
    antenna='*&*', field=field,
    averagedata = True, avgtime=avgtime, avgchannel = str(nchan),
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


def make_uvplt(msinfo):
    num_proc = 1
    plots_data_dir = './plots/plots_uvplt/'
    makedir(plots_data_dir)
    fields = msinfo['sources']['allsources'].split(',')
    logger.info('Producing in parallel uvplot for fields: {}'.format(msinfo['sources']['allsources']))
    # UVplot all sources
    pool = multiprocessing.Pool(num_proc)
    input_args = [(msinfo, field, plots_data_dir) for field in fields]
    pool.map(single_uvplt, input_args)
    pool.close()
    pool.join()
    # UVplot model, calibrators
    calsources = msinfo['sources']['calsources'].split(',')
    pool = multiprocessing.Pool(num_proc)
    input_args = [(msinfo, field, plots_data_dir) for field in fields]
    pool.map(single_uvplt_model, input_args)
    pool.close()
    pool.join()
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
    plots_obs_dir = './plots/plots_observation/'
    makedir(plots_obs_dir)
    freqs = get_freqs(msfile, allfreqs=True)
    tb.open(msfile+'/FIELD')
    fields_ms = tb.getcol('NAME')
    tb.close()
    #for f in msinfo['sources']['allsources'].split(','):
    for f in fields_ms:
        if f in fields_ms:
            plot_file = plots_obs_dir+'{0}_uvcov_{1}.png'.format(msinfo['msfilename'],f)
            logger.info('Plotting uvcov for {0}: {1}'.format(f, plot_file))
            u, v = read_uvw(msfile, f)
            single_uvcov(f, u, v, freqs, plot_file)
        else:
            logger.info('Cannot plot uvcov for {0}. Source not in ms.'.format(f))


def make_elevation(msfile, msinfo):
    plots_obs_dir = './plots/plots_observation/'
    makedir(plots_obs_dir)
    plot_file = plots_obs_dir+'{0}_elevation.png'.format(msinfo['msfilename'])
    logger.info('Plotting elevation to: {}'.format(plot_file))
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


def plot_flagstatistics(flag_stats, msinfo):
    # Different colors for each field
    scan_summary = read_scan_summary(msinfo['msfile'])
    scan_fieldID_dict =  {si:scan_summary[str(si)]['0']['FieldId'] for si in scan_summary.keys()}
    vis_fields = vishead(msinfo['msfile'],mode='list',listitems='field')['field'][0]

    # Compute % statistics
    i_scan, f_scan = count_flags(flag_stats, 'scan')
    i_field, f_field = count_flags(flag_stats, 'field', list_order=vis_fields)
    i_corr, f_corr = count_flags(flag_stats, 'correlation', list_order=['RR','LL','RL','LR'])
    i_spw, f_spw = count_flags(flag_stats, 'spw')
    i_ant, f_ant = count_flags(flag_stats, 'antenna', list_order=msinfo['antennas'])

    fig = plt.figure(figsize=(12,14))
    plt.subplots_adjust(hspace=0.3)
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(323)
    ax3 = fig.add_subplot(324)
    ax4 = fig.add_subplot(325)
    ax5 = fig.add_subplot(326)

    scan_fieldID = np.array([scan_fieldID_dict[str(si)] for si in i_scan])
    for i, fi in enumerate(vis_fields):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i))
        ax2.bar(i, f_field[np.argwhere(i_field==fi)[0][0]], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i), align='center')

    ax3.bar(range(len(i_corr)), f_corr, alpha=0.5, color='k', width=1, align='center')
    ax4.bar(range(len(i_spw)), f_spw, alpha=0.5, color='k', width=1, align='center')
    ax5.bar(range(len(i_ant)), f_ant, alpha=0.5, color='k', width=1, align='center')

    ax2.axes.set_xticks(range(len(i_field)))
    ax3.axes.set_xticks(range(len(i_corr)))
    ax5.axes.set_xticks(range(len(i_ant)))
    ax2.axes.set_xticks(range(len(i_field)))
    ax2.axes.set_xticklabels([])
    ax3.axes.set_xticks(range(len(i_corr)))
    ax3.axes.set_xticklabels(i_corr)
    ax4.axes.set_xticks(range(len(i_spw)))
    ax4.axes.set_xticklabels(range(len(i_spw)))
    ax5.axes.set_xticks(range(len(i_ant)))
    ax5.axes.set_xticklabels(i_ant)

    [ax2.annotate('{0} ({1})'.format(si, i), (i+0.1, 0.95), va='top', ha='right', rotation=90) for i, si in enumerate(i_field)]
    #[ax3.annotate(si, (i+0.5, f_corr[i])) for i, si in enumerate(i_corr)]
    #[ax5.annotate(si, (i+0.3, f_ant[i])) for i, si in enumerate(i_ant)]

    ax1.set_title('Scan')
    ax2.set_title('Field')
    ax3.set_title('Correlation')
    ax4.set_title('spw')
    ax5.set_title('Antenna')

    ax1.set_xlabel('Scan')
    ax4.set_xlabel('spw')
    ax1.xaxis.set_major_locator(MultipleLocator(10))
    ax1.set_ylabel('Flagged fraction')
    ax2.set_ylabel('Flagged fraction')
    ax4.set_ylabel('Flagged fraction')

    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)
    ax3.set_ylim(0,1)
    ax4.set_ylim(0,1)
    ax5.set_ylim(0,1)

    ax1.set_xlim(0.5, len(i_scan)+0.5)
    ax2.set_xlim(-0.5, len(i_field)-0.5)
    ax3.set_xlim(-0.5, len(i_corr)-0.5)
    ax4.set_xlim(-0.5, len(i_spw)-0.5)
    ax5.set_xlim(-0.5, len(i_ant)-0.5)

    #ax1.legend(loc=0)
    #ax2.legend(loc=0)

    plots_obs_dir = './plots/plots_flagstats/'
    plot_file1 = plots_obs_dir+'{0}_flagstats_all.png'.format(msinfo['msfilename'])
    logger.info('Plotting flagstats all: {0}'.format(plot_file1))
    fig.savefig(plot_file1, bbox_inches='tight')

    # Plot only scans:
    fig = plt.figure(figsize=(50,8))
    ax1 = fig.add_subplot(111)

    for i, fi in enumerate(np.unique(vis_fields[scan_fieldID])):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond]-0.5, f_scan[cond], alpha=1.0,
                color=plt.cm.Set1(1.0*i/len(i_field)), width=1,
                label='{0} ({1})'.format(fi, i))

    ax1.legend(loc=0)
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xlim(0.5, len(i_scan)+0.5)
    ax1.set_ylim(0,1)
    ax1.set_xlabel('Scan number')
    ax1.set_ylabel('Flagged fraction')

    plot_file2 = plots_obs_dir+'{0}_flagstats_scans.png'.format(msinfo['msfilename'])
    logger.info('Plotting flagstats scans: {0}'.format(plot_file2))
    fig.savefig(plot_file2, bbox_inches='tight')


def flag_statistics(msinfo):
    plots_obs_dir = './plots/plots_flagstats/'
    makedir(plots_obs_dir)
    logger.info('Start flagstatistics')
    logger.info('Running flagdata on {0}, mode="summary", action="calculate", antenna="*&*"'.format(msinfo['msfile']))
    flag_stats = flagdata(vis=msinfo['msfile'], mode='summary', action='calculate', display='none', antenna='*&*')
    save_obj(flag_stats, './plots/plots_flagstats/flagstats')
    logger.info('flagstats file saved to: ./plots/plots_flagstats/flagstats.pkl')
    # For testing (read instead of producing):
#    try:
#        flag_stats = load_obj('./plots/plots_flagstats/flagstats')
#    except:
#       flag_stats = flagdata(vis=msinfo['msfile'], mode='summary', action='calculate', display='none', antenna='*&*') 
#        save_obj(flag_stats, 'flagstats')
    logger.info('Flag statistics ready. Now plotting.')
    plot_flagstatistics(flag_stats, msinfo)
    logger.info('End flagstatistics')



