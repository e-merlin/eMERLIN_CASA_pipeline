#!/usr/local/python
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.ioff()

# CASA imports
from taskinit import *
from tasks import *

import logging
logger = logging.getLogger('logger')


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


def single_4plot(msfile, field, baseline, datacolumn, plot_file):
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


def make_4plots(msfile, msinfo, datacolumn='data'):
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
            single_4plot(msfile, f, baseline, datacolumn, plot_file)

def single_uvplt(msinfo, field, plot_file):
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


def make_uvplt(msinfo):
    plots_data_dir = './plots/plots_uvplt/'
    makedir(plots_data_dir)
    for f in msinfo['sources']['allsources'].split(','):
        plot_file =  plots_data_dir+'{0}_uvplt_{1}.png'.format(msinfo['msfilename'], f)
        logger.info('Plotting uvplt for {0}: {1}'.format(f, plot_file))
        single_uvplt(msinfo, f, plot_file)


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
    ax.set_xlabel('U [Mlambda]')
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
    overwrite=True,  showgui=showgui)

