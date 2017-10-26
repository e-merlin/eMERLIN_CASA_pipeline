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

def single_4plot(msfile, field, baseline, datacolumn, plot_file):
    logger.info('Baseline: {0:6s}, Output: {1}'.format(baseline, plot_file))
    tb.open(msfile+'/SPECTRAL_WINDOW')
    nchan = str(tb.getcol('NUM_CHAN')[0])
    tb.close()
    gridrows = gridcols = 2
    showgui=False
    avgtime = '300'

    plotms(vis=msfile, xaxis='time', yaxis='amp', title='Amp vs Time (color spw)',
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=0, plotindex=0,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis = 'spw',
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4)

    plotms(vis=msfile, xaxis='time', yaxis='phase', title='Amp vs Time (color spw)',

    gridrows=gridrows, gridcols=gridcols, rowindex=1, colindex=0, plotindex=1,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgchannel = nchan, #avgtime='16',
    xselfscale = True, xsharedaxis = True, coloraxis   = 'spw', clearplots=False,
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4)

    plotms(vis=msfile, xaxis='freq', yaxis='amp', title='Amp vs Frequency (color corr)',
    gridrows=gridrows, gridcols=gridcols, rowindex=0, colindex=1, plotindex=2,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', clearplots=False,
    plotfile = '', expformat = 'png', customsymbol = True, symbolshape = 'circle',
    overwrite=True,  showgui=showgui, symbolsize=4)

    plotms(vis=msfile, xaxis='freq', yaxis='phase', title='Phase vs Frequency (color corr)',
    gridrows=gridrows, gridcols=gridcols, rowindex=1, colindex=1, plotindex=3,
    xdatacolumn=datacolumn, ydatacolumn=datacolumn,correlation = 'RR, LL',
    antenna=baseline, field=field,
    averagedata = True, avgtime=avgtime,
    xselfscale = True, xsharedaxis = True, coloraxis   = 'corr', clearplots=False,
    plotfile = plot_file, expformat = 'png', customsymbol = True, symbolshape = 'circle',
    width=2600, height=1200, symbolsize=4,
    overwrite=True,  showgui=showgui)

def make_4plots(msfile, msinfo, datacolumn='data'):
    plots_data_dir = msinfo['plots_dir']+'plots_data/'
    makedir(plots_data_dir)
    for f in msinfo['sources']['allsources'].split(','):
        logger.info('Generating plot for msfile: {0}, field: {1}, column: {2}'.format(
                                                    (msfile), f, datacolumn))
        for baseline in msinfo['baselines']:
            plot_file = plots_data_dir+'{0}_4plot_{1}_{2}_{3}.png'.format(msinfo['run'], f,
                                                         baseline.replace('&','-'),
                                                         datacolumn)
            single_4plot(msfile, f, baseline, datacolumn, plot_file)




def single_uvcov(msfile, field, plot_file, freqs):
    c = 299792458.
    # Color scale:
    norm = matplotlib.colors.Normalize(vmin=np.min(freqs/1e9),
                                       vmax=np.max(freqs/1e9))
    cmap = matplotlib.cm.get_cmap('Spectral')

    tb.open(msfile)
    uvw = tb.getcol('UVW')
    u = uvw[0,:]; v = uvw[1,:];
    fields = tb.getcol('FIELD_ID')
    tb.close()
    tb.open(msfile+'/FIELD')
    unique_field = tb.getcol('NAME')
    tb.close()
    field_id = np.argwhere(unique_field==field)[0][0]
    sel = (fields == field_id)
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax1 = fig.add_axes([0.8, 0.1, 0.02, 0.8])
    for i, freqi in enumerate(freqs):
        col= cmap(norm(freqi/1e9)) #color
        fc = freqi/c/1e6
        ax.plot(+u[sel]*fc, +v[sel]*fc, marker='.', ms=0.01, ls='',
                 color=col, mec=col, label='{0:4.1f}GHz'.format(freqi/1e9))
        ax.plot(-u[sel]*fc, -v[sel]*fc, marker='.', ms=0.01, ls='',
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

def get_freqs(msfile, allfreqs=False):
    ms.open(msfile)
    axis_info = ms.getdata2(['axis_info'],ifraxis=True)
    ms.close()
    if allfreqs:
        channels = axis_info['axis_info']['freq_axis']['chan_freq']
        freqs=np.sort(axis_info['axis_info']['freq_axis']['chan_freq'].flatten())[::len(channels)/4]
    else:
        freqs = axis_info['axis_info']['freq_axis']['chan_freq'].mean(axis=0)
    return freqs

def make_uvcov(msfile, msinfo):
    plots_uvcov_dir = msinfo['plots_dir']+'plots_uvcov/'
    makedir(plots_uvcov_dir)
    freqs = get_freqs(msfile, allfreqs=True)
    for f in msinfo['sources']['allsources'].split(','):
        plot_file = plots_uvcov_dir+'{0}_uvcov_{1}.png'.format(msinfo['run'],f)
        single_uvcov(msfile, f, plot_file, freqs)



