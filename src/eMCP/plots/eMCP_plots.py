#!/usr/local/python
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from astropy.time import Time
from astropy.io import fits
from astropy.wcs import WCS
from astropy.constants import c as light_speed

import aplpy

from matplotlib.ticker import MultipleLocator, MaxNLocator
from matplotlib.ticker import ScalarFormatter
import datetime
import shutil
import glob

from ..weblog import eMCP_weblog as emwlog
from ..utils import eMCP_utils as emutils
from ..functions import eMCP_functions as em

import logging

from casaplotms import plotms
from casatools import ms as my_ms

plt.ioff()

ms = my_ms()

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


def add_step_time(step, eMCP, msg, t0, doweblog=True):
    t1 = datetime.datetime.utcnow()
    timestamp = t1.strftime('%Y-%m-%d %H:%M:%S')
    delta_t_min = (t1 - t0).total_seconds() / 60.
    eMCP['steps'][step] = [timestamp, delta_t_min, msg]
    emutils.save_obj(eMCP, info_dir + 'eMCP_info.pkl')
    os.system('cp eMCP.log {}eMCP.log.txt'.format(info_dir))
    if doweblog:
        emwlog.start_weblog(eMCP)
    return eMCP


def simple_plot_name(plot_file, i):
    try:
        actual_name = glob.glob('{0}{1}_*.png'.format(plot_file, i))[0]
        shutil.move(actual_name, '{0}{1}.png'.format(plot_file, i))
    except:
        pass


def single_4plot(msinfo, field, datacolumn, plots_data_dir):
    logger.info('Visibility plots for field: {0}, datacolumn: {1}'.format(
        field, datacolumn))
    msfile = msinfo['msfile']
    nchan = str(msinfo['nchan'])

    plot_file = os.path.join(
        plots_data_dir, f"{msinfo['msfilename']}_4plot_{field}_{datacolumn}")
    num_baselines = len(msinfo['baselines'])

    gridrows = num_baselines
    gridcols = 1
    showgui = False
    avgtime = '300'
    baseline = '*&*'
    # Find min and max times per field
    x_min_time, x_max_time = emutils.find_source_timerange(msfile, field)

    w, h = 1200, num_baselines * 200
    plotms(vis=msfile,
           xaxis='time',
           yaxis='amp',
           title='Amp vs Time {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR, LL',
           antenna=baseline,
           field=field,
           iteraxis='baseline',
           averagedata=True,
           avgchannel=nchan,
           avgtime='4',
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotrange=[x_min_time, x_max_time, 0, -1],
           plotfile=plot_file + '0.png',
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           width=w,
           height=h,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()
    simple_plot_name(plot_file, 0)

    plotms(vis=msfile,
           xaxis='time',
           yaxis='phase',
           title='Phase vs Time {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR, LL',
           antenna=baseline,
           field=field,
           iteraxis='baseline',
           averagedata=True,
           avgchannel=nchan,
           avgtime='4',
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotrange=[x_min_time, x_max_time, -180, 180],
           plotfile=plot_file + '1.png',
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           width=w,
           height=h,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()
    simple_plot_name(plot_file, 1)

    plotms(vis=msfile,
           xaxis="freq",
           yaxis="amp",
           title=f"Amp vs freq {field} (color=spw)",
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation="RR,LL",
           antenna=baseline,
           field=field,
           iteraxis="baseline",
           averagedata=True,
           avgchannel=nchan,
           avgtime=avgtime,
           xselfscale=True,
           xsharedaxis=True,
           coloraxis="corr",
           plotfile=plot_file + '1.png',
           expformat="png",
           customsymbol=True,
           symbolshape="circle",
           width=w,
           height=h,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)

    em.find_casa_problems()
    simple_plot_name(plot_file, 2)

    plotms(vis=msfile,
           xaxis='freq',
           yaxis='phase',
           title='Phase vs Frequency {0} (color=corr)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR, LL',
           antenna=baseline,
           field=field,
           iteraxis='baseline',
           averagedata=True,
           avgtime=avgtime,
           avgchannel='4',
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='corr',
           plotrange=[-1, -1, -180, 180],
           plotfile=plot_file + '3.png',
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           width=w,
           height=h,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)

    em.find_casa_problems()
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

    logger.info('Visibility plots finished')
    if datacolumn == 'corrected':
        make_uvplt(eMCP)
    logger.info('End plot_{}'.format(datacolumn))
    msg = ''
    eMCP = add_step_time('plot_' + datacolumn, eMCP, msg, t0)
    return eMCP


def single_uvplt(msinfo, field, plots_data_dir):
    logger.info('uvplt for field: {}'.format(field))
    plot_file_p = plots_data_dir + '{0}_uvplt_p_{1}.png'.format(
        msinfo['msfilename'], field)
    plot_file_a = plots_data_dir + '{0}_uvplt_a_{1}.png'.format(
        msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn = 'corrected'
    avgtime = '16'
    showgui = False
    gridrows = 1
    gridcols = 2
    # Amp

    plotms(vis=msfile,
           xaxis='UVwave',
           yaxis='amp',
           title='Amplitude vs UVWave {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR,LL',
           antenna='*&*',
           field=field,
           averagedata=True,
           avgtime=avgtime,
           avgchannel=str(nchan),
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotfile=plot_file_a,
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           symbolsize=4,
           clearplots=True,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()

    # Phase
    plotms(vis=msfile,
           xaxis='UVwave',
           yaxis='phase',
           title='Phase vs UVWave {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=1,
           plotindex=1,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR,LL',
           antenna='*&*',
           field=field,
           averagedata=True,
           avgtime=avgtime,
           avgchannel=str(nchan),
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotrange=[-1, -1, -180, 180],
           plotfile=plot_file_p,
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           width=1200,
           height=573,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()


def single_uvplt_model(msinfo, field, plots_data_dir):
    logger.info('uvplt (model) for field: {}'.format(field))
    plot_file_p = plots_data_dir + '{0}_uvpltmodel_p_{1}.png'.format(
        msinfo['msfilename'], field)
    plot_file_a = plots_data_dir + '{0}_uvpltmodel_a_{1}.png'.format(
        msinfo['msfilename'], field)
    msfile = msinfo['msfile']
    nchan = msinfo['nchan']
    datacolumn = 'model'
    avgtime = '600'
    showgui = False
    gridrows = 1
    gridcols = 2
    # Amp

    plotms(vis=msfile,
           xaxis='UVwave',
           yaxis='amp',
           title='Model Amplitude vs UVWave {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=0,
           plotindex=0,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR,LL',
           antenna='*&*',
           field=field,
           averagedata=True,
           avgtime=avgtime,
           avgchannel=str(int(nchan / 16)),
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotfile=plot_file_a,
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           symbolsize=4,
           clearplots=True,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()

    plotms(vis=msfile,
           xaxis='UVwave',
           yaxis='phase',
           title='Model Phase vs UVWave {0} (color=spw)'.format(field),
           gridrows=gridrows,
           gridcols=gridcols,
           rowindex=0,
           colindex=1,
           plotindex=1,
           xdatacolumn=datacolumn,
           ydatacolumn=datacolumn,
           correlation='RR,LL',
           antenna='*&*',
           field=field,
           averagedata=True,
           avgtime=avgtime,
           avgchannel=str(nchan),
           xselfscale=True,
           xsharedaxis=True,
           coloraxis='spw',
           plotrange=[-1, -1, -180, 180],
           plotfile=plot_file_p,
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           width=1200,
           height=573,
           symbolsize=4,
           clearplots=False,
           overwrite=True,
           showgui=showgui)
    em.find_casa_problems()


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


# UVplot model, calibrators
    calsources = msinfo['sources']['calsources'].split(',')
    for field in calsources:
        if field in mssources:
            single_uvplt_model(msinfo, field, plots_data_dir)
        else:
            logger.warning('Cannot plot {0}. Source not in ms.'.format(field))
    logger.info('uvplts finished')


def make_uvcov(msfile, msinfo):
    plots_obs_dir = './weblog/plots/plots_observation/'
    emutils.makedir(plots_obs_dir)
    # CASA 5.4 has a bug, it selects the uv limits of the first spw
    # I create this manual limit as a compromise
    max_freq = emutils.read_keyword(msfile,
                                    'CHAN_FREQ',
                                    subtable='SPECTRAL_WINDOW').max()
    #'#    msmd.open(msfile)
    c = light_speed.value
    #'#    max_freq = np.max(np.array([msmd.chanfreqs(spw) for spw in
    #'#                                msmd.datadescids()]))
    max_uvdist = 217000.0 / c * max_freq  # For 217 km baseline
    #freqs = get_freqs(msfile, allfreqs=True)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    logger.info('Plotting uvcov for:')
    for f in allsources:
        if f in mssources:
            plot_file = os.path.join(
                plots_obs_dir,
                '{0}_uvcov_{1}.png'.format(msinfo['msfilename'], f))
            logger.info('{0}'.format(f))
            avgtime = '32'
            nchan = msinfo['nchan']
            plotms(
                vis=msfile,
                xaxis='Uwave',
                yaxis='Vwave',
                field=f,
                title=f,
                correlation='RR',
                spw='',
                coloraxis='spw',
                width=900,
                height=900,
                symbolsize=1,
                plotrange=[-max_uvdist, +max_uvdist, -max_uvdist, +max_uvdist],
                averagedata=True,
                avgtime=avgtime,
                avgchannel=str(int(nchan / 8)),
                plotfile=plot_file,
                expformat='png',
                customsymbol=True,
                symbolshape='circle',
                overwrite=True,
                showlegend=False,
                showgui=False)
            em.find_casa_problems()

        else:
            logger.info(
                'Cannot plot uvcov for {0}. Source not in ms.'.format(f))


def make_elevation(msfile, msinfo):
    plots_obs_dir = './weblog/plots/plots_observation/'
    emutils.makedir(plots_obs_dir)
    plot_file = plots_obs_dir + '{0}_elevation.png'.format(
        msinfo['msfilename'])
    logger.info('Plotting elevation to:')
    logger.info('{}'.format(plot_file))
    avgtime = '16'
    showgui = False
    plotms(vis=msfile,
           xaxis='time',
           yaxis='elevation',
           correlation='RR',
           spw='',
           coloraxis='field',
           width=900,
           symbolsize=5,
           plotrange=[-1, -1, 0, 90],
           averagedata=True,
           avgtime=avgtime,
           plotfile=plot_file,
           expformat='png',
           customsymbol=True,
           symbolshape='circle',
           overwrite=True,
           showlegend=True,
           showgui=showgui)
    em.find_casa_problems()


# Flag statistics
def fperc(x):
    return 1.0 * x['flagged'] / x['total']


def sort_list(item, flagged, list_order):
    order = {a: i for i, a in enumerate(list_order)}
    item_sorted, flagged_sorted = np.asarray(
        sorted(zip(item, flagged), key=lambda d: order[d[0]])).T
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
        item_sorted, flagged_sorted = sort_list(item,
                                                flagged,
                                                list_order=list_order)
    return item_sorted, flagged_sorted


def plot_flagstatistics(flag_stats, msinfo, step):
    # Different colors for each field

    msfile = msinfo['msfile']
    scan_number = emutils.read_keyword(msfile, 'SCAN_NUMBER')
    field_id = emutils.read_keyword(msfile, 'FIELD_ID')
    scan_fieldID_dict = {}
    for scan in np.unique(scan_number):
        scan_fieldID_dict[str(scan)] = np.unique(
            field_id[np.where(scan_number == scan)[0]])[0]
    vis_fields = emutils.read_keyword(msfile, 'NAME', 'FIELD')

    # Compute % statistics
    i_scan, f_scan = count_flags(flag_stats, 'scan')
    i_field, f_field = count_flags(flag_stats, 'field', list_order=vis_fields)
    i_corr, f_corr = count_flags(flag_stats,
                                 'correlation',
                                 list_order=['RR', 'LL', 'RL', 'LR'])
    i_spw, f_spw = count_flags(flag_stats, 'spw')
    i_ant, f_ant = count_flags(flag_stats,
                               'antenna',
                               list_order=msinfo['antennas'])

    fig = plt.figure(figsize=(25, 4))
    plt.subplots_adjust(wspace=0.01)
    ax1 = fig.add_subplot(1, 5, (1, 2))
    ax2 = fig.add_subplot(153)
    #    ax3 = fig.add_subplot(153)
    ax4 = fig.add_subplot(154, sharey=ax2)
    ax5 = fig.add_subplot(155, sharey=ax2)

    scan_fieldID = np.array([scan_fieldID_dict[str(si)] for si in i_scan])
    for i, fi in enumerate(i_field):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond] - 0.5,
                f_scan[cond],
                alpha=1.0,
                color=plt.cm.Set1(1.0 * i / len(i_field)),
                width=1,
                label='{0} ({1})'.format(fi, i),
                zorder=10)
        field_value = f_field[np.argwhere(i_field == fi)[0][0]]
        ax2.bar(i,
                field_value,
                alpha=1.0,
                color=plt.cm.Set1(1.0 * i / len(i_field)),
                width=1,
                label='{0} ({1})'.format(fi, i),
                align='center',
                zorder=10)
        ax2.text(i - 0.1,
                 0.9 * field_value,
                 "{0:2.0f}".format(field_value * 100.),
                 color='k',
                 va='center',
                 zorder=12)


#    ax3.bar(range(len(i_corr)), f_corr, alpha=0.5, color='k', width=1, align='center')
    ax4.bar(range(len(i_spw)),
            f_spw,
            alpha=1.0,
            color='0.5',
            width=1,
            align='center',
            zorder=10)
    ax5.bar(range(len(i_ant)),
            f_ant,
            alpha=1.0,
            color='0.5',
            width=1,
            align='center',
            zorder=10)

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

    [
        ax2.annotate('{0} ({1})'.format(si, i), (i + 0.1, 0.95),
                     va='top',
                     ha='right',
                     rotation=90,
                     zorder=100) for i, si in enumerate(i_field)
    ]
    #[ax3.annotate(si, (i+0.5, f_corr[i])) for i, si in enumerate(i_corr)]
    #[ax5.annotate(si, (i+0.3, f_ant[i])) for i, si in enumerate(i_ant)]

    for i, v in enumerate(f_spw):
        ax4.text(i - 0.1,
                 0.9 * v,
                 "{0:2.0f}".format(v * 100.),
                 color='k',
                 va='center',
                 zorder=12)
    for i, v in enumerate(f_ant):
        ax5.text(i - 0.1,
                 0.9 * v,
                 "{0:2.0f}".format(v * 100.),
                 color='k',
                 va='center',
                 zorder=12)

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

    ax1.grid(axis='y', zorder=-1000, ls='-', color='0.6')
    ax2.grid(axis='y', zorder=-1000, ls='-', color='0.6')
    ax4.grid(axis='y', zorder=-1000, ls='-', color='0.6')
    ax5.grid(axis='y', zorder=-1000, ls='-', color='0.6')

    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 1)
    #    ax3.set_ylim(0,1)
    ax4.set_ylim(0, 1)
    ax5.set_ylim(0, 1)

    #'#    ax1.set_xlim(np.min(i_scan)-0.5, np.max(i_scan)+0.5)
    ax1.set_xlim(np.min(i_scan) - 1.0, np.max(i_scan) + 0.)
    ax2.set_xlim(-0.5, len(i_field) - 0.5)
    #    ax3.set_xlim(-0.5, len(i_corr)-0.5)
    ax4.set_xlim(-0.5, len(i_spw) - 0.5)
    ax5.set_xlim(-0.5, len(i_ant) - 0.5)

    #ax1.legend(loc=0)
    #ax2.legend(loc=0)

    plots_obs_dir = './weblog/plots/plots_flagstats/'
    plot_file1 = plots_obs_dir + '{0}_flagstats_{1}.png'.format(
        msinfo['msfilename'], step)
    #logger.info('Plot flagstats: {0}'.format(plot_file1))
    fig.savefig(plot_file1, bbox_inches='tight')

    # Plot only scans:
    fig = plt.figure(figsize=(50, 8))
    ax1 = fig.add_subplot(111)

    for i, fi in enumerate(np.unique(i_field)):
        cond = scan_fieldID == i
        ax1.bar(i_scan[cond] - 0.5,
                f_scan[cond],
                alpha=1.0,
                color=plt.cm.Set1(1.0 * i / len(i_field)),
                width=1,
                label='{0} ({1})'.format(fi, i),
                zorder=10)

    try:
        ax1.legend(loc=0)
    except:
        pass
    ax1.grid(axis='y', zorder=-1000, ls='-', color='0.6')

    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.set_xlim(np.min(i_scan) - 0.5, np.max(i_scan) + 0.5)
    ax1.set_ylim(0, 1)
    ax1.set_xlabel('Scan number')
    ax1.set_ylabel('Flagged fraction')

    plot_file2 = plots_obs_dir + '{0}_flagstats_scans_{1}.png'.format(
        msinfo['msfilename'], step)
    #logger.info('Plot flagstats scans: {0}'.format(plot_file2))
    fig.savefig(plot_file2, bbox_inches='tight')


def plot_Lo_drops(phscal_scans, scans, amp_mean, lo_dropout_scans, phscal,
                  eMCP):
    plots_obs_dir = './weblog/plots/plots_flagstats/'
    emutils.makedir(plots_obs_dir)
    msinfo = eMCP['msinfo']
    drops = np.array([scan in lo_dropout_scans for scan in phscal_scans])
    fig = plt.figure(figsize=(30, 8))
    ax1 = fig.add_subplot(111)

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
    if lo_dropout_scans != []:
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


def read_calfluxes(calfluxes, k, eMfactor):
    freq = calfluxes['freq']
    spws = calfluxes['spwID']
    fieldName = calfluxes[k]['fieldName']
    spindex = calfluxes[k]['spidx'][1]
    espindex = calfluxes[k]['spidxerr'][1]
    S0 = calfluxes[k]['fitFluxd'] * eMfactor
    eS0 = calfluxes[k]['fitFluxdErr'] * eMfactor
    freq0 = calfluxes[k]['fitRefFreq']
    flux = np.ones(len(spws)) * np.nan
    eflux = np.ones(len(spws)) * np.nan
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
    return freq, spws, fieldName, spindex, espindex, S0, eS0, freq0, flux, eflux


def fluxscale_models(calfluxes, eMfactor, msinfo):
    factor_unit = 1e-9
    units = 'Jy'
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)

    freq_min, freq_max = 1e20, 0.
    for k in calfluxes.keys():
        if type(calfluxes[k]) is dict:
            freq, spws, fieldName, spindex, espindex, S0, eS0, freq0, flux, eflux = read_calfluxes(
                calfluxes, k, eMfactor)
            freq_min = np.nanmin([freq_min, np.nanmin(freq)])
            freq_max = np.nanmax([freq_max, np.nanmax(freq)])
            freqspan = np.nanmax(freq) - np.nanmin(freq)
            freqlim = np.array([
                np.nanmin(freq) - 0.025 * freqspan,
                np.nanmax(freq) + 0.1 * freqspan
            ])
            fluxfit = S0 * (freqlim / freq0)**spindex
            ff_min = (S0 - eS0) * (freqlim / freq0)**(spindex - espindex)
            ff_max = (S0 + eS0) * (freqlim / freq0)**(spindex + espindex)
            ff_min2 = (S0 - eS0) * (freqlim / freq0)**(spindex + espindex)
            ff_max2 = (S0 + eS0) * (freqlim / freq0)**(spindex - espindex)
            label = '{0:>9s}: Flux density = {1:6.3f} +/-{2:6.3f}, '\
                    'spidx ={3:5.2f}+/-{4:5.2f}'.format(fieldName,
                        calfluxes[k]['fitFluxd']*eMfactor, calfluxes[k]['fitFluxdErr']*eMfactor,
                        calfluxes[k]['spidx'][1], calfluxes[k]['spidxerr'][1])
            p, = ax1.plot(freqlim * factor_unit,
                          fluxfit,
                          '-',
                          linewidth=2,
                          zorder=-5,
                          label=label)
            color1 = str(p.get_color())
            #ax1.errorbar(freq*factor_unit, flux, eflux, fmt = 'o', color =color1, mec = color1, zorder = 10)
            ax1.plot(freq * factor_unit,
                     flux,
                     marker='o',
                     ls='',
                     color='k',
                     mec='k',
                     zorder=10)
            ax1.fill_between(freqlim * factor_unit,
                             ff_min,
                             ff_max,
                             facecolor=color1,
                             color=color1,
                             alpha=1.0,
                             linewidth=0,
                             zorder=-32)
            ax1.fill_between(freqlim * factor_unit,
                             ff_max2,
                             ff_min2,
                             facecolor=color1,
                             color=color1,
                             alpha=1.0,
                             linewidth=0,
                             zorder=-32)

    freq = calfluxes['freq']
    freqspan = freq_max - freq_min
    freqlim = np.array(
        [freq_min - 0.025 * freqspan, freq_max + 0.1 * freqspan])
    ax1.set_xlabel("Frequency [GHz]")
    ax1.set_ylabel("Flux density [Jy]")
    ax1.set_xlim(freqlim[0] * factor_unit, freqlim[1] * factor_unit)
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
    plot_file = plots_obs_dir + '{0}_fluxscale.png'.format(
        msinfo['msfilename'])
    fig.savefig(plot_file, bbox_inches='tight')


### Plot caltables with matplotlib


def plot_gaintable(data, antenna, ax, calmode='ap', field_id=None, s=60):
    t = data['TIME']
    tm = Time(t / 60 / 60 / 24., format='mjd')
    antenna_id, antenna_name = antenna
    if calmode == 'p':
        value = np.angle(data['CPARAM']) * 180 / np.pi
        ax.set_ylim(-180, 180)
        ax.set_ylabel('Phase [deg]')
    elif calmode == 'ap':
        value = np.abs(data['CPARAM'])
        ax.set_ylabel('Amplitude')


#        ax.set_ylim(bottom=0)
    else:
        value = 0
    if field_id == None:
        cond1 = True
    else:
        cond1 = data['FIELD_ID'] == field_id
    cond2 = data['ANTENNA1'] == antenna_id
    cond3 = ~data['FLAG'][:, 0, 0]
    cond = cond1 * cond2 * cond3
    if len(np.unique(data['SPECTRAL_WINDOW_ID'])) > 1:
        color1 = color2 = data['SPECTRAL_WINDOW_ID'][cond]
    else:
        color1, color2 = '#0067cb', '#c67d50'
    logger.debug(f'Num points in plot: {len(value[cond][:,0,0])}')
    if np.count_nonzero(cond) > 1:
        ax.scatter(tm[cond].datetime64,
                   value[cond][:, 0, 0],
                   marker='.',
                   s=s,
                   c=color1,
                   ec='None',
                   alpha=0.5,
                   cmap=plt.get_cmap('winter_r'))
        ax.scatter(tm[cond].datetime64,
                   value[cond][:, 0, 1],
                   marker='.',
                   s=s,
                   c=color2,
                   ec='None',
                   alpha=0.5,
                   cmap=plt.get_cmap('copper'))
    ax.annotate(antenna_name, (0.01, 0.9), xycoords='axes fraction')
    return ax


def plot_delaytable(data, antenna, ax, calmode='p', field_id=None, s=120):
    t = data['TIME']
    tm = Time(t / 60 / 60 / 24., format='mjd')
    antenna_id, antenna_name = antenna
    value = data['FPARAM']
    if field_id == None:
        cond1 = True
    else:
        cond1 = data['FIELD_ID'] == field_id
    cond2 = data['ANTENNA1'] == antenna_id
    cond3 = ~data['FLAG'][:, 0, 0]
    cond = cond1 * cond2 * cond3
    if len(np.unique(data['SPECTRAL_WINDOW_ID'])) > 1:
        color1 = color2 = data['SPECTRAL_WINDOW_ID'][cond]
    else:
        color1, color2 = '#0067cb', '#c67d50'
    logger.debug(antenna_name)
    logger.debug(tm[cond].datetime64)
    logger.debug(cond)
    if np.count_nonzero(cond) > 1:
        ax.scatter(tm[cond][0:10].datetime64, value[cond][:, 0, 0][0:10])
        ax.scatter(tm[cond].datetime64,
                   value[cond][:, 0, 0],
                   marker='.',
                   s=s,
                   c=color1,
                   ec='None',
                   alpha=1.0,
                   cmap=plt.get_cmap('winter_r'))
        ax.scatter(tm[cond].datetime64,
                   value[cond][:, 0, 1],
                   marker='.',
                   s=s,
                   c=color2,
                   ec='None',
                   alpha=1.0,
                   cmap=plt.get_cmap('copper'))
    ax.annotate(antenna_name, (0.01, 0.9), xycoords='axes fraction')
    ax.set_ylabel('Delay [ns]')
    return ax


def plot_bptable(data, caltable, antenna, ax, calmode='p', field_id=None):
    antenna_id, antenna_name = antenna
    if calmode == 'p':
        value = np.angle(data['CPARAM']) * 180 / np.pi
        value_err = data['PARAMERR']
        ax.set_ylim(-180, 180)
        ax.set_ylabel('Phase [deg]')
    elif calmode == 'ap':
        value = np.abs(data['CPARAM'])
        value_err = np.abs(data['PARAMERR'])
        ax.set_ylabel('Amplitude')


#        ax.set_ylim(bottom=0)
    if field_id == None:
        cond1 = True
    else:
        cond1 = data['FIELD_ID'] == field_id
    cond2 = data['ANTENNA1'] == antenna_id
    cond3 = ~data['FLAG'][:, 0, 0]
    cond = cond1 * cond2 * cond3
    value[data['FLAG']] = np.nan
    s = 80
    spws = np.unique(emutils.read_keyword(caltable, 'SPECTRAL_WINDOW_ID'))
    all_freqs = emutils.read_keyword(caltable,
                                     'CHAN_FREQ',
                                     subtable='SPECTRAL_WINDOW')
    for spw in spws:
        cond4 = data['SPECTRAL_WINDOW_ID'] == spw
        cond = cond1 * cond2 * cond4
        freq = all_freqs[spw] / 1e9
        ax.scatter(freq, value[cond][0, :, 0], marker='.', s=s, c='#0067cb')
        ax.scatter(freq, value[cond][0, :, 1], marker='.', s=s, c='#c67d50')
        ax.errorbar(freq,
                    value[cond][0, :, 0],
                    value_err[cond][0, :, 0],
                    marker='.',
                    ls='',
                    color='#0067cb',
                    ms=1,
                    alpha=0.5)
        ax.errorbar(freq,
                    value[cond][0, :, 1],
                    value_err[cond][0, :, 1],
                    marker='.',
                    ls='',
                    color='#c67d50',
                    ms=1,
                    alpha=0.5)
    ax.annotate(antenna_name, (0.01, 0.9), xycoords='axes fraction')
    ax.set_xlabel('Freq [GHz]')
    return ax


def plot_caltable(caltable, filename, gaintype='G', calmode=''):
    logger.debug(f'Plotting {caltable}')
    data = emutils.read_caltable_data(caltable)
    antenna_names = em.get_antennas(caltable)
    num_antennas = len(antenna_names)
    fig, axes = plt.subplots(nrows=num_antennas,
                             ncols=1,
                             sharex=True,
                             figsize=(10, 14))
    fig.subplots_adjust(hspace=0)
    logger.debug(f"Points in table: {len(data['TIME'])}")
    points_in_table = len(data['TIME'])
    s = 120 + 10 * (30000 / points_in_table)**0.3
    s = np.min([np.max([s, 50]), 200])
    for i, ax in enumerate(axes):
        if gaintype == 'G':
            ax = plot_gaintable(data,
                                antenna=[i, antenna_names[i]],
                                ax=ax,
                                calmode=calmode,
                                field_id=None,
                                s=s)
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter('%Y/%m/%d %H:%M'))
        elif gaintype == 'K':
            plot_delaytable(data,
                            antenna=[i, antenna_names[i]],
                            ax=ax,
                            calmode=calmode,
                            field_id=None,
                            s=s)
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter('%Y/%m/%d %H:%M'))
        elif gaintype == 'B':
            plot_bptable(data,
                         caltable=caltable,
                         antenna=[i, antenna_names[i]],
                         ax=ax,
                         calmode=calmode,
                         field_id=None)
    if gaintype != 'B':
        fig.autofmt_xdate()
    print(f"Saving {filename}")
    fig.savefig(filename, bbox_inches='tight')


import matplotlib

cdict = {
    'red': ((0.00, 1.00, 1.00), (0.20, 0.00, 0.00), (0.40, 0.00,
                                                     0.00), (0.50, 0.40, 0.40),
            (0.60, 0.80, 0.80), (0.80, 1.00, 1.00), (1.00, 0.95, 0.95)),
    'green': ((0.00, 1.00, 1.00), (0.20, 0.50, 0.50), (0.40, 0.85, 0.85),
              (0.50, 0.80, 0.80), (0.60, 0.95, 0.95), (0.80, 0.65, 0.65),
              (1.00, 0.00, 0.00)),
    'blue': ((0.00, 1.00, 1.00), (0.20, 0.95, 0.95), (0.40, 0.20, 0.20),
             (0.50, 0.40, 0.40), (0.60, 0.15, 0.15), (0.80, 0.00, 0.00),
             (1.00, 0.00, 0.00))
}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

import cmasher as cmr


def fits2png(fits_name,
             rms,
             scaling,
             plot_title=None,
             cmap_name='viridis',
             colorbar=True,
             contour=True,
             zoom=False):
    """Make a PNG plot out of a FITS file

    Args:
        fits_name (str): path of fits file
        plot_title (str): plot title, default is name of the fits file
        cmap_name (str): name of colormap, default is viridis
        colorbar (bool): include colorbar, default is True
        contour (bool): include contour, default is True
    """
    # This is a trick because aplpy cannot deal with extra dimensions in fits produced by wsclean
    fits_name_tmp = fits_name + '_tmp'
    hdu = fits.open(fits_name)[0]
    wcs_celestial = WCS(hdu.header).celestial
    img_tmp = fits.PrimaryHDU(hdu.data[0, 0] * 1000.,
                              header=wcs_celestial.to_header())
    img_tmp.header['BUNIT'] = hdu.header['BUNIT']
    img_tmp.header['BMAJ'] = hdu.header['BMAJ']
    img_tmp.header['BMIN'] = hdu.header['BMIN']
    img_tmp.header['BPA'] = hdu.header['BPA']
    img_tmp.writeto(fits_name_tmp, overwrite=True)

    f = aplpy.FITSFigure(fits_name_tmp, figsize=(10, 8))
    if zoom:
        x_center = hdu.header['CRVAL1']
        y_center = hdu.header['CRVAL2']
        new_size_asec = hdu.data.shape[-1] * hdu.header[
            'CDELT2'] * 0.25  # 25% of total image
        f.recenter(x_center,
                   y_center,
                   width=new_size_asec,
                   height=new_size_asec)
        ext = '_zoom'
    else:
        ext = ''
    if plot_title == None:
        plot_title = fits_name.replace('.fits', '')
    plt.title(plot_title)
    vmin = hdu.data.min() * 1000
    vmax = hdu.data.max() * 1000.
    #    vmax = np.min([20*rms, vmax])
    #    rms = np.sqrt(np.mean(np.square(hdu.data)))
    #    logger.debug(f'rms: {rms}')
    #    f.show_colorscale(cmap=cmap_name, vmin=-1*rms, vmax=vmax, stretch='log', vmid=-5.1*rms)
    #    f.show_colorscale(cmap=cmap_name, stretch='log',vmin=vmin, vmid=vmin - (vmax - vmin) / 30. )
    #    f.show_colorscale(cmap=cmap_name, vmax=vmax)
    f.show_colorscale(cmap=plt.get_cmap('cmr.rainforest'),
                      vmin=np.max([vmin, -2 * rms * 1000]),
                      vmax=np.min([vmax, 20 * rms * 1000]))
    f.show_colorscale(cmap=cmr.get_sub_cmap('cmr.rainforest', 0.3, 0.85),
                      vmin=np.max([vmin, -2 * rms * 1000]),
                      vmax=np.min([vmax, 20 * rms * 1000]))
    #    f.show_colorscale(cmap=my_cmap, vmin=vmin, vmax=vmax, stretch='log', vmid=1.1*vmin)
    f.ticks.set_color('k')
    if colorbar:
        f.add_colorbar()
        f.colorbar.set_axis_label_text('mJy')
    if 'BMAJ' in fits.open(fits_name_tmp)[0].header:
        f.add_beam()
        f.beam.set_facecolor('0.5')
        f.beam.set_edgecolor('k')
        f.beam.set_linewidth(1.5)
        f.beam.set_alpha(0.5)
    if contour:
        levels = 3. * rms * 1000. * np.sqrt(3)**np.arange(1, 25, 1)
        logger.debug(f"levels: {levels*1000}")
        f.show_contour(levels=levels, alpha=0.4)
    else:
        levels = 3. * rms * 1000. * np.sqrt(3)**np.arange(1, 3, 1)
        logger.debug(f"levels: {levels*1000}")
        f.show_contour(levels=levels, alpha=0.4)
    output_name = fits_name.replace('.fits', ext + '.png')
    #    logger.info(f'Converting to png fits file {fits_name}')
    plt.savefig(output_name, dpi=200, bbox_inches='tight')
    emutils.rmfile(fits_name_tmp)
