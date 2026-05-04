#!/usr/local/python
import os
import subprocess
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

from ..utils.weblog_config import get_weblog_function

# Get the appropriate weblog function (original or modern)
start_weblog = get_weblog_function()
from ..utils import eMCP_utils as emutils
from ..utils import eMCP_paths as empaths
from ..functions import eMCP_functions as em

import logging

from casatools import measures
from casatools import ms as my_ms
from casatools import msmetadata

plt.ioff()

ms = my_ms()

logger = logging.getLogger('logger')

weblog_dir = empaths.WEBLOG_DIR
info_dir = empaths.INFO_DIR
calib_dir = empaths.CALIB_DIR
plots_dir = empaths.PLOTS_DIR
logs_dir = empaths.LOGS_DIR
images_dir = empaths.IMAGES_DIR

weblog_link = empaths.WEBLOG_LINK
info_link = empaths.INFO_LINK
calib_link = empaths.CALIB_LINK
plots_link = empaths.PLOTS_LINK
images_link = empaths.IMAGES_LINK

line0 = '-' * 15


def simple_plot_name(plot_file, i):
    try:
        actual_name = glob.glob('{0}{1}_*.png'.format(plot_file, i))[0]
        shutil.move(actual_name, '{0}{1}.png'.format(plot_file, i))
    except:
        pass


# --- shadems helper functions ---------------------------------------------------


def _run_shadems(cmd):
    """Run a shadems command via subprocess, logging output."""
    import shlex
    env = os.environ.copy()
    env.setdefault('MPLBACKEND', 'Agg')
    logger.info('shadems: %s', shlex.join(cmd))
    proc = subprocess.run(
        cmd, env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True, check=False)
    if proc.returncode:
        logger.warning('shadems exited with code %d', proc.returncode)
        lines = [line for line in proc.stdout.splitlines() if line.strip()]
        if lines:
            logger.warning('shadems output tail:\n%s',
                           '\n'.join(lines[-8:]))
        logger.debug(proc.stdout)
    else:
        logger.debug(proc.stdout)
    return proc.returncode


def _get_freq_limits(msfile):
    """Return (freq_min_Hz, freq_max_Hz) across all SPWs."""
    freqs = emutils.read_keyword(
        msfile, 'CHAN_FREQ', subtable='SPECTRAL_WINDOW')
    return float(freqs.min()), float(freqs.max())


def _compute_amp_max(msfile, field, datacolumn):
    """Return the maximum unflagged amplitude for *field* / *datacolumn*.

    Returns None on any failure so callers can fall back to auto-scale.
    """
    from casatools import table as tb_tool
    field_names = emutils.read_keyword(
        msfile, 'NAME', 'FIELD').tolist()
    try:
        field_id = field_names.index(field)
    except ValueError:
        return None
    tb = tb_tool()
    try:
        tb.open(msfile, nomodify=True)
        subtb = tb.query(f'FIELD_ID == {field_id}')
        if subtb.nrows() == 0:
            subtb.close()
            return None
        data = subtb.getcol(datacolumn)
        flags = subtb.getcol('FLAG')
        subtb.close()
        amp = np.abs(data)
        amp[flags] = np.nan
        return float(np.nanmax(amp)) * 1.05
    except Exception as exc:
        logger.warning('Cannot compute amp range for %s/%s: %s',
                       field, datacolumn, exc)
        return None
    finally:
        try:
            tb.close()
        except Exception:
            pass


def _ms_has_column(msfile, column):
    from casatools import table as tb_tool
    tb = tb_tool()
    try:
        tb.open(msfile, nomodify=True)
        return column in tb.colnames()
    except Exception as exc:
        logger.warning('Cannot inspect columns in %s: %s', msfile, exc)
        return False
    finally:
        try:
            tb.close()
        except Exception:
            pass


def _shadems_cmd(msfile, plot_dir, png_name, xaxis, yaxis, field,
                 corr='RR,LL', baseline='noautocorr',
                 xcanvas=1200, ycanvas=200, fontsize=14,
                 title='', norm='eq_hist', data_column=None, extra=None):
    """Build a shadems command list."""
    cmd = [
        'shadems', msfile,
        '--dir', plot_dir,
        '--png', png_name,
        '--title', title,
        '-x', xaxis, '-y', yaxis,
        '--field', field,
        '--corr', corr,
        '--baseline', baseline,
        '-X', str(xcanvas), '-Y', str(ycanvas),
        '--fontsize', str(fontsize),
        '--norm', norm,
        '-j', '1',
    ]
    if data_column:
        cmd.extend(['--col', data_column])
    if extra:
        cmd.extend(extra)
    return cmd


def _safe_name(text):
    """Sanitise *text* for use in a filename."""
    return ''.join(
        c if c.isalnum() or c in '._+-' else '_' for c in text)


# --- visibility 4-plots (shadems) ---------------------------------------------


def batch_4plots(msinfo, fields_list, datacolumn, plots_data_dir, num_proc):
    """Produce per-baseline visibility plots for multiple fields using shadems.

    For each of the four plot types (amp/phase × time/freq) shadems is
    called once with ``--iter-baseline`` and ``--iter-field``, producing
    PNGs for all selected fields and baselines in parallel.
    """
    logger.info('Batch visibility plots (shadems) for datacolumn: %s', datacolumn)
    msfile = msinfo['msfile']
    col = {'data': 'DATA', 'corrected': 'CORRECTED_DATA'}.get(
        datacolumn, datacolumn.upper())
    if not _ms_has_column(msfile, col):
        logger.warning('Cannot make %s plots. Column %s not found in %s.',
                       datacolumn, col, msfile)
        return

    freq_min, freq_max = _get_freq_limits(msfile)
    os.makedirs(plots_data_dir, exist_ok=True)

    plot_specs = [
        ('amp_time', 'TIME', 'amp', []),
        ('phase_time', 'TIME', 'phase', ['--ymin', '-180', '--ymax', '180']),
        ('amp_freq', 'FREQ', 'amp',
         ['--xmin', f'{freq_min:.12g}', '--xmax', f'{freq_max:.12g}',
          '--colour-by', 'CORR']),
        ('phase_freq', 'FREQ', 'phase',
         ['--ymin', '-180', '--ymax', '180',
          '--xmin', f'{freq_min:.12g}', '--xmax', f'{freq_max:.12g}',
          '--colour-by', 'CORR']),
    ]

    field_str = ','.join(fields_list)

    for ptype, xaxis, yaxis, extra in plot_specs:
        # shadems replaces {_field} with _fieldname and {_Baseline} with -Ant1-Ant2
        png = f"{msinfo['msfilename']}{{_field}}_{datacolumn}_{ptype}{{_Baseline}}.png"
        cmd = _shadems_cmd(
            msfile, plots_data_dir, png,
            xaxis=xaxis, yaxis=yaxis, field=field_str,
            data_column=col,
            extra=['--iter-baseline', '--iter-field'] + extra)

        # Override -j with num_proc
        try:
            j_idx = cmd.index('-j')
            cmd[j_idx + 1] = str(num_proc)
        except ValueError:
            cmd.extend(['-j', str(num_proc)])

        _run_shadems(cmd)

    logger.info('shadems 4-plots batch completed')


def make_4plots(eMCP, datacolumn='data'):
    logger.info(line0)
    msinfo = eMCP['msinfo']
    logger.info('Start plot_{}'.format(datacolumn))
    t0 = datetime.datetime.now(datetime.timezone.utc)
    if datacolumn == 'data':
        plots_data_dir = empaths.PLOTS_DATA_DIR
    elif datacolumn == 'corrected':
        plots_data_dir = empaths.PLOTS_CORRECTED_DIR
    else:
        plots_data_dir = empaths.PLOTS_DIR
    emutils.makedir(plots_data_dir)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    logger.info('Producing plots for: {}'.format(','.join(allsources)))
    import multiprocessing
    num_proc = max(1, multiprocessing.cpu_count() - 1)

    valid_fields = [f for f in allsources if f in mssources]
    if valid_fields:
        batch_4plots(msinfo, valid_fields, datacolumn, plots_data_dir, num_proc)

    for field in allsources:
        if field not in mssources:
            logger.warning('Cannot plot {0}. Source not in ms.'.format(field))

    logger.info('Visibility plots finished')
    if datacolumn == 'corrected':
        make_uvplt(eMCP)
    logger.info('End plot_{}'.format(datacolumn))
    msg = ''
    eMCP = em.add_step_time('plot_' + datacolumn, eMCP, msg, t0)
    return eMCP


# --- UV-distance plots (shadems) -----------------------------------------------


def single_uvplt(msinfo, field, plots_data_dir):
    """Corrected amplitude and phase vs uv-distance using shadems."""
    logger.info('uvplt (shadems) for field: %s', field)
    msfile = msinfo['msfile']
    if not _ms_has_column(msfile, 'CORRECTED_DATA'):
        logger.warning('Cannot make corrected uvplt. Column CORRECTED_DATA '
                       'not found in %s.', msfile)
        return
    prefix = f"{msinfo['msfilename']}_uvplt_{_safe_name(field)}"
    os.makedirs(plots_data_dir, exist_ok=True)

    # Corrected amplitude
    _run_shadems(_shadems_cmd(
        msfile, plots_data_dir,
        f'{prefix}_corrected_amp.png',
        xaxis='uv', yaxis='amp', field=field,
        data_column='CORRECTED_DATA',
        ycanvas=600, title=f'Corrected Amp vs UV-dist  {field}'))

    # Corrected phase
    _run_shadems(_shadems_cmd(
        msfile, plots_data_dir,
        f'{prefix}_corrected_phase.png',
        xaxis='uv', yaxis='phase', field=field,
        data_column='CORRECTED_DATA',
        ycanvas=600, title=f'Corrected Phase vs UV-dist  {field}',
        extra=['--ymin', '-180', '--ymax', '180']))


def single_uvplt_model(msinfo, field, plots_data_dir, amp_max=None):
    """Model amplitude and phase vs uv-distance using shadems.

    If *amp_max* is given the amplitude plot uses the same y-range as
    the corresponding corrected plot.
    """
    logger.info('uvplt model (shadems) for field: %s', field)
    msfile = msinfo['msfile']
    if not _ms_has_column(msfile, 'MODEL_DATA'):
        logger.warning('Cannot make model uvplt. Column MODEL_DATA not '
                       'found in %s.', msfile)
        return
    prefix = f"{msinfo['msfilename']}_uvplt_{_safe_name(field)}"
    os.makedirs(plots_data_dir, exist_ok=True)

    amp_extra = []
    if amp_max is not None:
        amp_extra = ['--ymin', '0', '--ymax', str(amp_max)]

    # Model amplitude (same y-range as corrected)
    _run_shadems(_shadems_cmd(
        msfile, plots_data_dir,
        f'{prefix}_model_amp.png',
        xaxis='uv', yaxis='amp', field=field,
        data_column='MODEL_DATA',
        ycanvas=600, title=f'Model Amp vs UV-dist  {field}',
        extra=amp_extra))

    # Model phase
    _run_shadems(_shadems_cmd(
        msfile, plots_data_dir,
        f'{prefix}_model_phase.png',
        xaxis='uv', yaxis='phase', field=field,
        data_column='MODEL_DATA',
        ycanvas=600, title=f'Model Phase vs UV-dist  {field}',
        extra=['--ymin', '-180', '--ymax', '180']))


def make_uvplt(eMCP):
    msinfo = eMCP['msinfo']
    msfile = msinfo['msfile']
    plots_data_dir = empaths.PLOTS_UVPLT_DIR
    emutils.makedir(plots_data_dir)
    if not _ms_has_column(msfile, 'CORRECTED_DATA'):
        logger.warning('Cannot make corrected uvplt plots. Column '
                       'CORRECTED_DATA not found in %s.', msfile)
        return
    has_model_data = _ms_has_column(msfile, 'MODEL_DATA')
    if not has_model_data:
        logger.warning('Cannot make model uvplt plots. Column MODEL_DATA '
                       'not found in %s.', msfile)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    calsources = [s.strip() for s in
                  msinfo['sources']['calsources'].split(',')]
    logger.info('Producing uvplot for: %s', ','.join(allsources))

    for field in allsources:
        if field not in mssources:
            logger.warning('Cannot plot %s. Source not in ms.', field)
            continue
        # --- corrected amp/phase ---
        single_uvplt(msinfo, field, plots_data_dir)

        # --- model amp/phase (calibrators only) ---
        if has_model_data and field in calsources:
            # Compute corrected amplitude range so the model plot matches
            amp_max = _compute_amp_max(msfile, field, 'CORRECTED_DATA')
            single_uvplt_model(msinfo, field, plots_data_dir,
                               amp_max=amp_max)

    logger.info('uvplts finished')

def make_uvcov(msfile, msinfo):
    """Produce V vs U coverage plots using shadems (one PNG per field)."""
    import multiprocessing
    plots_obs_dir = empaths.PLOTS_OBSERVATION_DIR
    emutils.makedir(plots_obs_dir)
    allsources = msinfo['sources']['allsources'].split(',')
    mssources = msinfo['sources']['mssources'].split(',')
    valid_fields = [f for f in allsources if f in mssources]

    if not valid_fields:
        logger.warning('No valid fields for uvcov plotting')
        return

    num_proc = max(1, multiprocessing.cpu_count() - 1)
    field_str = ','.join(valid_fields)

    logger.info('Plotting uvcov (shadems) for: %s', field_str)

    # shadems {_field} placeholder produces -fieldname in the filename
    png = f"{msinfo['msfilename']}_uvcov{{_field}}.png"
    cmd = _shadems_cmd(
        msfile, plots_obs_dir, png,
        xaxis='u', yaxis='v', field=field_str,
        corr='RR', xcanvas=900, ycanvas=900,
        title=f"UV coverage {{_field}}",
        extra=['--iter-field'])
    # Override -j
    try:
        j_idx = cmd.index('-j')
        cmd[j_idx + 1] = str(num_proc)
    except ValueError:
        cmd.extend(['-j', str(num_proc)])
    _run_shadems(cmd)

    for f in allsources:
        if f not in mssources:
            logger.warning(
                'Cannot plot uvcov for %s. Source not in ms.', f)


def _phasecenter_for_epoch(msmd, field_id, epoch):
    """Return the field phase centre, using epoch-aware metadata if available."""
    try:
        return msmd.phasecenter(field_id, epoch)
    except TypeError:
        return msmd.phasecenter(field_id)


def _plot_elevation_track(ax, msmd, me, field_id, field_name):
    times = np.asarray(msmd.timesforfield(field_id), dtype=float)
    if times.size == 0:
        logger.warning('No times found for field %s. Skipping elevation plot.',
                       field_name)
        return False

    elevations = []
    valid_times = []
    for timestamp in times:
        epoch = me.epoch('utc', f'{timestamp / 86400.0}d')
        me.doframe(epoch)
        direction = _phasecenter_for_epoch(msmd, field_id, epoch)
        azel = me.measure(direction, 'azel')
        elevations.append(np.degrees(azel['m1']['value']))
        valid_times.append(timestamp)

    plot_times = Time(np.asarray(valid_times) / 86400.0,
                      format='mjd', scale='utc').datetime
    ax.plot(plot_times, elevations, '.', ms=4, label=field_name)
    return True


def make_elevation(msfile, msinfo):
    plots_obs_dir = empaths.PLOTS_OBSERVATION_DIR
    emutils.makedir(plots_obs_dir)
    plot_file = plots_obs_dir + '{0}_elevation.png'.format(
        msinfo['msfilename'])
    logger.info('Plotting elevation to:')
    logger.info('{}'.format(plot_file))

    msmd = msmetadata()
    me = measures()
    plotted = False
    fig, ax = plt.subplots(figsize=(9, 6))

    try:
        msmd.open(msfile)
        me.doframe(msmd.observatoryposition())

        field_names = list(msmd.namesforfields())
        allsources = msinfo['sources']['allsources'].split(',')
        mssources = msinfo['sources']['mssources'].split(',')
        valid_sources = [source for source in allsources if source in mssources]

        for field_id, field_name in enumerate(field_names):
            if valid_sources and field_name not in valid_sources:
                continue
            plotted |= _plot_elevation_track(
                ax, msmd, me, field_id, field_name)

        for field_name in allsources:
            if field_name not in mssources:
                logger.warning(
                    'Cannot plot elevation for %s. Source not in ms.',
                    field_name)

    finally:
        try:
            msmd.done()
        except Exception:
            pass

    if not plotted:
        plt.close(fig)
        logger.warning('No valid elevation data found. Plot not created.')
        return

    ax.set_xlabel('UTC time')
    ax.set_ylabel('Elevation [deg]')
    ax.set_ylim(0, 90)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize='small')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d %H:%M'))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(plot_file, dpi=100)
    plt.close(fig)


# Flag statistics
def fperc(x):
    return 1.0 * x['flagged'] / x['total']


def sort_list(item, flagged, list_order):
    if len(item) == 0:
        return np.array([]), np.array([])
    order = {a: i for i, a in enumerate(list_order)}
    sorted_pairs = sorted(
        zip(item, flagged),
        key=lambda d: (order.get(d[0], len(order)), str(d[0])))
    item_sorted, flagged_sorted = np.asarray(sorted_pairs, dtype=object).T
    return item_sorted, np.asarray(flagged_sorted, dtype=float)


def read_scan_summary(datain):
    ms.open(datain)
    scan_summary = ms.getscansummary()
    ms.close()
    return scan_summary


def count_flags(flag_stats, label, list_order=None):
    if list_order is None:
        list_order = []
    item = []
    flagged = []
    for s in flag_stats.get(label, {}).keys():
        stats = flag_stats[label][s]
        if stats.get('total', 0) == 0:
            logger.debug('Skipping empty flag statistics item: %s/%s', label, s)
            continue
        flagged.append(fperc(stats))
        try:
            item.append(int(s))
        except:
            item.append(s)
    if len(item) == 0:
        return np.array([]), np.array([])
    if len(list_order) == 0:
        order = np.array(item).argsort()
        item_sorted = np.array(item)[order]
        flagged_sorted = np.array(flagged, dtype=float)[order]
    else:
        item_sorted, flagged_sorted = sort_list(item,
                                                flagged,
                                                list_order=list_order)
    return item_sorted, flagged_sorted


def _mark_empty_axis(ax, message):
    ax.text(0.5, 0.5, message, ha='center', va='center',
            transform=ax.transAxes)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def plot_flagstatistics(flag_stats, msinfo, step):
    """Create flag statistics plots for a processing step.

    Creates two separate plots:
    1. A plot showing only scan flags with detailed view
    2. A plot showing field, spw, and antenna flags

    Both plots have the same dimensions for consistency.
    """
    # Read MS information for field identification
    msfile = msinfo['msfile']
    scan_number = emutils.read_keyword(msfile, 'SCAN_NUMBER')
    field_id = emutils.read_keyword(msfile, 'FIELD_ID')
    scan_fieldID_dict = {}
    for scan in np.unique(scan_number):
        scan_fieldID_dict[str(scan)] = np.unique(
            field_id[np.where(scan_number == scan)[0]])[0]
    vis_fields = emutils.read_keyword(msfile, 'NAME', 'FIELD').tolist()

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
    f_ant = np.array(f_ant, dtype=float)
    f_field = np.array(f_field, dtype=float)

    # Create output directory
    plots_obs_dir = empaths.PLOTS_FLAGSTATS_DIR
    emutils.makedir(plots_obs_dir)

    # Define common figure size for both plots
    figsize = (25, 4)

    # 1. First plot - Only scan flags
    fig_scans = plt.figure(figsize=figsize)
    ax_scan = fig_scans.add_subplot(111)

    if len(i_scan) == 0:
        logger.warning('No scan flag statistics available for step %s', step)
        _mark_empty_axis(ax_scan, 'No scan flag statistics available')
    else:
        # Map scan to field colors
        scan_fieldID = np.array([
            scan_fieldID_dict.get(str(si), -1) for si in i_scan
        ])

        if len(i_field) == 0:
            ax_scan.bar(i_scan - 0.5,
                        f_scan,
                        alpha=1.0,
                        color='0.5',
                        width=1,
                        zorder=10)
        else:
            # Plot bars for each field with different colors
            for i, fi in enumerate(i_field):
                cond = scan_fieldID == i
                if np.any(cond):
                    ax_scan.bar(i_scan[cond] - 0.5,
                                f_scan[cond],
                                alpha=1.0,
                                color=plt.cm.Set1(1.0 * i / len(i_field)),
                                width=1,
                                label='{0} ({1})'.format(fi, i),
                                zorder=10)

    # Scan plot styling
    ax_scan.grid(axis='y', ls='-', color='0.6', zorder=-1000)
    handles, labels = ax_scan.get_legend_handles_labels()
    if handles:
        ax_scan.legend(loc=2, fontsize=7, ncol=4)
    ax_scan.xaxis.set_major_locator(MultipleLocator(10))
    if len(i_scan) > 0:
        ax_scan.set_xlim(np.min(i_scan) - 0.5, np.max(i_scan) + 0.5)
    ax_scan.set_ylim(0, 1)
    ax_scan.set_xlabel('Scan number')
    ax_scan.set_ylabel('Flagged fraction')
    ax_scan.set_title('Scan Flags - {0}'.format(step))

    # Save scan plot
    plot_file_scans = plots_obs_dir + '{0}_flagstats_scans_{1}.png'.format(
        msinfo['msfilename'], step)
    fig_scans.savefig(plot_file_scans, bbox_inches='tight')
    plt.close(fig_scans)

    # 2. Second plot - Field, SPW, and Antenna flags
    fig_other = plt.figure(figsize=figsize)
    plt.subplots_adjust(wspace=0.01)

    # Create three side-by-side subplots for field, spw, and antenna
    ax_field = fig_other.add_subplot(131)
    ax_spw = fig_other.add_subplot(132, sharey=ax_field)
    ax_ant = fig_other.add_subplot(133, sharey=ax_field)

    # Plot Field flags
    if len(i_field) == 0:
        logger.warning('No field flag statistics available for step %s', step)
        _mark_empty_axis(ax_field, 'No field flag statistics available')
    for i, (fi, field_value) in enumerate(zip(i_field, f_field)):
        ax_field.bar(i,
                  field_value,
                  alpha=1.0,
                  color=plt.cm.Set1(1.0 * i / len(i_field)),
                  width=1,
                  label='{0} ({1})'.format(fi, i),
                  align='center',
                  zorder=10)
        ax_field.text(i - 0.1,
                   0.9 * float(field_value),
                   "{0:2.0f}".format(float(field_value) * 100.),
                   color='k',
                   va='center',
                   zorder=12)

    # Plot SPW flags
    if len(i_spw) == 0:
        logger.warning('No SPW flag statistics available for step %s', step)
        _mark_empty_axis(ax_spw, 'No SPW flag statistics available')
    else:
        ax_spw.bar(range(len(i_spw)),
                  f_spw,
                  alpha=1.0,
                  color='0.5',
                  width=1,
                  align='center',
                  zorder=10)

    # Plot Antenna flags
    if len(i_ant) == 0:
        logger.warning('No antenna flag statistics available for step %s',
                       step)
        _mark_empty_axis(ax_ant, 'No antenna flag statistics available')
    else:
        ax_ant.bar(range(len(i_ant)),
                  f_ant,
                  alpha=1.0,
                  color='0.5',
                  width=1,
                  align='center',
                  zorder=10)

    # Add text annotations for values
    for i, v in enumerate(f_spw):
        ax_spw.text(i - 0.1,
                   0.9 * v,
                   "{0:2.0f}".format(v * 100.),
                   color='k',
                   va='center',
                   zorder=12)

    for i, v in enumerate(f_ant):
        v = float(v)
        ax_ant.text(i - 0.1,
                   0.9 * v,
                   "{0:2.0f}".format(v * 100.),
                   color='k',
                   va='center',
                   zorder=12)

    # Set axis properties for field plot
    ax_field.set_xticks(range(len(i_field)))
    ax_field.set_xticklabels([])
    ax_field.set_title('Field')
    ax_field.set_ylabel('Flagged fraction')
    ax_field.set_ylim(0, 1)
    if len(i_field) > 0:
        ax_field.set_xlim(-0.5, len(i_field) - 0.5)
    ax_field.grid(axis='y', zorder=-1000, ls='-', color='0.6')

    # Add field names as rotated annotations
    for i, fi in enumerate(i_field):
        ax_field.annotate('{0} ({1})'.format(fi, i),
                         (i + 0.1, 0.95),
                         va='top',
                         ha='right',
                         rotation=90,
                         zorder=100)

    # Set axis properties for spw plot
    ax_spw.set_xticks(range(len(i_spw)))
    ax_spw.set_xticklabels(range(len(i_spw)))
    ax_spw.set_title('SPW')
    ax_spw.set_xlabel('spw')
    ax_spw.set_yticklabels([])
    if len(i_spw) > 0:
        ax_spw.set_xlim(-0.5, len(i_spw) - 0.5)
    ax_spw.grid(axis='y', zorder=-1000, ls='-', color='0.6')

    # Set axis properties for antenna plot
    ax_ant.set_xticks(range(len(i_ant)))
    ax_ant.set_xticklabels(i_ant, rotation=90)
    ax_ant.set_title('Antenna')
    ax_ant.set_yticklabels([])
    if len(i_ant) > 0:
        ax_ant.set_xlim(-0.5, len(i_ant) - 0.5)
    ax_ant.grid(axis='y', zorder=-1000, ls='-', color='0.6')

    # Add overall title
    fig_other.suptitle('Other Flags - {0}'.format(step), fontsize=14)

    # Save other flags plot
    plot_file_other = plots_obs_dir + '{0}_flagstats_other_{1}.png'.format(
        msinfo['msfilename'], step)
    fig_other.savefig(plot_file_other, bbox_inches='tight')
    plt.close(fig_other)

    # Return both plot filenames for use in the weblog
    return plot_file_scans, plot_file_other


def plot_Lo_drops(phscal_scans, scans, amp_mean, lo_dropout_scans, phscal,
                  eMCP):
    plots_obs_dir = empaths.PLOTS_FLAGSTATS_DIR
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
    if len(lo_dropout_scans) > 0:
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

    plots_obs_dir = empaths.PLOTS_FLAGSTATS_DIR
    plot_file_Lo = plots_obs_dir + '{0}_Lo_dropout_scans{1}.png'.format(
        msinfo['msfilename'], phscal)
    fig.savefig(plot_file_Lo, bbox_inches='tight')


def read_calfluxes(calfluxes, k, eMfactor):
    freq = np.array(calfluxes['freq'])
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

    freq = np.array(calfluxes['freq'])
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
    cond3 = ~data['FLAG'][0, 0, :]
    cond = cond1 * cond2 * cond3
    if len(np.unique(data['SPECTRAL_WINDOW_ID'])) > 1:
        color1 = color2 = data['SPECTRAL_WINDOW_ID'][cond]
    else:
        color1, color2 = '#0067cb', '#c67d50'
    logger.debug(f'Num points in plot: {len(value[0,0][cond])}')
    if np.count_nonzero(cond) > 1:
        ax.scatter(tm[cond].datetime64,
                   value[0,0][cond],
                   marker='.',
                   s=s,
                   c=color1,
                   ec='None',
                   alpha=0.5,
                   cmap=plt.get_cmap('winter_r'))
        ax.scatter(tm[cond].datetime64,
                   value[1,0][cond],
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
    cond3 = ~data['FLAG'][0, 0,:]
    cond = cond1 * cond2 * cond3
    if len(np.unique(data['SPECTRAL_WINDOW_ID'])) > 1:
        color1 = color2 = data['SPECTRAL_WINDOW_ID'][cond]
    else:
        color1, color2 = '#0067cb', '#c67d50'
    logger.debug(antenna_name)
    logger.debug(tm[cond].datetime64)
    logger.debug(cond)
    if np.count_nonzero(cond) > 1:
        ax.scatter(tm[cond][0:10].datetime64, value[0,0][cond][0:10])
        ax.scatter(tm[cond].datetime64,
                   value[0,0][cond],
                   marker='.',
                   s=s,
                   c=color1,
                   ec='None',
                   alpha=1.0,
                   cmap=plt.get_cmap('winter_r'))
        ax.scatter(tm[cond].datetime64,
                   value[1,0][cond],
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
    cond3 = ~data['FLAG'][0, 0, :]
    cond = cond1 * cond2 * cond3
    value[data['FLAG']] = np.nan
    s = 80
    spws = np.unique(emutils.read_keyword(caltable, 'SPECTRAL_WINDOW_ID'))
    all_freqs = emutils.read_keyword(caltable,
                                     'CHAN_FREQ',
                                     subtable='SPECTRAL_WINDOW').T
    plotted = False
    for spw in spws:
        cond4 = data['SPECTRAL_WINDOW_ID'] == spw
        cond = cond1 * cond2 * cond4
        freq = all_freqs[spw] / 1e9
        idxs = np.where(cond)[0]
        if len(idxs) == 0:
            logger.debug('No bandpass rows for antenna %s spw %s in %s',
                         antenna_name, spw, caltable)
            continue
        idx = idxs[0]
        plotted = True

        ax.scatter(freq, value[0,:,idx], marker='.', s=s, c='#0067cb')
        ax.scatter(freq, value[1,:,idx], marker='.', s=s, c='#c67d50')
        ax.errorbar(freq,
                    value[0,:,idx],
                    value_err[0,:,idx],
                    marker='.',
                    ls='',
                    color='#0067cb',
                    ms=1,
                    alpha=0.5)
        ax.errorbar(freq,
                    value[1,:,idx],
                    value_err[1,:,idx],
                    marker='.',
                    ls='',
                    color='#c67d50',
                    ms=1,
                    alpha=0.5)
    if not plotted:
        _mark_empty_axis(ax, 'No bandpass solutions')
    ax.annotate(antenna_name, (0.01, 0.9), xycoords='axes fraction')
    ax.set_xlabel('Freq [GHz]')
    return ax


def plot_caltable(caltable, filename, gaintype='G', calmode=''):
    logger.debug(f'Plotting {caltable}')
    data = emutils.read_caltable_data(caltable)
    antenna_names = em.get_antennas(caltable)
    num_antennas = len(antenna_names)
    points_in_table = len(data.get('TIME', []))
    if num_antennas == 0 or points_in_table == 0:
        logger.warning('No plottable data in caltable %s', caltable)
        fig, ax = plt.subplots(figsize=(10, 4))
        _mark_empty_axis(ax, 'No plottable calibration data')
        fig.savefig(filename, bbox_inches='tight')
        plt.close(fig)
        return
    fig, axes = plt.subplots(nrows=num_antennas,
                             ncols=1,
                             sharex=True,
                             figsize=(10, 14))
    axes = np.atleast_1d(axes)
    fig.subplots_adjust(hspace=0)
    logger.debug(f"Points in table: {len(data['TIME'])}")
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
    #ToDo I need to make this work when running tclean
    # img_tmp.header['BUNIT'] = hdu.header['BUNIT']
    # img_tmp.header['BMAJ'] = hdu.header['BMAJ']
    # img_tmp.header['BMIN'] = hdu.header['BMIN']
    # img_tmp.header['BPA'] = hdu.header['BPA']
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
        logger.info(f"levels: {levels*1000}")
        f.show_contour(levels=levels, alpha=0.4)
    else:
        pass
        #levels = 3. * rms * 1000. * np.sqrt(3)**np.arange(1, 3, 1)
        #logger.debug(f"levels: {levels*1000}")
        #f.show_contour(levels=levels, alpha=0.4)
    output_name = fits_name.replace('.fits', ext + '.png')
    #    logger.info(f'Converting to png fits file {fits_name}')
    plt.savefig(output_name, dpi=200, bbox_inches='tight')
    plt.close()
    f.close()
    emutils.rmfile(fits_name_tmp)
