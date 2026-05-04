import os
import glob
import logging
import numpy as np
from .eMCP_weblog_modern import weblog_header, weblog_foot
from ..utils import eMCP_utils as emutils
from ..utils import eMCP_paths as empaths

logger = logging.getLogger('logger')
weblog_dir = empaths.WEBLOG_DIR

def weblog_flagstats(msinfo):
    """Create flag statistics page with modern layout and sticky jump-to sidebar (calib-style)."""

    wlog = open(weblog_dir + "flagstats.html", "w")
    weblog_header(wlog, 'Flag statistics', msinfo['run'])

    wlog.write('<div class="calib-layout">\n')
    wlog.write('<div class="calib-main">\n')

    # Define the flag stats steps in the correct order
    flagstats_steps = [
        'run_importfits', 'flag_aoflagger', 'flag_apriori', 'flag_manual',
        'restore_flags', 'flag_manual_avg', 'bandpass', 'initial_gaincal',
        'applycal_all', 'flag_target'
    ]

    # Get list of steps that have plots
    available_steps = []
    for step in flagstats_steps:
        scan_glob = glob.glob(
            f'{empaths.PLOTS_FLAGSTATS_DIR}*_flagstats_scans_{step}.png')
        other_glob = glob.glob(
            f'{empaths.PLOTS_FLAGSTATS_DIR}*_flagstats_other_{step}.png')
        if scan_glob or other_glob:
            available_steps.append(step)

    # Process each step in the defined order
    prev_perc_flagged = 0.0
    for step in flagstats_steps:
        flag_stats_file = '{}flagstats_{}.yaml'.format(
            empaths.PLOTS_FLAGSTATS_DIR, step)
        if not os.path.isfile(flag_stats_file):
            continue

        try:
            flag_stats = emutils.load_obj(flag_stats_file)
            total_flags = flag_stats.get('total', 0)
            if total_flags == 0:
                logger.warning('Skipping empty flag statistics file: %s',
                               flag_stats_file)
                continue
            perc_flagged = flag_stats['flagged'] / total_flags * 100.
            diff_flagged = perc_flagged - prev_perc_flagged
            prev_perc_flagged = perc_flagged

            scan_plot_list = glob.glob(
                '{}*_flagstats_scans_{}.png'.format(
                    empaths.PLOTS_FLAGSTATS_DIR, step))
            scan_plot = scan_plot_list[0] if scan_plot_list else None

            other_plot_list = glob.glob(
                '{}*_flagstats_other_{}.png'.format(
                    empaths.PLOTS_FLAGSTATS_DIR, step))
            other_plot = other_plot_list[0] if other_plot_list else None

            if scan_plot is not None or other_plot is not None:
                wlog.write('<div id="{0}" class="subsection">\n'.format(step))
                wlog.write('  <h3 class="collapsible-header">\n')
                wlog.write(f'    {step}')
                wlog.write('    <span class="stats-label">\n')
                wlog.write('      (Total: <span class="total-stat">{0:3.1f}%</span>\n'.format(perc_flagged))
                wlog.write('      Increase: <span class="increase-stat">{0:3.1f}%</span>)\n'.format(diff_flagged))
                wlog.write('    </span>\n')
                wlog.write('  </h3>\n')
                wlog.write('  <div>\n')

                if scan_plot and os.path.isfile(scan_plot):
                    wlog.write('<a href=".{0}" target="_blank">\n'.format(scan_plot))
                    wlog.write('  <img style="max-width:1200px" src=".{0}">\n'.format(scan_plot))
                    wlog.write('</a><br>\n')

                if other_plot and os.path.isfile(other_plot):
                    wlog.write('<a href=".{0}" target="_blank">\n'.format(other_plot))
                    wlog.write('  <img style="max-width:1200px" src=".{0}">\n'.format(other_plot))
                    wlog.write('</a><br>\n')

                wlog.write('<hr>\n')
                wlog.write('  </div>\n')
                wlog.write('</div>\n')
        except Exception as e:
            logger.warning(f'Error processing flag stats for {step}: {e}')
            wlog.write('  </div>\n')
            wlog.write('</div>\n')

    # Flag summary table if available
    flag_summary_file = os.path.join(empaths.FLAGSTATS_DIR,
                                     'flag_summary.txt')
    if os.path.isfile(flag_summary_file):
        wlog.write('<div class="card mt-4 mb-4">\n')
        wlog.write('  <div class="card-header">\n')
        wlog.write('    <h3 class="mb-0">Flag Summary Table</h3>\n')
        wlog.write('  </div>\n')
        wlog.write('  <div class="card-body">\n')

        try:
            with open(flag_summary_file, 'r') as f:
                flag_data = f.readlines()

            wlog.write('    <div class="table-responsive">\n')
            wlog.write('      <table class="table table-sm table-striped">\n')
            header_found = False

            for line in flag_data:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if not header_found:
                    headers = line.split()
                    wlog.write('        <thead>\n')
                    wlog.write('          <tr>\n')
                    for header in headers:
                        wlog.write(f'            <th>{header}</th>\n')
                    wlog.write('          </tr>\n')
                    wlog.write('        </thead>\n')
                    wlog.write('        <tbody>\n')
                    header_found = True
                else:
                    values = line.split()
                    wlog.write('          <tr>\n')
                    for value in values:
                        wlog.write(f'            <td>{value}</td>\n')
                    wlog.write('          </tr>\n')

            wlog.write('        </tbody>\n')
            wlog.write('      </table>\n')
            wlog.write('    </div>\n')

        except Exception as e:
            logger.warning(f'Error parsing flag summary: {e}')
            wlog.write('    <p>Error reading flag summary data.</p>\n')

        wlog.write('  </div>\n')
        wlog.write('</div>\n')

    wlog.write('</div>')  # close calib-main

    # --- Sticky sidebar with Jump Links ---
    wlog.write('<nav class="calib-jump-sidebar">\n')
    wlog.write('<h4>Jump to section</h4>\n')
    wlog.write('<div class="calib-jump-links">\n')
    for step in available_steps:
        wlog.write(f'<a class="calib-jump-link" href="#{step}">{step}</a>\n')
    wlog.write('</div>\n')
    wlog.write('</nav>\n')

    wlog.write('</div>')  # close calib-layout

    weblog_foot(wlog)
    wlog.close()
