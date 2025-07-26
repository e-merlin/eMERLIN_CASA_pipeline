import os
import glob
import logging
import numpy as np
from .eMCP_weblog_modern import weblog_header, weblog_foot
from ..utils import eMCP_utils as emutils

logger = logging.getLogger('logger')
weblog_dir = './weblog/'

def weblog_flagstats(msinfo):
    """Create flag statistics page with modern layout and sticky jump-to sidebar (calib-style)."""

    wlog = open(weblog_dir + "flagstats.html", "w")
    weblog_header(wlog, 'Flag statistics', msinfo['run'])

    # Use calib-layout for consistency
    wlog.write('''
    <style>
    .calib-layout {
      display: flex;
      flex-direction: row;
      max-width: 2400px;
      margin: 0 auto;
      width: 100%;
    }
    .calib-main {
      flex: 1 1 0;
      padding: 28px 26px 28px 0;
      min-width: 0;
      max-width: 1200px;
    }
    .calib-jump-sidebar {
      width: 220px;
      position: sticky;
      top: 35px;
      height: fit-content;
      align-self: flex-start;
      background: #f8f9fa;
      border-left: 1.5px solid #e3e3e3;
      border-radius: 8px 0 0 8px;
      padding: 18px 16px 16px 16px;
      margin-left: 18px;
      z-index: 10;
    }
    .calib-jump-sidebar h4 {
      font-size: 1.06em;
      margin: 0 0 14px 0;
      color: #444;
      font-weight: 600;
      text-align: left;
    }
    .calib-jump-links {
      display: flex;
      flex-direction: column;
      gap: 0.48em;
    }
    .calib-jump-link {
      display: block;
      padding: 7px 12px;
      background: #f2f2f2;
      color: #34618c;
      border-left: 4px solid #e5e5e5;
      border-radius: 4px;
      text-decoration: none;
      font-size: 1em;
      transition: background .13s, color .13s, border .13s;
      margin-left: 0;
    }
    .calib-jump-link:hover, .calib-jump-link.active {
      background: #ddeefd;
      color: #1279bc;
      border-left: 4px solid #3498db;
    }
    .stats-label {
      display: inline-block;
      margin-left: 5px;
      font-weight: normal;
    }
    .total-stat {
      background-color: #0d6efd;
      color: white;
      padding: 2px 6px;
      border-radius: 4px;
      font-size: 0.9em;
    }
    .increase-stat {
      background-color: #dc3545;
      color: white;
      padding: 2px 6px;
      border-radius: 4px;
      font-size: 0.9em;
    }
    .plot-container {
      margin-bottom: 30px;
    }
    </style>
    ''')

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
        scan_glob = glob.glob(f'./weblog/plots/plots_flagstats/*_flagstats_scans_{step}.png')
        other_glob = glob.glob(f'./weblog/plots/plots_flagstats/*_flagstats_other_{step}.png')
        if scan_glob or other_glob:
            available_steps.append(step)

    # Process each step in the defined order
    prev_perc_flagged = 0.0
    for step in flagstats_steps:
        flag_stats_file = './weblog/plots/plots_flagstats/flagstats_{}.yaml'.format(step)
        if not os.path.isfile(flag_stats_file):
            continue

        try:
            flag_stats = emutils.load_obj(flag_stats_file)
            perc_flagged = flag_stats['flagged'] / flag_stats['total'] * 100.
            diff_flagged = perc_flagged - prev_perc_flagged
            prev_perc_flagged = perc_flagged

            scan_plot_list = glob.glob('./weblog/plots/plots_flagstats/*_flagstats_scans_{}.png'.format(step))
            scan_plot = scan_plot_list[0] if scan_plot_list else None

            other_plot_list = glob.glob('./weblog/plots/plots_flagstats/*_flagstats_other_{}.png'.format(step))
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
    if os.path.isfile('./weblog/flagstats/flag_summary.txt'):
        wlog.write('<div class="card mt-4 mb-4">\n')
        wlog.write('  <div class="card-header">\n')
        wlog.write('    <h3 class="mb-0">Flag Summary Table</h3>\n')
        wlog.write('  </div>\n')
        wlog.write('  <div class="card-body">\n')

        try:
            with open('./weblog/flagstats/flag_summary.txt', 'r') as f:
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

