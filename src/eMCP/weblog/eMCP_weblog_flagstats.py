import os
import glob
import logging
import numpy as np
from .eMCP_weblog_modern import weblog_header, weblog_foot
from ..utils import eMCP_utils as emutils

logger = logging.getLogger('logger')
weblog_dir = './weblog/'

def weblog_flagstats(msinfo):
    """Create flag statistics page following the original implementation with modern styling"""
    # Create the flag statistics page
    wlog = open(weblog_dir + "flagstats.html", "w")
    weblog_header(wlog, 'Flag statistics', msinfo['run'])
    
    # Add minimal styling for the page
    wlog.write('<style>\n')
    wlog.write('.sticky-nav {\n')
    wlog.write('  position: sticky;\n')
    wlog.write('  top: 0;\n')
    wlog.write('  background-color: #f8f9fa;\n')
    wlog.write('  border-bottom: 1px solid #ddd;\n')
    wlog.write('  padding: 10px;\n')
    wlog.write('  z-index: 1000;\n')
    wlog.write('  text-align: center;\n')
    wlog.write('  box-shadow: 0 2px 4px rgba(0,0,0,0.1);\n')
    wlog.write('}\n')
    wlog.write('.nav-link {\n')
    wlog.write('  display: inline-block;\n')
    wlog.write('  margin: 5px;\n')
    wlog.write('  padding: 5px 10px;\n')
    wlog.write('  background-color: #f1f1f1;\n')
    wlog.write('  border-radius: 4px;\n')
    wlog.write('  text-decoration: none;\n')
    wlog.write('  color: #333;\n')
    wlog.write('  font-size: 14px;\n')
    wlog.write('}\n')
    wlog.write('.nav-link:hover {\n')
    wlog.write('  background-color: #e0e0e0;\n')
    wlog.write('}\n')
    wlog.write('.stats-label {\n')
    wlog.write('  display: inline-block;\n')
    wlog.write('  margin-left: 5px;\n')
    wlog.write('  font-weight: normal;\n')
    wlog.write('}\n')
    wlog.write('.total-stat {\n')
    wlog.write('  background-color: #0d6efd;\n')
    wlog.write('  color: white;\n')
    wlog.write('  padding: 2px 6px;\n')
    wlog.write('  border-radius: 4px;\n')
    wlog.write('  font-size: 0.9em;\n')
    wlog.write('}\n')
    wlog.write('.increase-stat {\n')
    wlog.write('  background-color: #dc3545;\n')
    wlog.write('  color: white;\n')
    wlog.write('  padding: 2px 6px;\n')
    wlog.write('  border-radius: 4px;\n')
    wlog.write('  font-size: 0.9em;\n')
    wlog.write('}\n')
    wlog.write('.plot-container {\n')
    wlog.write('  margin-bottom: 30px;\n')
    wlog.write('}\n')
    wlog.write('</style>\n')
    
    # Introduction text
    wlog.write('<div style="background-color: #cff4fc; padding: 15px; border-radius: 4px; margin-bottom: 20px;">\n')
    wlog.write('  <p>This page shows flag statistics for each processing step. Click on the plots for high-resolution versions.</p>\n')
    wlog.write('</div>\n')
    
    # Define the flag stats steps in the correct order
    flagstats_steps = [
        'run_importfits', 'flag_aoflagger', 'flag_apriori', 'flag_manual',
        'restore_flags', 'flag_manual_avg', 'bandpass', 'initial_gaincal',
        'applycal_all', 'flag_target'
    ]
    
    # Create navigation links
    wlog.write('<div class="sticky-nav">\n')
    wlog.write('  <strong>Jump to: </strong>\n')
    
    # Get list of steps that have plots
    available_steps = []
    for step in flagstats_steps:
        scan_plot = './weblog/plots/plots_flagstats/{}_flagstats_scans_{}.png'.format(msinfo['msfilename'], step)
        other_plot = './weblog/plots/plots_flagstats/{}_flagstats_other_{}.png'.format(msinfo['msfilename'], step)
        if os.path.isfile(scan_plot) or os.path.isfile(other_plot):
            available_steps.append(step)
    
    # Create navigation links for each available step
    for step in available_steps:
        wlog.write('<a class="nav-link" href="#{0}">{0}</a>\n'.format(step))
    
    wlog.write('</div>\n')
    
    # Process each step in the defined order
    prev_perc_flagged = 0.0
    for step in flagstats_steps:
        # Check if the flag stats file exists for this step
        flag_stats_file = './weblog/plots/plots_flagstats/flagstats_{}.pkl'.format(step)
        if not os.path.isfile(flag_stats_file):
            continue
            
        try:
            # Load the flag statistics for this step
            flag_stats = emutils.load_obj(flag_stats_file)
            # Calculate flagging percentages
            perc_flagged = flag_stats['flagged'] / flag_stats['total'] * 100.
            diff_flagged = perc_flagged - prev_perc_flagged
            prev_perc_flagged = perc_flagged
            
            # Define paths for the two plot types
            scan_plot_list = glob.glob('./weblog/plots/plots_flagstats/*_flagstats_scans_{}.png'.format(step))
            scan_plot = scan_plot_list[0] if scan_plot_list else None

            other_plot_list = glob.glob('./weblog/plots/plots_flagstats/*_flagstats_other_{}.png'.format(step))
            other_plot = other_plot_list[0] if other_plot_list else None
            
            # Check if at least one of the plots exists
            if scan_plot is not None or other_plot is not None:
                # Create a section for this step
                wlog.write('<div id="{0}" class="plot-container">\n'.format(step))
                
                # Step header with name and statistics (simple format as requested)
                wlog.write('<h3>\n')
                link_target = scan_plot if os.path.isfile(scan_plot) else other_plot
                #wlog.write('  <a href=".{0}" target="_blank">{1}</a>\n'.format(link_target, step))
                wlog.write(f'  {step}')
                wlog.write('  <span class="stats-label">\n')
                wlog.write('    (Total: <span class="total-stat">{0:3.1f}%</span>\n'.format(perc_flagged))
                wlog.write('    Increase: <span class="increase-stat">{0:3.1f}%</span>)\n'.format(diff_flagged))
                wlog.write('  </span>\n')
                wlog.write('</h3>\n')
                
                # Display the scan flags plot if it exists
                if os.path.isfile(scan_plot):
                    wlog.write('<a href=".{0}" target="_blank">\n'.format(scan_plot))
                    wlog.write('  <img style="max-width:1200px" src=".{0}">\n'.format(scan_plot))
                    wlog.write('</a><br>\n')
                
                # Display the other flags plot if it exists
                if os.path.isfile(other_plot):
                    wlog.write('<a href=".{0}" target="_blank">\n'.format(other_plot))
                    wlog.write('  <img style="max-width:1200px" src=".{0}">\n'.format(other_plot))
                    wlog.write('</a><br>\n')
                
                wlog.write('<hr>\n')
                wlog.write('</div>\n')
        except Exception as e:
            logger.warning('Error processing flag stats for step {}: {}'.format(step, e))
    
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
            
            # Create a modern table from the flag summary data
            wlog.write('    <div class="table-responsive">\n')
            wlog.write('      <table class="table table-sm table-striped">\n')
            
            # Parse header and data from the text file
            header_found = False
            
            for line in flag_data:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if not header_found:
                    # This is the header line
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
                    # This is a data line
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
    
    # Add floating back-to-top button
    wlog.write('<style>\n')
    wlog.write('.back-to-top {\n')
    wlog.write('  position: fixed;\n')
    wlog.write('  bottom: 20px;\n')
    wlog.write('  right: 20px;\n')
    wlog.write('  background-color: #007bff;\n')
    wlog.write('  color: white;\n')
    wlog.write('  padding: 10px 15px;\n')
    wlog.write('  border-radius: 4px;\n')
    wlog.write('  text-decoration: none;\n')
    wlog.write('  box-shadow: 0 2px 5px rgba(0,0,0,0.2);\n')
    wlog.write('}\n')
    wlog.write('.back-to-top:hover {\n')
    wlog.write('  background-color: #0056b3;\n')
    wlog.write('  color: white;\n')
    wlog.write('}\n')
    wlog.write('</style>\n')
    wlog.write('<a href="#top" class="back-to-top">↑ Top</a>\n')
    
    # Close the page
    weblog_foot(wlog)
    wlog.close()
