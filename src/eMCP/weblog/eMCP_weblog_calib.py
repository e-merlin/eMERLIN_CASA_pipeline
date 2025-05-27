import os
import glob
import logging
import pickle
import datetime
import numpy as np
from ..utils import eMCP_utils as emutils
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'
calib_dir = './weblog/calib/'

def write_caltable(caltable, wlog):
    """Display calibration table information"""
    wlog.write('<table class="table table-sm table-bordered">\n')
    
    for name, value in caltable.items():
        if isinstance(value, str):
            value = value.replace(',', ', ')
        wlog.write(f'<tr><td><strong>{name}</strong></td><td>{value}</td></tr>\n')
    
    wlog.write('</table>\n')

def write_fluxscale(wlog, msinfo):
    # Check for fluxscale plot file
    fluxscale_files = glob.glob(calib_dir + '{0}_fluxscale.png'.format(msinfo['msfilename']))
    if fluxscale_files:
        p = fluxscale_files[0]
        wlog.write(
            '<td><a href = ".{0}"><img style="max-width:700px" src=".{0}"></a></td>\n'
            .format(p))
    else:
        wlog.write('<td><p>No fluxscale plot available.</p></td>\n')

    # Check for fluxscale text file
    fluxscale_txt = calib_dir + 'allcal_ap.G1_fluxes.txt'
    if os.path.isfile(fluxscale_txt):
        with open(fluxscale_txt) as f:
            lines = f.readlines()
        wlog.write('\n<br><br>\n')
        wlog.write('CASA fluxscale output (not corrected by eMfactor):')
        wlog.write('\n<pre>')
        for line in lines:
            wlog.write(line)
        wlog.write('\n</pre>\n<br>')
    else:
        wlog.write('\n<br><br>\n')
        wlog.write('No fluxscale output available.')

def write_applycal_table(wlog, eMCP, ap_dict='applycal_dict'):
    """Display applycal information"""
    if ap_dict in eMCP:
        applycal_dict = eMCP[ap_dict]
        
        # Display the title based on which applycal dictionary we're using
        if ap_dict == 'applycal_dict_sp':
            wlog.write('<h4>Tables applied to narrow band data:</h4>\n')
        
        # Create a modern table with centering
        wlog.write('<div class="table-responsive">\n')
        wlog.write('<table class="table table-sm table-bordered" style="width:90%; margin: 0 auto;">\n')
        wlog.write('<thead>\n')
        wlog.write('<tr>\n')
        wlog.write('<th>Source</th>\n')
        wlog.write('<th>Table</th>\n')
        wlog.write('<th>Gainfield</th>\n')
        wlog.write('<th>spwmap</th>\n')
        wlog.write('<th>Interpolation</th>\n')
        wlog.write('</tr>\n')
        wlog.write('</thead>\n')
        wlog.write('<tbody>\n')
        
        for source, tables in applycal_dict.items():
            for i, table in enumerate(tables.keys()):
                if table == 'spwmap':
                    continue
                    
                wlog.write('<tr>\n')
                
                # Only show source name in first row for this source
                if i == 0:
                    wlog.write(f'<td>{source}</td>\n')
                else:
                    wlog.write('<td></td>\n')
                
                # Table name
                wlog.write(f'<td>{table}</td>\n')
                
                # Parameters (match original order: gainfield, spwmap, interp)
                if len(tables[table]) >= 2:
                    gainfield = tables[table][1]
                else:
                    gainfield = ''
                    
                if 'spwmap' in tables and table in tables['spwmap']:
                    spwmap = str(tables['spwmap'][table])
                else:
                    spwmap = ''
                    
                if len(tables[table]) >= 1:
                    interp = tables[table][0]
                else:
                    interp = ''
                
                wlog.write(f'<td>{gainfield}</td>\n')
                wlog.write(f'<td><code>{spwmap}</code></td>\n')
                wlog.write(f'<td>{interp}</td>\n')
                wlog.write('</tr>\n')
        
        wlog.write('</tbody>\n')
        wlog.write('</table>\n')
        wlog.write('</div>\n')
        wlog.write('<br>\n')

def weblog_calibration(eMCP):
    """Create a calibration page with modern styling"""
    # Get msinfo from eMCP
    msinfo = eMCP['msinfo']
    
    # Create the calibration page
    wlog = open(weblog_dir + "calibration.html", "w")
    weblog_header(wlog, 'Calibration', msinfo['run'])
    
    # Define all calibration steps as in the original
    all_calsteps = [
        'bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G0', 'bpcal.BP0', 'allcal_d.K1',
        'allcal_p.G1', 'allcal_ap.G1', 'fluxscale', 'bpcal.BP2', 'allcal_p.G3',
        'allcal_ap.G3', 'phscal_p_scan.G3', 'phscal_ap_scan.G3'
    ]
    
    # Add mixed mode steps if applicable
    if eMCP['is_mixed_mode']:
        all_calsteps.append('narrow_p_offset.G3')
        all_calsteps.append('narrow_bpcal.BP2')
        
    # Add sticky navigation bar
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
    wlog.write('</style>\n')
    
    wlog.write('<div class="sticky-nav">\n')
    wlog.write('  <strong>Jump to: </strong>\n')
    
    # Load caltables to get the proper names
    caltable_names = {}
    if os.path.isfile('./weblog/calib/caltables.yaml'):
        caltables = emutils.load_obj('./weblog/calib/caltables.yaml')
        for calstep in all_calsteps:
            if calstep in caltables:
                if calstep == 'fluxscale':
                    caltable_names[calstep] = 'Fluxscale'
                else:
                    caltable_names[calstep] = caltables[calstep]['name']
    
    # Add navigation links for each calibration step
    for calstep in all_calsteps:
        if calstep in caltable_names:
            wlog.write(f'<a class="nav-link" href="#{calstep}">{caltable_names[calstep]}</a>\n')
    
    # Add applycal link
    wlog.write('<a class="nav-link" href="#applycal">Applycal</a>\n')
    wlog.write('</div>\n')
    
    # Load caltables if available
    if os.path.isfile('./weblog/calib/caltables.yaml'):
        caltables = emutils.load_obj('./weblog/calib/caltables.yaml')
        
        # Process each calibration step
        for calstep in all_calsteps:
            logger.debug('calstep {}'.format(calstep))
            try:
                # Special handling for fluxscale
                if calstep == 'fluxscale':
                    if os.path.isfile(calib_dir + 'allcal_ap.G1_fluxes.txt'):
                        wlog.write(f'<div id="{calstep}">\n')
                        wlog.write('<h4 style="text-align:center">fluxscale</h4>\n')
                        write_fluxscale(wlog, msinfo)
                        wlog.write('</div>\n')
                        wlog.write('<hr>\n')
                # Process other calibration steps
                elif calstep in caltables:  # Only process if the calstep exists in caltables
                    wlog.write(f'<div id="{calstep}">\n')
                    wlog.write('<h4 style="text-align:center">{}</h4>\n'.format(caltables[calstep]['name']))
                    
                    # Create a table with caltable on left and plots on right
                    wlog.write('<table cellspacing="20" cellpadding="4px" style="width:90%; margin: 0 auto;">\n')
                    wlog.write('<tr>\n')
                    
                    # Left column for caltable
                    wlog.write('<td valign="top" style="width:40%">\n')
                    write_caltable(caltables[calstep], wlog)
                    wlog.write('</td>\n')
                    
                    # Right column for plots
                    all_plots = np.sort(glob.glob('./weblog/plots/caltables/*{}*.png'.format(calstep)))
                    
                    if len(all_plots) > 0:
                        for p in all_plots:
                            wlog.write('<td align="center">\n')
                            wlog.write('<a href=".{0}" target="_blank">\n'.format(p))
                            wlog.write('<img style="max-width:700px" src=".{0}" alt="{1}">\n'.format(p, os.path.basename(p)))
                            wlog.write('</a>\n')
                            wlog.write('</td>\n')
                    else:
                        wlog.write('<td><p>No plots available for this calibration step.</p></td>\n')
                    
                    wlog.write('</tr>\n')
                    wlog.write('</table>\n')
                    wlog.write('</div>\n')
                    wlog.write('<hr>\n')
            except Exception as e:
                logger.warning('Error processing calstep {}: {}'.format(calstep, e))
    
    # Applycal section
    try:
        wlog.write('<div id="applycal">\n')
        wlog.write('<h2 style="text-align:center">Applycal</h2>\n')
        wlog.write('<p style="text-align:center">List of tables and apply parameters used to correct each source:</p>\n')
        write_applycal_table(wlog, eMCP, ap_dict='applycal_dict')
        
        # Handle mixed mode data
        if eMCP['is_mixed_mode']:
            wlog.write('<p>Tables applied to narrow band data:</p>\n')
            write_applycal_table(wlog, eMCP, ap_dict='applycal_dict_sp')
        
        wlog.write('</div>\n')
    except Exception as e:
        logger.warning('Error in applycal section: {}'.format(e))
    
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
