import os
import glob
import logging
import numpy as np
from ..utils import eMCP_utils as emutils
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'
calib_dir = './weblog/calib/'

def write_caltable(caltable, wlog):
    """Display calibration table information"""
    wlog.write('<table class="table table-sm" style="width: 100%; max-width: 580px; ">\n')
    for name, value in caltable.items():
        if isinstance(value, str):
            value = value.replace(',', ', ')
        wlog.write(f'<tr><td><strong>{name}</strong></td><td>{value}</td></tr>\n')
    wlog.write('</table>\n')

def write_fluxscale(wlog, msinfo):
    fluxscale_files = glob.glob(calib_dir + '{0}_fluxscale.png'.format(msinfo['msfilename']))
    if fluxscale_files:
        p = fluxscale_files[0]
        wlog.write(
            '<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>\n'
            .format(p))
    else:
        wlog.write('<td><p>No fluxscale plot available.</p></td>\n')
    fluxscale_txt = calib_dir + 'allcal_ap.G1_fluxes.txt'
    if os.path.isfile(fluxscale_txt):
        with open(fluxscale_txt) as f:
            lines = f.readlines()
        wlog.write('\n<br><br>\n')
        wlog.write('CASA fluxscale output (not corrected by eMfactor):')
        wlog.write('\n<pre style="text-align: left; max-width: 700px; font-size: 0.96em;">')
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
        if ap_dict == 'applycal_dict_sp':
            wlog.write('<h4>Tables applied to narrow band data:</h4>\n')
        wlog.write('<div class="table-responsive">\n')
        wlog.write('<table class="table table-sm" style="width: 98%; max-width: 1000px; ">\n')
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
                if i == 0:
                    wlog.write(f'<td>{source}</td>\n')
                else:
                    wlog.write('<td></td>\n')
                wlog.write(f'<td>{table}</td>\n')
                gainfield = tables[table][1] if len(tables[table]) >= 2 else ''
                spwmap = str(tables['spwmap'][table]) if 'spwmap' in tables and table in tables['spwmap'] else ''
                interp = tables[table][0] if len(tables[table]) >= 1 else ''
                wlog.write(f'<td>{gainfield}</td>\n')
                wlog.write(f'<td><code>{spwmap}</code></td>\n')
                wlog.write(f'<td>{interp}</td>\n')
                wlog.write('</tr>\n')
        wlog.write('</tbody>\n')
        wlog.write('</table>\n')
        wlog.write('</div>\n')
        wlog.write('<br>\n')

def weblog_calibration(eMCP):
    """Create a calibration page with sidebar layout, sticky Jump-to and wide body."""
    msinfo = eMCP['msinfo']

    wlog = open(weblog_dir + "calibration.html", "w")
    weblog_header(wlog, 'Calibration', msinfo.get('run', msinfo.get('project', 'eMERLIN')))

    # Calibration steps
    all_calsteps = [
        'bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G0', 'bpcal.BP0', 'allcal_d.K1',
        'allcal_p.G1', 'allcal_ap.G1', 'fluxscale', 'bpcal.BP2', 'allcal_p.G3',
        'allcal_ap.G3', 'phscal_p_scan.G3', 'phscal_ap_scan.G3'
    ]
    if eMCP.get('is_mixed_mode', False):
        all_calsteps += ['narrow_p_offset.G3', 'narrow_bpcal.BP2']

    # Load caltable names
    caltable_names = {}
    if os.path.isfile('./weblog/calib/caltables.yaml'):
        caltables = emutils.load_obj('./weblog/calib/caltables.yaml')
        for calstep in all_calsteps:
            if calstep in caltables:
                caltable_names[calstep] = caltables[calstep]['name'] if calstep != 'fluxscale' else 'Fluxscale'
    
    wlog.write('<div class="calib-layout">\n')
    wlog.write('<div class="calib-main">\n')
    
    # --- Calibration Steps Content ---
    if os.path.isfile('./weblog/calib/caltables.yaml'):
        caltables = emutils.load_obj('./weblog/calib/caltables.yaml')
        for calstep in all_calsteps:
            logger.debug('calstep {}'.format(calstep))
            try:
                # Special handling for fluxscale
                if calstep == 'fluxscale':
                    if os.path.isfile(calib_dir + 'allcal_ap.G1_fluxes.txt'):
                        wlog.write(f'<div id="{calstep}" class="subsection">\n')
                        wlog.write('  <h3 class="5ollapsible-header">Fluxscale</h3>\n')
                        wlog.write('  <div>\n')
                        write_fluxscale(wlog, msinfo)
                        wlog.write('  </div>\n')
                        wlog.write('</div>\n')
                elif calstep in caltables:
                    wlog.write(f'<div id="{calstep}" class="subsection">\n')
                    wlog.write('  <h3 class="collapsible-header">{}</h3>\n'.format(caltables[calstep]['name']))
                    wlog.write('  <div>\n')
                    wlog.write('    <table cellspacing="18" cellpadding="4px" style="width: 100%; max-width: 1000px; ">\n')
                    wlog.write('    <tr>\n')
                    # Left column for caltable
                    wlog.write('<td valign="top" style="width:44%; min-width:220px;">\n')
                    write_caltable(caltables[calstep], wlog)
                    wlog.write('</td>\n')
                    # Right column for plots
                    all_plots = np.sort(glob.glob('./weblog/plots/caltables/*{}*.png'.format(calstep)))
                    if len(all_plots) > 0:
                        for p in all_plots:
                            wlog.write('<td align="center">\n')
                            wlog.write('<a href=".{0}" target="_blank">\n'.format(p))
                            wlog.write('<img style="max-width:500px; margin-bottom:8px;" src=".{0}" alt="{1}">\n'.format(p, os.path.basename(p)))
                            wlog.write('</a>\n')
                            wlog.write('</td>\n')
                    else:
                        wlog.write('<td><p>No plots available for this calibration step.</p></td>\n')
                    wlog.write('    </tr>\n')
                    wlog.write('    </table>\n')
                    wlog.write('  </div>\n')
                    wlog.write('</div>\n')
            except Exception as e:
                logger.warning('Error processing calstep {}: {}'.format(calstep, e))
    
    # --- Applycal section ---
    try:
        wlog.write('<div id="applycal" class="section">\n')
        wlog.write('  <h2 class="section-title">Applycal</h2>\n')
        wlog.write('  <p style="text-align:center">List of tables and apply parameters used to correct each source:</p>\n')
        write_applycal_table(wlog, eMCP, ap_dict='applycal_dict')
        if eMCP.get('is_mixed_mode', False):
            wlog.write('<p>Tables applied to narrow band data:</p>\n')
            write_applycal_table(wlog, eMCP, ap_dict='applycal_dict_sp')
        wlog.write('</div>\n')
    except Exception as e:
        logger.warning('Error in applycal section: {}'.format(e))
    
    wlog.write('</div>')  # close calib-main
    
    # --- Sticky Sidebar with Jump Links ---
    wlog.write('<nav class="calib-jump-sidebar">\n')
    wlog.write('<h4>Jump to section</h4>\n')
    wlog.write('<div class="calib-jump-links">\n')
    for calstep in all_calsteps:
        if calstep in caltable_names:
            wlog.write(f'<a class="calib-jump-link" href="#{calstep}">{caltable_names[calstep]}</a>\n')
    wlog.write('<a class="calib-jump-link" href="#applycal">Applycal</a>\n')
    wlog.write('</div>\n')
    wlog.write('</nav>\n')
    
    wlog.write('</div>')  # close calib-layout

    weblog_foot(wlog)
    wlog.close()
