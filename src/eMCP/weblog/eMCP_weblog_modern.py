import os
import numpy as np
import glob
import datetime
from ..utils import eMCP_utils as emutils

import logging

logger = logging.getLogger('logger')

weblog_dir = './weblog/'
info_dir = './weblog/info/'
calib_dir = './weblog/calib/'
plots_dir = './weblog/plots/'
logs_dir = './logs/'
images_dir = './weblog/images/'
flagstats_dir = './weblog/flagstats/'

weblog_link = './'
info_link = './info/'
calib_link = './calib/'
plots_link = './plots/'
images_link = './images/'

line0 = '-' * 15

def weblog_nav_item(weblog_link, name, link, active=False):
    """Creates a navigation item for the modern weblog."""
    active_class = ' active' if active else ''
    return f'<li class="nav-item"><a class="nav-link{active_class}" href="{weblog_link}{link}.html">{name}</a></li>\n'

def weblog_header(wlog, section, project):
    """Creates a modern header with responsive navigation for the weblog."""
    wlog.write('<!DOCTYPE html>\n')
    wlog.write('<html lang="en">\n')
    wlog.write('<head>\n')
    wlog.write('  <meta charset="UTF-8">\n')
    wlog.write('  <meta name="viewport" content="width=device-width, initial-scale=1.0">\n')
    wlog.write(f'  <title>{project} - eMCP</title>\n')
    wlog.write(f'  <link rel="stylesheet" type="text/css" href="{weblog_link}eMCP_modern.css"/>\n')
    wlog.write(f'  <link rel="icon" href="{weblog_link}eMCP_logo.png">\n')
    wlog.write('</head>\n')
    wlog.write('<body>\n')
    
    # Header with logo
    wlog.write('  <header class="header">\n')
    wlog.write('    <div class="container">\n')
    wlog.write('      <div class="header-content">\n')
    wlog.write('        <div class="header-logo">\n')
    wlog.write(f'          <a href="./index.html"><img src="{weblog_link}eMCP_logo.png" alt="eMCP Logo"></a>\n')
    wlog.write(f'          <h1 class="header-title">e-MERLIN Pipeline Web Log<br>{project}</h1>\n')
    wlog.write('        </div>\n')
    wlog.write('      </div>\n')
    wlog.write('    </div>\n')
    wlog.write('  </header>\n')
    
    # Navigation
    wlog.write('  <nav class="nav">\n')
    wlog.write('    <div class="container">\n')
    wlog.write('      <ul class="nav-list">\n')
    
    # Navigation items
    nav_items = [
        ('Home', 'index'),
        ('Observation summary', 'obs_summary'),
        ('Pipeline info', 'pipelineinfo'),
        ('Calibration', 'calibration'),
        ('Plots', 'plots'),
        ('Flag statistics', 'flagstats'),
        ('Images', 'images'),
        ('Download data', 'download')
    ]
    
    for name, link in nav_items:
        is_active = name.lower() == section.lower() or (section == 'Home' and link == 'index')
        wlog.write(weblog_nav_item(weblog_link, name, link, is_active))
    
    wlog.write('      </ul>\n')
    wlog.write('    </div>\n')
    wlog.write('  </nav>\n')
    
    # Main content container
    wlog.write('  <main class="main">\n')
    wlog.write('    <div class="container">\n')
    wlog.write(f'      <section class="section">\n')

def weblog_foot(wlog):
    """Creates a modern footer for the weblog."""
    wlog.write('    </section>\n')  # Close the main section
    wlog.write('  </div>\n')  # Close the container
    
    # Footer
    wlog.write('  <footer class="footer">\n')
    wlog.write('    <div class="container">\n')
    wlog.write('      <div class="footer-content">\n')
    wlog.write('        <div class="footer-logo">\n')
    wlog.write(f'          <a href="http://www.e-merlin.ac.uk/" target="_blank"><img src="{weblog_link}emerlin-2.gif" alt="e-MERLIN"></a>\n')
    wlog.write('        </div>\n')
    wlog.write('        <ul class="footer-links">\n')
    
    links = {
        'e-MERLIN Pipeline Github': 'https://github.com/e-merlin/eMERLIN_CASA_pipeline',
        'e-MERLIN': 'http://www.e-merlin.ac.uk/',
        'User support and data reduction for e-MERLIN': 'http://www.e-merlin.ac.uk/data_red/',
        'CASA': 'https://casa.nrao.edu/'
    }
    
    for name, link in links.items():
        wlog.write(f'          <li><a href="{link}" target="_blank">{name}</a></li>\n')
    
    wlog.write('        </ul>\n')
    wlog.write('      </div>\n')
    wlog.write('    </div>\n')
    wlog.write('  </footer>\n')
    
    # JavaScript for interactivity
    wlog.write(f'  <script src="{weblog_link}eMCP_modern.js"></script>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext, text='txt'):
    """Creates a link to a text file."""
    wlog.write(f'<p>{intext}: <a href="{infile}" target="_blank">{text}</a></p>\n')

def create_data_summary(data_dict, title=None):
    """Creates a data summary box with key-value pairs."""
    html = '<div class="data-summary">\n'
    if title:
        html += f'  <h3>{title}</h3>\n'
    html += '  <dl>\n'
    
    for key, value in data_dict.items():
        html += f'    <dt>{key}</dt>\n'
        html += f'    <dd>{value}</dd>\n'
    
    html += '  </dl>\n'
    html += '</div>\n'
    
    return html

def weblog_index(msinfo):
    """Create a home page for the weblog in table format like the original"""
    wlog = open("./weblog/index.html", "w")
    weblog_header(wlog, 'Home', msinfo['run'])
    
    # Main info section with table format like the original
    wlog.write('<div class="section centered">\n')
    wlog.write('  <h2 class="section-title">Project Information</h2>\n')
    
    # Display key information in table format like the original
    wlog.write('  <table class="table" style="width:50%; margin: 0 auto;">\n')
    
    # Project info
    wlog.write(f'    <tr><td>Project</td><td>{msinfo.get("project", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>Run</td><td>{msinfo.get("run", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>MS file</td><td>{msinfo.get("msfile", "Unknown")}</td></tr>\n')
    
    # Time information
    if 't_ini' in msinfo and 't_end' in msinfo:
        wlog.write(f'    <tr><td>Start</td><td>{msinfo["t_ini"].strftime("%Y-%m-%d %H:%M")}</td></tr>\n')
        wlog.write(f'    <tr><td>End</td><td>{msinfo["t_end"].strftime("%Y-%m-%d %H:%M")}</td></tr>\n')
    
    # Observation details
    wlog.write(f'    <tr><td>Band</td><td>{msinfo.get("band", "Unknown")}</td></tr>\n')
    
    if 'antennas' in msinfo:
        wlog.write(f'    <tr><td>Antennas</td><td>{", ".join(msinfo["antennas"])}</td></tr>\n')
    
    if 'sources' in msinfo and 'mssources' in msinfo['sources']:
        num_sources = len(msinfo['sources']['mssources'].split(','))
        wlog.write(f'    <tr><td>Number of sources</td><td>{num_sources}</td></tr>\n')
    
    wlog.write(f'    <tr><td>Integration time</td><td>{msinfo.get("int_time", "Unknown")}s</td></tr>\n')
    
    # Frequency information
    if 'freq_ini' in msinfo and 'freq_end' in msinfo:
        wlog.write(f'    <tr><td>Frequency</td><td>{msinfo["freq_ini"]:5.2f} - {msinfo["freq_end"]:5.2f} GHz</td></tr>\n')
    
    wlog.write(f'    <tr><td>Num. spw</td><td>{msinfo.get("num_spw", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>Channels/spw</td><td>{msinfo.get("nchan", "Unknown")}</td></tr>\n')
    
    if 'chan_res' in msinfo:
        chan_width = msinfo["chan_res"] * 1000  # Convert to MHz
        wlog.write(f'    <tr><td>Channel width</td><td>{chan_width:6.2f} MHz</td></tr>\n')
        
        if 'nchan' in msinfo:
            spw_bw = chan_width * msinfo['nchan']
            wlog.write(f'    <tr><td>spw bandwidth</td><td>{spw_bw:6.0f} MHz</td></tr>\n')
            
            if 'num_spw' in msinfo:
                total_bw = spw_bw * msinfo['num_spw']
                wlog.write(f'    <tr><td>Total bandwidth</td><td>{total_bw:6.0f} MHz</td></tr>\n')
    
    wlog.write(f'    <tr><td>Polarizations</td><td>{msinfo.get("polarizations", "Unknown")}</td></tr>\n')
    wlog.write('  </table>\n')
    
    # Add notes if available
    notes_file = './{0}.notes.txt'.format(msinfo['run'])
    if os.path.isfile(notes_file):
        wlog.write('  <div class="subsection centered">\n')
        wlog.write('    <h3>Notes and Comments</h3>\n')
        write_link_txt(wlog, notes_file, 'Notes and observing comments')
        wlog.write('  </div>\n')
    
    wlog.write('</div>\n')
    
    weblog_foot(wlog)
    wlog.close()

def weblog_obssum(msinfo):
    """Creates a modernized observation summary page."""
    wlog = open(weblog_dir + "obs_summary.html", "w")
    weblog_header(wlog, 'Observation summary', msinfo['run'])
    
    # Summary section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Summary</h3>\n')
    wlog.write('  <div>\n')
    
    # Listobs link
    listobs_file = os.path.join(info_link, msinfo['msfile'] + '.listobs.txt')
    write_link_txt(wlog, listobs_file, 'Summary of current observation (listobs)')
    
    # Other listobs files if available
    all_listobs = [
        os.path.basename(l)
        for l in glob.glob(info_dir + msinfo['run'] + '*listobs.txt')
    ]
    try:
        all_listobs.remove(os.path.basename(listobs_file))
    except ValueError:
        pass
    
    if len(all_listobs) > 0:
        wlog.write('<p>Other available listobs files:</p>\n')
        wlog.write('<div class="links-list">\n')
        for listobs in all_listobs:
            write_link_txt(wlog, info_link + listobs, os.path.basename(listobs))
        wlog.write('</div>\n')
    
    wlog.write('  </div>\n')
    wlog.write('</div>\n')
    
    # Sources section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Sources</h3>\n')
    wlog.write('  <div>\n')
    
    # Source separation table
    sepfile = info_link + 'source_separations.txt'
    if msinfo['sources']['targets'] != '':
        wlog.write('<div class="table-responsive">\n')
        wlog.write('  <table class="table">\n')
        wlog.write('    <thead>\n')
        wlog.write('      <tr>\n')
        wlog.write('        <th>Target</th>\n')
        wlog.write('        <th>Phase cal</th>\n')
        wlog.write('        <th>Separation [deg]</th>\n')
        wlog.write('      </tr>\n')
        wlog.write('    </thead>\n')
        wlog.write('    <tbody>\n')
        
        separations = msinfo['separations']
        for s1, s2 in zip(msinfo['sources']['targets'].split(','),
                         msinfo['sources']['phscals'].split(',')):
            wlog.write('      <tr>\n')
            
            if s1 not in msinfo['sources']['mssources'].split(','):
                logger.warning('{} not in MS'.format(s1))
                separation = 'Not in MS'
            elif s2 not in msinfo['sources']['mssources'].split(','):
                logger.warning('{} not in MS'.format(s2))
                separation = 'Not in MS'
            else:
                try:
                    # Match the original format exactly for the separations
                    separation = '{0:5.2f}'.format(separations[s1 + '-' + s2])
                except:
                    try:
                        separation = '{0:5.2f}'.format(separations[s2 + '-' + s1])
                    except:
                        separation = '0.0'
            
            wlog.write(f'        <td>{s1}</td>\n')
            wlog.write(f'        <td>{s2}</td>\n')
            wlog.write(f'        <td>{separation}</td>\n')
            wlog.write('      </tr>\n')
        
        wlog.write('    </tbody>\n')
        wlog.write('  </table>\n')
        wlog.write('</div>\n')
        
        # Add link to the separations file
        if os.path.isfile('./weblog/info/source_separations.txt'):
            write_link_txt(wlog, sepfile, 'View all source separations')
    else:
        wlog.write('<p>No target sources found</p>\n')
    
    wlog.write('  </div>\n')
    wlog.write('</div>\n')
    
    # Sources in MS section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Sources in MS</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table" style="width:80%; margin: 0 auto;">\n')
    wlog.write('        <thead>\n')
    wlog.write('          <tr>\n')
    wlog.write('            <th>Source</th>\n')
    wlog.write('            <th>Intent</th>\n')
    wlog.write('            <th>Coordinates (h:m:s d:m:s)</th>\n')
    wlog.write('            <th>MJD range</th>\n')
    wlog.write('          </tr>\n')
    wlog.write('        </thead>\n')
    wlog.write('        <tbody>\n')
    
    # Add source information like in the original
    if 'sources' in msinfo and 'mssources' in msinfo['sources']:
        for source in msinfo['sources']['mssources'].split(','):
            if 'source_timerange_mjd' in msinfo['sources'] and source in msinfo['sources']['source_timerange_mjd']:
                mjd_ini, mjd_end = msinfo['sources']['source_timerange_mjd'][source]
                mjd_range = f"{mjd_ini:11.5f} - {mjd_end:11.5f}"
            else:
                mjd_range = "Unknown"
                
            intent = msinfo['sources'].get('source_intent', {}).get(source, 'Unknown')
            
            # Get coordinates from directions in msinfo if available
            if 'directions' in msinfo and source in msinfo['directions']:
                coords = msinfo['directions'][source]
            else:
                coords = "Unknown"
            
            wlog.write(f'          <tr>\n')
            wlog.write(f'            <td>{source}</td>\n')
            wlog.write(f'            <td>{intent}</td>\n')
            wlog.write(f'            <td>{coords}</td>\n')
            wlog.write(f'            <td>{mjd_range}</td>\n')
            wlog.write('          </tr>\n')
    
    wlog.write('        </tbody>\n')
    wlog.write('      </table>\n')
    wlog.write('    </div>\n')
    
    # Add missing sources info
    if 'sources' in msinfo and 'allsources' in msinfo['sources'] and 'mssources' in msinfo['sources']:
        missing_sources = ', '.join([
            s for s in msinfo['sources']['allsources'].split(',')
            if s not in msinfo['sources']['mssources']
        ])
        
        if missing_sources != '':
            wlog.write(f'    <p><em>Sources specified in the inputs file but not in MS: {missing_sources}</em></p>\n')
        else:
            wlog.write('    <p><em>All sources specified in the inputs file are in the MS.</em></p>\n')
    
    wlog.write('  </div>\n')
    wlog.write('</div>\n')
    
    # Antennas section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Antennas</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <p class="centered">Reference antenna: {}</p>\n'.format(msinfo.get('refant', 'Not specified')))
    
    # List antennas in a centered table
    if 'antennas' in msinfo:
        wlog.write('    <table class="table" style="width:50%; margin: 0 auto;">\n')
        wlog.write('      <tr><th>Antennas</th></tr>\n')
        for ant in msinfo['antennas']:
            wlog.write(f'      <tr><td style="text-align:center;">{ant}</td></tr>\n')
        wlog.write('    </table>\n')
    wlog.write('  </div>\n')
    
    wlog.write('</div>\n')
    
    # Source elevation section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Source elevation</h3>\n')
    wlog.write('  <div>\n')
    
    # Look for elevation plots using the msfilename pattern like the original
    if 'msfilename' in msinfo:
        elev_plots = glob.glob(f'./weblog/plots/plots_observation/{msinfo["msfilename"]}_elevation.png')
        if elev_plots:
            plot_path = elev_plots[0]
            rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
            wlog.write('    <div class="image-container">\n')
            wlog.write(f'      <img src="{rel_path}" class="centered" style="max-width:700px;" alt="Source elevation plot">\n')
            wlog.write('    </div>\n')
        else:
            # Try to find any elevation plot as fallback
            all_elev_plots = glob.glob('./weblog/plots/plots_observation/*elevation*.png')
            if all_elev_plots:
                plot_path = all_elev_plots[0]
                rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
                wlog.write('    <div class="image-container">\n')
                wlog.write(f'      <img src="{rel_path}" class="centered" style="max-width:700px;" alt="Source elevation plot">\n')
                wlog.write('    </div>\n')
            else:
                wlog.write('    <p class="centered">No elevation plot available</p>\n')
    else:
        wlog.write('    <p class="centered">No elevation plot available - msfilename not specified</p>\n')
    wlog.write('  </div>\n')
    
    wlog.write('</div>\n')
    
    # UV coverage section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3>UV coverage</h3>\n')
    
    # Look for UV coverage plots using the msfilename pattern like the original
    if 'msfilename' in msinfo:
        uvcov_plots = sorted(glob.glob(f'./weblog/plots/plots_observation/{msinfo["msfilename"]}_uvcov_*.png'))
        if uvcov_plots:
            wlog.write('  <div class="image-container">\n')
            for plot_path in uvcov_plots:
                rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
                wlog.write(f'    <img src="{rel_path}" class="centered" style="max-width:700px; margin-bottom:20px;" alt="UV coverage plot">\n')
            wlog.write('  </div>\n')
        else:
            # Try to find any UV plot as fallback
            all_uv_plots = glob.glob('./weblog/plots/plots_observation/*uv*.png') + glob.glob('./weblog/plots/plots_observation/*UV*.png')
            if all_uv_plots:
                wlog.write('  <div class="image-container">\n')
                for plot_path in sorted(all_uv_plots):
                    rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
                    wlog.write(f'    <img src="{rel_path}" class="centered" style="max-width:700px; margin-bottom:20px;" alt="UV coverage plot">\n')
                wlog.write('  </div>\n')
            else:
                wlog.write('  <p class="centered">No UV coverage plots available</p>\n')
    else:
        wlog.write('  <p class="centered">No UV coverage plots available - msfilename not specified</p>\n')
    
    wlog.write('</div>\n')
    
    weblog_foot(wlog)
    wlog.close()
