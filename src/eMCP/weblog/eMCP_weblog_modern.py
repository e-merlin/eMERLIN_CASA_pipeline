import os
import numpy as np
import glob
import datetime
from ..utils import eMCP_utils as emutils
from ..utils import eMCP_paths as empaths
import html

import logging

logger = logging.getLogger('logger')

weblog_dir = empaths.WEBLOG_DIR
info_dir = empaths.INFO_DIR
calib_dir = empaths.CALIB_DIR
plots_dir = empaths.PLOTS_DIR
logs_dir = empaths.LOGS_DIR
images_dir = empaths.IMAGES_DIR
flagstats_dir = empaths.FLAGSTATS_DIR

weblog_link = empaths.WEBLOG_LINK
info_link = empaths.INFO_LINK
calib_link = empaths.CALIB_LINK
plots_link = empaths.PLOTS_LINK
images_link = empaths.IMAGES_LINK

line0 = '-' * 15

def weblog_nav_item(weblog_link, name, link, active=False):
    """Creates a navigation item for the modern weblog."""
    active_class = ' active' if active else ''
    safe_name = html.escape(str(name))
    safe_href = html.escape(f'{weblog_link}{link}.html', quote=True)
    return f'<li class="nav-item"><a class="nav-link{active_class}" href="{safe_href}">{safe_name}</a></li>\n'

def weblog_header(wlog, section, project):
    weblog_link = empaths.WEBLOG_LINK
    safe_project = html.escape(str(project))
    wlog.write('<!DOCTYPE html>\n')
    wlog.write('<html lang="en">\n')
    wlog.write('<head>\n')
    wlog.write('  <meta charset="UTF-8">\n')
    wlog.write('  <meta name="viewport" content="width=device-width, initial-scale=1.0">\n')
    wlog.write(f'  <title>{safe_project} - eMCP</title>\n')
    wlog.write(f'  <link rel="stylesheet" type="text/css" href="{weblog_link}eMCP_modern.css"/>\n')
    wlog.write(f'  <link rel="icon" href="{weblog_link}eMCP_logo.png">\n')
    wlog.write('</head>\n')
    wlog.write('<body>\n')

    wlog.write('<div class="sidebar">\n')
    wlog.write('  <div class="header-logo" style="text-align:center; margin-bottom:1em;">\n')
    wlog.write(f'    <a href="./index.html"><img src="{weblog_link}eMCP_logo.png" alt="eMCP Logo" style="max-width:100px;"></a>\n')
#    wlog.write(f'    <h2 class="header-title" style="font-size:1.1em; margin:0.5em 0 0 0;">{project}</h2>\n')
    wlog.write('  </div>\n')
    wlog.write('  <ul class="nav-list">\n')
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

    wlog.write('  </ul>\n')
    wlog.write('</div>\n')

    # Main content container to the right of the sidebar
    wlog.write('<main class="main">\n')
    wlog.write(f'  <div class="project-title">{safe_project}</div>\n')
    wlog.write('  <section class="section">\n')


def weblog_foot(wlog):
    """Creates a modern footer for the weblog."""
    weblog_link = empaths.WEBLOG_LINK
    wlog.write('    </section>\n')  # Close the main section
    wlog.write('  </main>\n')  # Close the main content container

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
    wlog.write('<a href="#top" class="back-to-top">↑ Top</a>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext, text='txt'):
    """Creates a link to a text file."""
    safe_intext = html.escape(str(intext))
    safe_infile = html.escape(str(infile), quote=True)
    safe_text = html.escape(str(text))
    wlog.write(f'<p>{safe_intext}: <a href="{safe_infile}" target="_blank">{safe_text}</a></p>\n')
