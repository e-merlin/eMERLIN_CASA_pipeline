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
    weblog_link = './'
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
    wlog.write(f'  <div class="project-title">{project}</div>\n')
    wlog.write('  <section class="section">\n')


def weblog_foot(wlog):
    """Creates a modern footer for the weblog."""
    weblog_link = './'
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
    wlog.write('<a href="#top" class="back-to-top">↑ Top</a>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext, text='txt'):
    """Creates a link to a text file."""
    wlog.write(f'<p>{intext}: <a href="{infile}" target="_blank">{text}</a></p>\n')

