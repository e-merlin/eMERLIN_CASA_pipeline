import os
import glob
import logging
import numpy as np
from .eMCP_weblog_modern import weblog_header, weblog_foot
from .eMCP_weblog import create_pnghtml_baselines

logger = logging.getLogger('logger')
weblog_dir = './weblog/'


def plots_data(msinfo, wlog):
    wlog.write('<div id="uncalibrated" class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Uncalibrated visibilities</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table table-sm table-bordered" style="width:80%; margin: 0 auto;">\n')
    wlog.write('      <tbody>\n')
    
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines(
            'plots_data', source,
            'Uncalibrated amplitude and phase against time and frequency.',
            msinfo, 'data')
        wlog.write('        <tr>\n')
        wlog.write('          <td><strong>{0}</strong></td>\n'.format(source))
        wlog.write('          <td><a href=".{0}" target="_blank" class="btn btn-primary btn-sm">View Plots</a></td>\n'.format(page_path))
        wlog.write('        </tr>\n')
    
    wlog.write('      </tbody>\n')
    wlog.write('      </table>\n')
    wlog.write('    </div>\n')
    wlog.write('  </div>\n')
    wlog.write('</div>\n')


def plots_corrected(msinfo, wlog):
    wlog.write('<div id="calibrated" class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Calibrated visibilities</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table table-sm table-bordered" style="width:80%; margin: 0 auto;">\n')
    wlog.write('      <tbody>\n')
    
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines(
            'plots_corrected', source,
            'Calibrated amplitude and phase against time and frequency.',
            msinfo, 'corrected')
        wlog.write('        <tr>\n')
        wlog.write('          <td><strong>{0}</strong></td>\n'.format(source))
        wlog.write('          <td><a href=".{0}" target="_blank" class="btn btn-primary btn-sm">View Plots</a></td>\n'.format(page_path))
        wlog.write('        </tr>\n')
    
    wlog.write('      </tbody>\n')
    wlog.write('      </table>\n')
    wlog.write('    </div>\n')
    wlog.write('  </div>\n')
    wlog.write('</div>\n')


def plots_uvplt(msinfo, wlog):
    wlog.write('<div id="uvplots" class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Calibrated UVplots</h3>\n')
    wlog.write('  <div>\n')
    all_plots = np.sort(glob.glob('./weblog/plots/plots_uvplt/*_uvplt_*png'))
    
    for p in all_plots:
        source_name = os.path.splitext(p)[0].split('_')[-1]
        wlog.write('    <div id="{0}-uvplot" class="subsection">\n'.format(source_name))
        wlog.write('      <h4 class="collapsible-header">{0}</h4>\n'.format(source_name))
        wlog.write('      <div>\n')
        
        wlog.write('<div class="row">\n')
        # Amplitude plot
        wlog.write('<div class="col-md-6 text-center">\n')
        wlog.write('<a href=".{0}" target="_blank">\n'.format(p))
        wlog.write('<img class="img-fluid" style="max-width:100%" src=".{0}" alt="Amplitude UVplot for {1}">\n'.format(p, source_name))
        wlog.write('</a>\n')
        wlog.write('</div>\n')
        
        # Phase plot
        wlog.write('<div class="col-md-6 text-center">\n')
        wlog.write('<a href=".{0}" target="_blank">\n'.format(p.replace('_a_', '_p_')))
        wlog.write('<img class="img-fluid" style="max-width:100%" src=".{0}" alt="Phase UVplot for {1}">\n'.format(p.replace('_a_', '_p_'), source_name))
        wlog.write('</a>\n')
        wlog.write('</div>\n')
        wlog.write('</div>\n')
        
        # Model plots for calibration sources
        if source_name in msinfo['sources']['calsources'].split(','):
            p_model = './weblog/plots/plots_uvplt/{0}_uvpltmodel_a_{1}.png'.format(
                msinfo['msfilename'], source_name)
                
            wlog.write('<div class="row mt-3">\n')
            # Model amplitude plot
            wlog.write('<div class="col-md-6 text-center">\n')
            wlog.write('<a href=".{0}" target="_blank">\n'.format(p_model))
            wlog.write('<img class="img-fluid" style="max-width:100%" src=".{0}" alt="Model Amplitude UVplot for {1}">\n'.format(p_model, source_name))
            wlog.write('</a>\n')
            wlog.write('</div>\n')
            
            # Model phase plot
            wlog.write('<div class="col-md-6 text-center">\n')
            wlog.write('<a href=".{0}" target="_blank">\n'.format(p_model.replace('_a_', '_p_')))
            wlog.write('<img class="img-fluid" style="max-width:100%" src=".{0}" alt="Model Phase UVplot for {1}">\n'.format(p_model.replace('_a_', '_p_'), source_name))
            wlog.write('</a>\n')
            wlog.write('</div>\n')
            wlog.write('</div>\n')
        
        wlog.write('      </div>\n')
        wlog.write('    </div>\n')
    
    wlog.write('  </div>\n')
    wlog.write('</div>\n')


def weblog_plots(weblog, weblog_link, plots_link, msinfo):
    
    # Create the plots page
    wlog = open(weblog_dir + "plots.html", "w")
    weblog_header(wlog, 'Plots', msinfo['run'])
    
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
    
    # Add navigation links
    plot_sections = []
    
    if (os.path.isdir('./weblog/plots/plots_data/')) and \
       (os.listdir('./weblog/plots/plots_data/')):
        plot_sections.append(('uncalibrated', 'Uncalibrated Visibilities'))
        
    if (os.path.isdir('./weblog/plots/plots_corrected/')) and \
       (os.listdir('./weblog/plots/plots_corrected/')):
        plot_sections.append(('calibrated', 'Calibrated Visibilities'))
        
    if (os.path.isdir('./weblog/plots/plots_uvplt/')) and \
       (os.listdir('./weblog/plots/plots_uvplt/')):
        plot_sections.append(('uvplots', 'UV Plots'))
    
    # Create navigation links
    for section_id, section_name in plot_sections:
        wlog.write('<a class="nav-link" href="#{0}">{1}</a>\n'.format(section_id, section_name))
    
    wlog.write('</div>\n')
    
    # Display available plots
    if (os.path.isdir('./weblog/plots/plots_data/')) and \
       (os.listdir('./weblog/plots/plots_data/')):
        plots_data(msinfo, wlog)
        
    if (os.path.isdir('./weblog/plots/plots_corrected/')) and \
       (os.listdir('./weblog/plots/plots_corrected/')):
        plots_corrected(msinfo, wlog)
        
    if (os.path.isdir('./weblog/plots/plots_uvplt/')) and \
       (os.listdir('./weblog/plots/plots_uvplt/')):
        plots_uvplt(msinfo, wlog)
    
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
