import os
import glob
import logging
import numpy as np
from html import escape
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'


def create_pnghtml_baselines(plots_path, source, subtitle, msinfo, datacolumn):
    page_path = weblog_dir + plots_path + '_' + source + ".html"
    plot_base = weblog_dir + 'plots/' + plots_path + '/{0}_4plot_{1}_{2}'.format(
        msinfo['msfilename'], source, datacolumn)
    plot_labels = (
        ('Amp vs Time', plot_base + '0.png'),
        ('Phase vs Time', plot_base + '1.png'),
        ('Amp vs Freq', plot_base + '2.png'),
        ('Phase vs Freq', plot_base + '3.png'),
    )

    wlog = open(page_path, "w")
    weblog_header(wlog, 'Plots', msinfo['run'])
    wlog.write(f'<div class="subsection">\n')
    wlog.write(f'  <h3 class="collapsible-header">{escape(source)}</h3>\n')
    wlog.write(f'  <p>{escape(subtitle)}</p>\n')
    wlog.write('  <div class="visibility-plot-grid">\n')
    for label, plot_path in plot_labels:
        wlog.write('    <div class="visibility-plot-cell">\n')
        wlog.write(f'      <h4>{label}</h4>\n')
        if os.path.isfile(plot_path):
            wlog.write(f'      <a href=".{plot_path}" target="_blank">\n')
            wlog.write(f'        <img src=".{plot_path}" alt="{label} for {escape(source)}">\n')
            wlog.write('      </a>\n')
        else:
            wlog.write(f'      <div class="missing-plot">Missing plot: {escape(os.path.basename(plot_path))}</div>\n')
        wlog.write('    </div>\n')
    wlog.write('  </div>\n')
    wlog.write('</div>\n')
    weblog_foot(wlog)
    wlog.close()
    return page_path


def plots_data(msinfo, wlog):
    wlog.write('<div id="uncalibrated" class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Uncalibrated visibilities</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table table-sm table-bordered" style="width:80%;">\n')
    wlog.write('      <tbody>\n')
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines(
            'plots_data', source,
            'Uncalibrated amplitude and phase against time and frequency.',
            msinfo, 'data')
        wlog.write(f'        <tr><td><strong>{source}</strong></td>')
        wlog.write(f'<td><a href=".{page_path}" target="_blank" class="btn btn-primary btn-sm">View Plots</a></td></tr>\n')
    wlog.write('      </tbody></table></div></div></div>\n')


def plots_corrected(msinfo, wlog):
    wlog.write('<div id="calibrated" class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Calibrated visibilities</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table table-sm table-bordered" style="width:80%;">\n')
    wlog.write('      <tbody>\n')
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines(
            'plots_corrected', source,
            'Calibrated amplitude and phase against time and frequency.',
            msinfo, 'corrected')
        wlog.write(f'        <tr><td><strong>{source}</strong></td>')
        wlog.write(f'<td><a href=".{page_path}" target="_blank" class="btn btn-primary btn-sm">View Plots</a></td></tr>\n')
    wlog.write('      </tbody></table></div></div></div>\n')

def plots_uvplt(msinfo, wlog):
    all_uvplt = np.sort(glob.glob('./weblog/plots/plots_uvplt/*_uvplt_*.png'))
    calsources = [x.strip() for x in msinfo['sources']['calsources'].split(',')]
    msfilename = msinfo['msfilename']

    for p in all_uvplt:
        source_name = os.path.splitext(p)[0].split('_')[-1]
        wlog.write(f'<div id="uvplot-{source_name}" class="subsection">\n')
        wlog.write(f'  <h3 class="collapsible-header">UV Plot: {source_name}</h3>\n')
        wlog.write('  <div class="text-left">\n')
        wlog.write(f'    <a href=".{p}" target="_blank">\n')
        wlog.write(f'      <img class="img-fluid" style="max-width:100%" src=".{p}" alt="UVplot for {source_name}">\n')
        wlog.write('    </a>\n')
        wlog.write('  </div>\n')
        # Model plot for calibration sources
        if source_name in calsources:
            model_pattern = f'./weblog/plots/plots_uvplt/{msfilename}_uvpltmodel_{source_name}.png'
            if os.path.isfile(model_pattern):
                wlog.write('  <div class="text-left mt-2">\n')
                wlog.write(f'    <a href=".{model_pattern}" target="_blank">\n')
                wlog.write(f'      <img class="img-fluid" style="max-width:100%" src=".{model_pattern}" alt="Model UVplot for {source_name}">\n')
                wlog.write('    </a>\n')
                wlog.write('  </div>\n')
        wlog.write('</div>\n')



def weblog_plots(weblog, weblog_link, plots_link, msinfo):
    wlog = open(weblog_dir + "plots.html", "w")
    weblog_header(wlog, 'Plots', msinfo['run'])

    wlog.write('<div class="plots-layout">\n')
    wlog.write('<div class="plots-main">\n')
    
    # Display available plots
    if os.path.isdir('./weblog/plots/plots_data/') and os.listdir('./weblog/plots/plots_data/'):
        plots_data(msinfo, wlog)
    if os.path.isdir('./weblog/plots/plots_corrected/') and os.listdir('./weblog/plots/plots_corrected/'):
        plots_corrected(msinfo, wlog)
    if os.path.isdir('./weblog/plots/plots_uvplt/') and os.listdir('./weblog/plots/plots_uvplt/'):
        plots_uvplt(msinfo, wlog)
    
    wlog.write('</div>')  # close plots-main
    
    # --- Sidebar navigation for jump-to-section (MUST be inside plots-layout) ---
    plot_sections = []
    if os.path.isdir('./weblog/plots/plots_data/') and os.listdir('./weblog/plots/plots_data/'):
        plot_sections.append(('uncalibrated', 'Uncalibrated Visibilities'))
    if os.path.isdir('./weblog/plots/plots_corrected/') and os.listdir('./weblog/plots/plots_corrected/'):
        plot_sections.append(('calibrated', 'Calibrated Visibilities'))
    if os.path.isdir('./weblog/plots/plots_uvplt/') and os.listdir('./weblog/plots/plots_uvplt/'):
        all_uvplt = np.sort(glob.glob('./weblog/plots/plots_uvplt/*_uvplt_*.png'))
        for p in all_uvplt:
            source_name = os.path.splitext(p)[0].split('_')[-1]
            plot_sections.append((f'uvplot-{source_name}', f'UV Plot: {source_name}'))
    
    wlog.write('<nav class="plots-jump-sidebar">\n')
    wlog.write('<h4>Jump to section</h4>\n')
    wlog.write('<div class="plots-jump-links">\n')
    for section_id, section_name in plot_sections:
        wlog.write(f'<a class="plots-jump-link" href="#{section_id}">{section_name}</a>\n')
    wlog.write('</div>\n')
    wlog.write('</nav>\n')
    
    wlog.write('</div>')  # close plots-layout

    # Close the page
    weblog_foot(wlog)
    wlog.close()
