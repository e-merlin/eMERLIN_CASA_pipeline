import os
import glob
import logging
import numpy as np
from html import escape
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'


def _safe_name(text):
    """Sanitise *text* for use in a filename (must match eMCP_plots._safe_name)."""
    return ''.join(c if c.isalnum() or c in '._+-' else '_' for c in text)


def _discover_baselines(plots_dir, prefix, plot_type):
    """Find per-baseline PNGs and return sorted list of (baseline, filepath).

    Files are expected as ``{prefix}_{plot_type}-Ant1-Ant2.png`` (the
    shadems ``{_Baseline}`` placeholder produces ``-Ant1-Ant2``).
    """
    pattern = os.path.join(plots_dir, f'{prefix}_{plot_type}-*.png')
    files = sorted(glob.glob(pattern))
    results = []
    tag = f'{prefix}_{plot_type}-'
    for fpath in files:
        bname = os.path.basename(fpath)
        if bname.startswith(tag):
            baseline = bname[len(tag):].replace('.png', '')
            results.append((baseline, fpath))
    return results


def _plot_cell(wlog, fpath, alt_text):
    """Write one table cell containing an image or a 'no data' placeholder."""
    if os.path.isfile(fpath):
        wlog.write(f'      <td style="padding:2px;">'
                   f'<a href=".{fpath}" target="_blank">'
                   f'<img src=".{fpath}" alt="{escape(alt_text)}" '
                   f'style="width:100%; height:auto;">'
                   f'</a></td>\n')
    else:
        wlog.write('      <td style="padding:2px;color:#999;text-align:center;">'
                   'no data</td>\n')


def create_pnghtml_baselines(plots_path, source, subtitle, msinfo, datacolumn):
    """Create a per-source HTML page with a baseline × plot-type table.

    Each row is one baseline; columns are amp vs time, phase vs time,
    amp vs freq, phase vs freq.  Missing plots (e.g. fully flagged
    baselines) are shown as 'no data'.
    """
    page_path = weblog_dir + plots_path + '_' + source + ".html"
    plots_dir = weblog_dir + 'plots/' + plots_path
    # shadems {_field} placeholder produces a leading dash (-fieldname)
    prefix = f"{msinfo['msfilename']}-{_safe_name(source)}_{datacolumn}"

    plot_types = [
        ('amp_time',   'Amp vs Time'),
        ('phase_time', 'Phase vs Time'),
        ('amp_freq',   'Amp vs Freq'),
        ('phase_freq', 'Phase vs Freq'),
    ]

    # Collect all baselines that appear for any plot type
    all_baselines = set()
    for ptype, _ in plot_types:
        for bl, _ in _discover_baselines(plots_dir, prefix, ptype):
            all_baselines.add(bl)
    baselines_sorted = sorted(all_baselines)

    wlog = open(page_path, "w")
    weblog_header(wlog, 'Plots', msinfo['run'])
    wlog.write(f'<div class="subsection">\n')
    wlog.write(f'  <h3 class="collapsible-header">{escape(source)}</h3>\n')
    wlog.write(f'  <p>{escape(subtitle)}</p>\n')

    if not baselines_sorted:
        wlog.write('  <p>No per-baseline plots found.</p>\n')
    else:
        wlog.write('  <div class="table-responsive">\n')
        wlog.write('  <table class="table table-sm table-bordered"'
                   ' style="width:100%; table-layout:fixed;">\n')
        # Header row
        wlog.write('    <thead><tr><th style="width:70px;">Baseline</th>')
        for _, label in plot_types:
            wlog.write(f'<th style="width:calc(25% - 17.5px);">{label}</th>')
        wlog.write('</tr></thead>\n')
        # Body – one row per baseline
        wlog.write('    <tbody>\n')
        for bl in baselines_sorted:
            wlog.write(f'    <tr>\n      <td><strong>{escape(bl)}</strong></td>\n')
            for ptype, label in plot_types:
                fpath = os.path.join(plots_dir,
                                     f'{prefix}_{ptype}-{bl}.png')
                _plot_cell(wlog, fpath, f'{label} {source} {bl}')
            wlog.write('    </tr>\n')
        wlog.write('    </tbody>\n')
        wlog.write('  </table>\n')
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
        wlog.write(f'<td><a href=".{page_path}" target="_blank" '
                   f'class="btn btn-primary btn-sm">View Plots</a></td></tr>\n')
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
        wlog.write(f'<td><a href=".{page_path}" target="_blank" '
                   f'class="btn btn-primary btn-sm">View Plots</a></td></tr>\n')
    wlog.write('      </tbody></table></div></div></div>\n')


def plots_uvplt(msinfo, wlog):
    """Render UV-distance plots section.

    New naming convention from shadems:
      {msfilename}_uvplt_{field}_corrected_amp.png
      {msfilename}_uvplt_{field}_corrected_phase.png
      {msfilename}_uvplt_{field}_model_amp.png
      {msfilename}_uvplt_{field}_model_phase.png
    """
    msfilename = msinfo['msfilename']
    calsources = [x.strip() for x in
                  msinfo['sources']['calsources'].split(',')]
    mssources = [x.strip() for x in
                 msinfo['sources']['mssources'].split(',')]
    uvplt_dir = './weblog/plots/plots_uvplt/'

    for source in mssources:
        safe_src = _safe_name(source)
        prefix = f'{msfilename}_uvplt_{safe_src}'

        corr_amp   = os.path.join(uvplt_dir, f'{prefix}_corrected_amp.png')
        corr_phase = os.path.join(uvplt_dir, f'{prefix}_corrected_phase.png')
        model_amp  = os.path.join(uvplt_dir, f'{prefix}_model_amp.png')
        model_phase= os.path.join(uvplt_dir, f'{prefix}_model_phase.png')

        # Only show section if at least one corrected plot exists
        if not (os.path.isfile(corr_amp) or os.path.isfile(corr_phase)):
            continue

        wlog.write(f'<div id="uvplot-{escape(source)}" class="subsection">\n')
        wlog.write(f'  <h3 class="collapsible-header">UV Plot: {escape(source)}</h3>\n')
        wlog.write('  <div class="text-left">\n')
        wlog.write('  <table style="width:100%; table-layout:fixed; border-collapse:collapse;">\n')
        wlog.write('    <thead><tr><th style="width:50%;">Amplitude</th><th style="width:50%;">Phase</th></tr></thead>\n')

        # --- Corrected row ---
        wlog.write('  <tr>\n')
        for fpath, alt in [(corr_amp, 'Corrected Amp'),
                           (corr_phase, 'Corrected Phase')]:
            if os.path.isfile(fpath):
                wlog.write(f'    <td style="padding:4px;">'
                           f'<a href=".{fpath}" target="_blank">'
                           f'<img class="img-fluid" style="width:100%; height:auto;" '
                           f'src=".{fpath}" alt="{alt} {source}"></a></td>\n')
            else:
                wlog.write('    <td></td>\n')
        wlog.write('  </tr>\n')

        # --- Model row (calibrators only) ---
        if source in calsources and (os.path.isfile(model_amp)
                                     or os.path.isfile(model_phase)):
            wlog.write('  <tr>\n')
            for fpath, alt in [(model_amp, 'Model Amp'),
                               (model_phase, 'Model Phase')]:
                if os.path.isfile(fpath):
                    wlog.write(
                        f'    <td style="padding:4px;">'
                        f'<a href=".{fpath}" target="_blank">'
                        f'<img class="img-fluid" style="width:100%; height:auto;" '
                        f'src=".{fpath}" alt="{alt} {source}"></a></td>\n')
                else:
                    wlog.write('    <td></td>\n')
            wlog.write('  </tr>\n')

        wlog.write('  </table>\n')
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

    # --- Sidebar navigation ---
    plot_sections = []
    if os.path.isdir('./weblog/plots/plots_data/') and os.listdir('./weblog/plots/plots_data/'):
        plot_sections.append(('uncalibrated', 'Uncalibrated Visibilities'))
    if os.path.isdir('./weblog/plots/plots_corrected/') and os.listdir('./weblog/plots/plots_corrected/'):
        plot_sections.append(('calibrated', 'Calibrated Visibilities'))
    if os.path.isdir('./weblog/plots/plots_uvplt/') and os.listdir('./weblog/plots/plots_uvplt/'):
        mssources = [s.strip() for s in msinfo['sources']['mssources'].split(',')]
        for source in mssources:
            safe_src = _safe_name(source)
            prefix = f"{msinfo['msfilename']}_uvplt_{safe_src}"
            uvplt_dir = './weblog/plots/plots_uvplt/'
            if (os.path.isfile(os.path.join(uvplt_dir, f'{prefix}_corrected_amp.png'))
                    or os.path.isfile(os.path.join(uvplt_dir, f'{prefix}_corrected_phase.png'))):
                plot_sections.append((f'uvplot-{source}', f'UV Plot: {source}'))

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
