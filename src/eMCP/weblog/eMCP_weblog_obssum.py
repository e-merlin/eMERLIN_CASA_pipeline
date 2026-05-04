import os
import glob
from ..utils import eMCP_paths as empaths
from .eMCP_weblog_modern import weblog_header, weblog_foot, write_link_txt

def weblog_obssum(msinfo):
    """Creates the Observation summary page (obs_summary.html) for the weblog."""
    weblog_dir = empaths.WEBLOG_DIR
    info_dir = empaths.INFO_DIR
    info_link = empaths.INFO_LINK

    wlog = open(os.path.join(weblog_dir, "obs_summary.html"), "w")
    weblog_header(wlog, 'Observation summary', msinfo.get('run', msinfo.get('project', 'eMERLIN')))

    # Summary section
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Summary</h3>\n')
    wlog.write('  <div>\n')

    # Listobs link
    listobs_file = os.path.join(info_link, msinfo.get('msfile', 'UNKNOWN') + '.listobs.txt')
    write_link_txt(wlog, listobs_file, 'Summary of current observation (listobs)')

    # Other listobs files if available
    all_listobs = [
        os.path.basename(l)
        for l in glob.glob(os.path.join(info_dir, msinfo.get('run', '') + '*listobs.txt'))
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
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Sources</h3>\n')
    wlog.write('  <div>\n')

    # Source separation table
    sepfile = info_link + 'source_separations.txt'
    sources = msinfo.get('sources', {})
    targets = sources.get('targets', '')
    phscals = sources.get('phscals', '')
    mssources = sources.get('mssources', '')
    separations = msinfo.get('separations', {})

    if targets != '':
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

        target_list = targets.split(',') if isinstance(targets, str) else []
        phscal_list = phscals.split(',') if isinstance(phscals, str) else []

        for s1, s2 in zip(target_list, phscal_list):
            wlog.write('      <tr>\n')
            if s1 not in mssources:
                separation = 'Not in MS'
            elif s2 not in mssources:
                separation = 'Not in MS'
            else:
                try:
                    separation = '{0:5.2f}'.format(separations[s1 + '-' + s2])
                except Exception:
                    try:
                        separation = '{0:5.2f}'.format(separations[s2 + '-' + s1])
                    except Exception:
                        separation = '0.0'
            wlog.write(f'        <td>{s1}</td>\n')
            wlog.write(f'        <td>{s2}</td>\n')
            wlog.write(f'        <td>{separation}</td>\n')
            wlog.write('      </tr>\n')

        wlog.write('    </tbody>\n')
        wlog.write('  </table>\n')
        wlog.write('</div>\n')

        # Add link to the separations file
        if os.path.isfile(os.path.join(info_dir, 'source_separations.txt')):
            write_link_txt(wlog, sepfile, 'View all source separations')
    else:
        wlog.write('<p>No target sources found</p>\n')

    wlog.write('  </div>\n')
    wlog.write('</div>\n')

    # Sources in MS section
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Sources in MS</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <div class="table-responsive">\n')
    wlog.write('      <table class="table">\n')
    wlog.write('        <thead>\n')
    wlog.write('          <tr>\n')
    wlog.write('            <th>Source</th>\n')
    wlog.write('            <th>Intent</th>\n')
    wlog.write('            <th>Coordinates (h:m:s d:m:s)</th>\n')
    wlog.write('            <th>MJD range</th>\n')
    wlog.write('          </tr>\n')
    wlog.write('        </thead>\n')
    wlog.write('        <tbody>\n')

    if 'sources' in msinfo and 'mssources' in msinfo['sources']:
        for source in msinfo['sources']['mssources'].split(','):
            # MJD range
            if 'source_timerange_mjd' in msinfo['sources'] and source in msinfo['sources']['source_timerange_mjd']:
                mjd_ini, mjd_end = msinfo['sources']['source_timerange_mjd'][source]
                mjd_range = f"{mjd_ini:11.5f} - {mjd_end:11.5f}"
            else:
                mjd_range = "Unknown"
            # Intent
            intent = msinfo['sources'].get('source_intent', {}).get(source, 'Unknown')
            # Coordinates
            coords = msinfo.get('directions', {}).get(source, "Unknown")
            wlog.write('          <tr>\n')
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
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Antennas</h3>\n')
    wlog.write('  <div>\n')
    refant = msinfo.get('refant') or 'not defined yet'
    first_refant = refant.split(',')[0].strip()
    if 'antennas' in msinfo:
        wlog.write('    <div class="antenna-list" aria-label="Antennas">\n')
        for ant in msinfo['antennas']:
            ant_class = 'antenna-code antenna-code-ref' if ant == first_refant else 'antenna-code'
            wlog.write(f'      <span class="{ant_class}">{ant}</span>\n')
        wlog.write('    </div>\n')
    wlog.write('    <p>Reference antenna: {}</p>\n'.format(refant))
    wlog.write('  </div>\n')
    wlog.write('</div>\n')

    # Source elevation section
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3 class="collapsible-header">Source elevation</h3>\n')
    wlog.write('  <div>\n')

    # Elevation plots (try msfilename, fallback to any elevation plot)
    msfilename = msinfo.get('msfilename', None)
    plots_dir = empaths.PLOTS_OBSERVATION_DIR
    elev_plots = []
    if msfilename:
        elev_plots = glob.glob(f'{plots_dir}{msfilename}_elevation.png')
    if not elev_plots:
        elev_plots = glob.glob(f'{plots_dir}*elevation*.png')
    if elev_plots:
        plot_path = elev_plots[0]
        rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
        wlog.write('    <div >\n')
        wlog.write(f'      <img src="{rel_path}" style="width:100%; max-width:850px;" alt="Source elevation plot">\n')
        wlog.write('    </div>\n')
    else:
        wlog.write('    <p>No elevation plot available</p>\n')
    wlog.write('  </div>\n')
    wlog.write('</div>\n')

    # UV coverage section
    wlog.write('<div class="subsection">\n')
    wlog.write('  <h3>UV coverage</h3>\n')

    uvcov_plots = []
    if msfilename:
        uvcov_plots = sorted(glob.glob(f'{plots_dir}{msfilename}_uvcov*.png'))
    if not uvcov_plots:
        uvcov_plots = sorted(glob.glob(f'{plots_dir}*uv*.png') + glob.glob(f'{plots_dir}*UV*.png'))
    if uvcov_plots:
        wlog.write('  <div >\n')
        for plot_path in uvcov_plots:
            rel_path = os.path.join('plots', 'plots_observation', os.path.basename(plot_path))
            wlog.write(f'    <img src="{rel_path}" style="width:30%; max-width:600px; margin-bottom:20px;" alt="UV coverage plot">\n')
        wlog.write('  </div>\n')
    else:
        wlog.write('  <p>No UV coverage plots available</p>\n')
    wlog.write('</div>\n')

    weblog_foot(wlog)
    wlog.close()
