import os
from ..utils import eMCP_paths as empaths
from .eMCP_weblog_modern import weblog_header, weblog_foot, write_link_txt

def weblog_index(msinfo):
    """Create the Home (index.html) page for the weblog."""
    wlog = open(os.path.join(empaths.WEBLOG_DIR, "index.html"), "w")
    weblog_header(wlog, 'Home', msinfo.get('run', msinfo.get('project', 'eMERLIN')))

    wlog.write('<div class="section centered">\n')
    wlog.write('  <h2 class="section-title">Project Information</h2>\n')

    # Info table
    wlog.write('  <table class="table" style="width:50%; margin: 0 auto;">\n')
    wlog.write(f'    <tr><td>Project</td><td>{msinfo.get("project", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>Run</td><td>{msinfo.get("run", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>MS file</td><td>{msinfo.get("msfile", "Unknown")}</td></tr>\n')
    if 't_ini' in msinfo and 't_end' in msinfo:
        wlog.write(f'    <tr><td>Start</td><td>{msinfo["t_ini"].strftime("%Y-%m-%d %H:%M")}</td></tr>\n')
        wlog.write(f'    <tr><td>End</td><td>{msinfo["t_end"].strftime("%Y-%m-%d %H:%M")}</td></tr>\n')
    wlog.write(f'    <tr><td>Band</td><td>{msinfo.get("band", "Unknown")}</td></tr>\n')
    if 'antennas' in msinfo:
        wlog.write(f'    <tr><td>Antennas</td><td>{", ".join(msinfo["antennas"])}</td></tr>\n')
    if 'sources' in msinfo and 'mssources' in msinfo['sources']:
        num_sources = len(msinfo['sources']['mssources'].split(','))
        wlog.write(f'    <tr><td>Number of sources</td><td>{num_sources}</td></tr>\n')
    wlog.write(f'    <tr><td>Integration time</td><td>{msinfo.get("int_time", "Unknown")}s</td></tr>\n')
    if 'freq_ini' in msinfo and 'freq_end' in msinfo:
        wlog.write(f'    <tr><td>Frequency</td><td>{msinfo["freq_ini"]:5.2f} - {msinfo["freq_end"]:5.2f} GHz</td></tr>\n')
    wlog.write(f'    <tr><td>Num. spw</td><td>{msinfo.get("num_spw", "Unknown")}</td></tr>\n')
    wlog.write(f'    <tr><td>Channels/spw</td><td>{msinfo.get("nchan", "Unknown")}</td></tr>\n')
    if 'chan_res' in msinfo:
        chan_width = msinfo["chan_res"] * 1000  # MHz
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
    notes_file = f'./{msinfo.get("run", "")}.notes.txt'
    if os.path.isfile(notes_file):
        wlog.write('  <div class="subsection centered">\n')
        wlog.write('    <h3>Notes and Comments</h3>\n')
        write_link_txt(wlog, notes_file, 'Notes and observing comments')
        wlog.write('  </div>\n')
    wlog.write('</div>\n')

    weblog_foot(wlog)
    wlog.close()
