import os
import glob
import logging
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'

def weblog_download(msinfo):
    """Create a download page with links to data files"""
    # Create the download page
    wlog = open(weblog_dir + "download.html", "w")
    weblog_header(wlog, 'Download data', msinfo['run'])
    
    # Main archive section
    wlog.write('<div class="subsection centered">\n')
    wlog.write('  <h3 class="collapsible-header">Main Archive</h3>\n')
    wlog.write('  <div>\n')
    wlog.write('    <p>This tar file contains the MS and all the plots in the weblog:</p>\n')
    
    filepath = '../{}.tar'.format(msinfo['run'])
    file_size = ""
    if os.path.exists(filepath):
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
        file_size = f" ({size_mb:.1f} MB)"
    
    wlog.write('  <div class="mb-3">\n')
    wlog.write(f'    <strong>{msinfo["run"]}</strong>\n')
    wlog.write(f'    <a href="../{filepath}" target="_blank" class="download-link">Download TAR</a>\n')
    wlog.write(f'    <div class="file-info">Archive file{file_size}</div>\n')
    wlog.write('  </div>\n')
    wlog.write('  </div>\n')
    wlog.write('</div>\n')
    
    # FITS images section
    fits_files = glob.glob(weblog_dir + 'images/**/*.fits', recursive=True)
    if fits_files:
        wlog.write('<div class="subsection centered">\n')
        wlog.write('  <h3 class="collapsible-header">FITS Images</h3>\n')
        wlog.write('  <div>\n')
        wlog.write('    <p>Individual FITS image files:</p>\n')
        wlog.write('  <div class="table-responsive">\n')
        wlog.write('    <table class="table table-striped table-sm">\n')
        wlog.write('      <thead>\n')
        wlog.write('        <tr>\n')
        wlog.write('          <th>Source</th>\n')
        wlog.write('          <th>Image</th>\n')
        wlog.write('          <th>Size</th>\n')
        wlog.write('          <th>Download</th>\n')
        wlog.write('        </tr>\n')
        wlog.write('      </thead>\n')
        wlog.write('      <tbody>\n')
        
        for fits_file in sorted(fits_files):
            rel_path = fits_file.replace(weblog_dir, './')
            filename = os.path.basename(fits_file)
            source = os.path.basename(os.path.dirname(fits_file))
            
            # Get file size
            size_kb = os.path.getsize(fits_file) / 1024
            size_str = f"{size_kb:.1f} KB"
            if size_kb > 1024:
                size_mb = size_kb / 1024
                size_str = f"{size_mb:.1f} MB"
            
            wlog.write('        <tr>\n')
            wlog.write(f'          <td>{source}</td>\n')
            wlog.write(f'          <td>{filename}</td>\n')
            wlog.write(f'          <td>{size_str}</td>\n')
            wlog.write(f'          <td><a href="{rel_path}" class="btn btn-sm btn-primary" download>Download</a></td>\n')
            wlog.write('        </tr>\n')
        
        wlog.write('      </tbody>\n')
        wlog.write('    </table>\n')
        wlog.write('  </div>\n')
        wlog.write('  </div>\n')
        wlog.write('</div>\n')
    
    # Optional measurement sets section
    ms_files = glob.glob('./*.ms')
    if ms_files:
        wlog.write('<div class="subsection centered">\n')
        wlog.write('  <h3 class="collapsible-header">Measurement Sets</h3>\n')
        wlog.write('  <div>\n')
        wlog.write('    <p>Note: These files are typically large and not directly downloadable through the browser.</p>\n')
        wlog.write('    <ul class="list-group">\n')
        
        for ms_file in sorted(ms_files):
            ms_name = os.path.basename(ms_file)
            wlog.write(f'      <li class="list-group-item">{ms_name}</li>\n')
        
        wlog.write('    </ul>\n')
        wlog.write('  </div>\n')
        wlog.write('</div>\n')
    
    # Close the page
    weblog_foot(wlog)
    wlog.close()
