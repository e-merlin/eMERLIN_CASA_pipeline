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
    
    # Add section styling
    wlog.write('<style>\n')
    wlog.write('.download-section {\n')
    wlog.write('  background-color: #f8f9fa;\n')
    wlog.write('  border-radius: 8px;\n')
    wlog.write('  padding: 20px;\n')
    wlog.write('  margin: 20px 0;\n')
    wlog.write('  box-shadow: 0 2px 4px rgba(0,0,0,0.1);\n')
    wlog.write('}\n')
    wlog.write('.download-title {\n')
    wlog.write('  color: #007bff;\n')
    wlog.write('  margin-bottom: 15px;\n')
    wlog.write('  border-bottom: 1px solid #dee2e6;\n')
    wlog.write('  padding-bottom: 10px;\n')
    wlog.write('}\n')
    wlog.write('.download-link {\n')
    wlog.write('  display: inline-block;\n')
    wlog.write('  padding: 8px 16px;\n')
    wlog.write('  margin: 5px 0;\n')
    wlog.write('  background-color: #007bff;\n')
    wlog.write('  color: white;\n')
    wlog.write('  text-decoration: none;\n')
    wlog.write('  border-radius: 4px;\n')
    wlog.write('  transition: background-color 0.2s;\n')
    wlog.write('}\n')
    wlog.write('.download-link:hover {\n')
    wlog.write('  background-color: #0056b3;\n')
    wlog.write('}\n')
    wlog.write('.file-info {\n')
    wlog.write('  font-size: 14px;\n')
    wlog.write('  color: #6c757d;\n')
    wlog.write('  margin-top: 5px;\n')
    wlog.write('}\n')
    wlog.write('</style>\n')
    
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
