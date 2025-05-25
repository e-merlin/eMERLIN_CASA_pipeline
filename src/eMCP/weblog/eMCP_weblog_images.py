import os
import glob
import logging
import numpy as np
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'

def check_any_image(eMCP):
    msinfo = eMCP['msinfo']
    images = False
    for i in range(len(msinfo['sources']['targets'].split(','))):
        img_dir = weblog_dir + 'images/{0}/'.format(
            msinfo['sources']['targets'].split(',')[i])
        if (os.path.isdir(img_dir)) and (os.listdir(img_dir)):
            images = True
    return images

def write_img_zoom(img):
    my_string = """
    <td>
    <div class="imageBox">
      <div class="imageInn">
        <a href="{0}" target="_blank">
          <img src="{0}" alt="Image" class="img-fluid" style="max-width:600px">
        </a>
      </div>
      <div class="hoverImg">
        <img src="{1}" alt="Zoom" class="img-fluid" style="max-width:800px">
      </div>
    </div>
    </td>\n""".format(img, img.replace('.png', '_zoom.png'))
    return my_string

def show_image(eMCP, wlog, i):
    msinfo = eMCP['msinfo']
    num = 0
    ext = ''
    target = msinfo['sources']['targets'].split(',')[i]
    phscal = msinfo['sources']['phscals'].split(',')[i]
    
    # Get image statistics if available
    try:
        peak_target, noise_target, scaling_target = eMCP['img_stats'][target]
        peak_phscal, noise_phscal, scaling_phscal = eMCP['img_stats'][phscal]
    except:
        peak_target, noise_target, scaling_target = 0.0, 0.0, 0.0
        peak_phscal, noise_phscal, scaling_phscal = 0.0, 0.0, 0.0
    
    # Get image paths
    img_target = weblog_dir + 'images/{0}/{1}_{0}_img{2:02d}'.format(
        target, msinfo['msfilename'], num)
    img_phscal = weblog_dir + 'images/{0}/{1}_{0}_img{2:02d}'.format(
        phscal, msinfo['msfilename'], num)
    
    # Create image section with modern styling
    wlog.write('<div id="{0}" class="image-section mb-4">\n'.format(target))
    wlog.write('<h3 style="text-align:center">{0} <a href="#top" class="small">(up)</a></h3>\n'.format(target))
    
    # Create card with image comparison
    wlog.write('<div class="card">\n')
    wlog.write('<div class="card-body p-2">\n')
    
    # Create table with statistics and images
    wlog.write('<table class="table table-bordered" style="width:100%">\n')
    
    # Table header with statistics
    wlog.write('<thead>\n')
    wlog.write('<tr>\n')
    wlog.write('<th><strong>{0}</strong> (Target)<br>Peak: {1:3.3f} mJy (RMS: {2:3.3f} mJy/b)</th>\n'.format(
        target, peak_target*1000., noise_target*1000.))
    wlog.write('<th>Target residual</th>\n')
    wlog.write('<th><strong>{0}</strong> (Phasecal)<br>Peak: {1:5.1f} mJy (RMS: {2:5.3f} mJy/b)</th>\n'.format(
        phscal, peak_phscal*1000., noise_phscal*1000.))
    wlog.write('<th>Phasecal residual</th>\n')
    wlog.write('</tr>\n')
    wlog.write('</thead>\n')
    
    # Image row
    wlog.write('<tbody>\n')
    wlog.write('<tr>\n')
    
    # Add the images using the write_img_zoom function
    wlog.write(write_img_zoom(f'.{img_target}-image{ext}.png'))
    wlog.write(write_img_zoom(f'.{img_target}-residual{ext}.png'))
    wlog.write(write_img_zoom(f'.{img_phscal}-image{ext}.png'))
    wlog.write(write_img_zoom(f'.{img_phscal}-residual{ext}.png'))
    
    wlog.write('</tr>\n')
    wlog.write('</tbody>\n')
    wlog.write('</table>\n')
    
    wlog.write('</div>\n')
    wlog.write('</div>\n')
    wlog.write('</div>\n')
    wlog.write('<hr>\n')

def weblog_images(eMCP):
    msinfo = eMCP['msinfo']
    
    # Create the images page
    wlog = open(weblog_dir + "images.html", "w")
    weblog_header(wlog, 'Crude images', msinfo['run'])
    
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
    wlog.write('.imageBox {\n')
    wlog.write('  position: relative;\n')
    wlog.write('  width: 100%;\n')
    wlog.write('  height: 100%;\n')
    wlog.write('}\n')
    wlog.write('.imageBox .imageInn {\n')
    wlog.write('  width: 100%;\n')
    wlog.write('  height: 100%;\n')
    wlog.write('}\n')
    wlog.write('.imageBox .hoverImg {\n')
    wlog.write('  position: absolute;\n')
    wlog.write('  top: 0;\n')
    wlog.write('  left: 0;\n')
    wlog.write('  right: 0;\n')
    wlog.write('  bottom: 0;\n')
    wlog.write('  background: rgba(0, 0, 0, 0.8);\n')
    wlog.write('  display: none;\n')
    wlog.write('  justify-content: center;\n')
    wlog.write('  align-items: center;\n')
    wlog.write('  z-index: 100;\n')
    wlog.write('}\n')
    wlog.write('.imageBox:hover .hoverImg {\n')
    wlog.write('  display: flex;\n')
    wlog.write('}\n')
    wlog.write('</style>\n')
    
    # Check if any images exist
    if check_any_image(eMCP):
        # Add navigation links to targets
        wlog.write('<div class="sticky-nav">\n')
        wlog.write('  <strong>Jump to: </strong>\n')
        
        for target in msinfo['sources']['targets'].split(','):
            wlog.write('<a class="nav-link" href="#{0}">{0}</a>\n'.format(target))
        
        wlog.write('</div>\n')
        
        # Add information note
        note_text = 'Note: These are crude images of targets and phase '\
        'calibrators. The images are produced automatically with wsclean '\
        'without human intervention. These images should not be used for production.'
        
        wlog.write('<div class="alert alert-info mt-3 mb-4" style="max-width:90%; margin:0 auto;">\n')
        wlog.write('  <p>{0}</p>\n'.format(note_text))
        wlog.write('</div>\n')
        
        # Create the images for each target
        for i in range(len(msinfo['sources']['targets'].split(','))):
            img_dir = weblog_dir + 'images/{0}/'.format(
                msinfo['sources']['targets'].split(',')[i])
            if (os.path.isdir(img_dir)) and (os.listdir(img_dir)):
                show_image(eMCP, wlog, i)
    else:
        # No images found
        wlog.write('<div class="alert alert-warning" style="max-width:80%; margin:20px auto;">\n')
        wlog.write('  <h4>No images found</h4>\n')
        wlog.write('  <p>Images will appear here after they are created during the imaging steps.</p>\n')
        wlog.write('</div>\n')
    
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
