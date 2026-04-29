import os
import logging
from .eMCP_weblog_modern import weblog_header, weblog_foot

logger = logging.getLogger('logger')
weblog_dir = './weblog/'

def check_any_image(eMCP):
    msinfo = eMCP['msinfo']
    for target in msinfo['sources']['targets'].split(','):
        img_dir = os.path.join(weblog_dir, 'images', target)
        if os.path.isdir(img_dir) and os.listdir(img_dir):
            return True
    return False

def show_image(eMCP, wlog, i):
    msinfo = eMCP['msinfo']
    num = 0
    ext = ''
    target = msinfo['sources']['targets'].split(',')[i]
    phscal = msinfo['sources']['phscals'].split(',')[i]
    try:
        peak_target, noise_target, scaling_target = eMCP['img_stats'][target]
        peak_phscal, noise_phscal, scaling_phscal = eMCP['img_stats'][phscal]
    except:
        peak_target, noise_target, scaling_target = 0.0, 0.0, 0.0
        peak_phscal, noise_phscal, scaling_phscal = 0.0, 0.0, 0.0
    img_target = f"{weblog_dir}images/{target}/{msinfo['msfilename']}_{target}_img{num:02d}-image{ext}.png"
    img_phscal = f"{weblog_dir}images/{phscal}/{msinfo['msfilename']}_{phscal}_img{num:02d}-image{ext}.png"
    img_target_resid = f"{weblog_dir}images/{target}/{msinfo['msfilename']}_{target}_img{num:02d}-residual{ext}.png"
    img_phscal_resid = f"{weblog_dir}images/{phscal}/{msinfo['msfilename']}_{phscal}_img{num:02d}-residual{ext}.png"
    img_target_zoom = f"{weblog_dir}images/{target}/{msinfo['msfilename']}_{target}_img{num:02d}-image{ext}_zoom.png"
    img_phscal_zoom = f"{weblog_dir}images/{phscal}/{msinfo['msfilename']}_{phscal}_img{num:02d}-image{ext}_zoom.png"
    img_target_resid_zoom = f"{weblog_dir}images/{target}/{msinfo['msfilename']}_{target}_img{num:02d}-residual{ext}_zoom.png"
    img_phscal_resid_zoom = f"{weblog_dir}images/{phscal}/{msinfo['msfilename']}_{phscal}_img{num:02d}-residual{ext}_zoom.png"
    # Section (anchor) for this target
    wlog.write(f'<div id="{target}" class="subsection">')
    wlog.write(f'<h3 class="collapsible-header">{target}</h3>')
    wlog.write('<div>')
    wlog.write('<table class="table table-bordered" style="width:100%">\n<tbody><tr>')
    # Target/Phscal images
    wlog.write('<td style="vertical-align:top">')
    wlog.write(f'<div style="font-weight:bold;color:#337ab7;">Target</div>')
    wlog.write(f'<div style="font-size:90%">{target} (Peak: {peak_target*1000:.3f} mJy, RMS: {noise_target*1000:.3f} mJy/b, S/N: {peak_target/noise_target:.1f})</div>')
    wlog.write(f'<a href=".{img_target}" target="_blank"><img src=".{img_target}" alt="Target image" style="width:100%;max-width:760px;margin:4px 0;"></a>')
    wlog.write(f'<a href=".{img_target_zoom}" target="_blank">Zoom</a>')
    wlog.write('<br>')
    wlog.write(f'<div style="font-size:90%">Residual</div>')
    wlog.write(f'<a href=".{img_target_resid}" target="_blank"><img src=".{img_target_resid}" alt="Target residual image" style="width:100%;max-width:760px;margin:4px 0;"></a>')
    wlog.write(f'<a href=".{img_target_resid_zoom}" target="_blank">Zoom</a>')
    wlog.write('</td>')
    wlog.write('<td style="vertical-align:top">')
    wlog.write(f'<div style="font-weight:bold;color:#337ab7;">Phase Calibrator</div>')
    wlog.write(f'<div style="font-size:90%">{phscal} (Peak: {peak_phscal*1000:.1f} mJy, RMS: {noise_phscal*1000:.3f} mJy/b, S/N: {peak_phscal/noise_phscal:.1f})</div>')
    wlog.write(f'<a href=".{img_phscal}" target="_blank"><img src=".{img_phscal}" alt="Phasecal image" style="width:100%;max-width:760px;margin:4px 0;"></a>')
    wlog.write(f'<a href=".{img_phscal_zoom}" target="_blank">Zoom</a>')
    wlog.write('<br>')
    wlog.write(f'<div style="font-size:90%">Residual</div>')
    wlog.write(f'<a href=".{img_phscal_resid}" target="_blank"><img src=".{img_phscal_resid}" alt="Phasecal residual image" style="width:100%;max-width:760px;margin:4px 0;"></a>')
    wlog.write(f'<a href=".{img_phscal_resid_zoom}" target="_blank">Zoom</a>')
    wlog.write('</td>')
    wlog.write('</tr></tbody></table>')
    wlog.write('</div>')
    wlog.write('</div>')

def weblog_images(eMCP):
    msinfo = eMCP['msinfo']
    wlog = open(weblog_dir + "images.html", "w")
    weblog_header(wlog, 'Crude images', msinfo['run'])
    wlog.write('<div class="calib-layout">\n')
    wlog.write('<div class="calib-main">\n')
    if check_any_image(eMCP):
        targets = msinfo['sources']['targets'].split(',')
        for i, target in enumerate(targets):
            img_dir = os.path.join(weblog_dir, 'images', target)
            if os.path.isdir(img_dir) and os.listdir(img_dir):
                show_image(eMCP, wlog, i)
    else:
        wlog.write('<div class="alert alert-warning" style="max-width:80%;margin:20px auto;">')
        wlog.write('<h4>No images found</h4>')
        wlog.write('<p>Images will appear here after they are created during the imaging steps.</p>')
        wlog.write('</div>')
    wlog.write('</div>')  # close calib-main
    # Sticky Jump-to Sidebar (right)
    wlog.write('<nav class="calib-jump-sidebar">\n')
    wlog.write('<h4>Jump to section</h4>\n')
    wlog.write('<div class="calib-jump-links">\n')
    if check_any_image(eMCP):
        for target in msinfo['sources']['targets'].split(','):
            wlog.write(f'<a class="calib-jump-link" href="#{target}">{target}</a>\n')
    wlog.write('</div>\n')
    wlog.write('</nav>\n')
    wlog.write('</div>')  # close calib-layout
    weblog_foot(wlog)
    wlog.close()
