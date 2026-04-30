import os
import logging
from html import escape
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


def image_stats_line(peak, noise, peak_precision=3):
    if noise:
        snr = f'{peak/noise:.1f}'
    else:
        snr = '&ndash;'
    return (
        '<dl class="image-stats">'
        f'<div><dt>Peak</dt><dd>{peak*1000:.{peak_precision}f} mJy</dd></div>'
        f'<div><dt>RMS</dt><dd>{noise*1000:.3f} mJy/b</dd></div>'
        f'<div><dt>S/N</dt><dd>{snr}</dd></div>'
        '</dl>'
    )


def image_block(label, image, zoom, alt):
    return (
        '<div class="image-block">'
        '<div class="image-block-header">'
        f'<span>{escape(label)}</span>'
        f'<a class="image-zoom-link" href=".{escape(zoom, quote=True)}" target="_blank">Zoom</a>'
        '</div>'
        f'<a href=".{escape(image, quote=True)}" target="_blank"><img src=".{escape(image, quote=True)}" alt="{escape(alt, quote=True)}" class="weblog-image"></a>'
        '</div>'
    )


def image_source_header(kind, source, peak, noise, peak_precision=3):
    return (
        '<div class="image-source-header">'
        f'<div class="image-source-title">{escape(kind)}: {escape(source)}</div>'
        f'{image_stats_line(peak, noise, peak_precision=peak_precision)}'
        '</div>'
    )


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
    wlog.write(f'<div id="{escape(target)}" class="subsection">')
    wlog.write(f'<h3 class="collapsible-header">{escape(target)}</h3>')
    wlog.write('<div>')
    wlog.write('<table class="table table-bordered image-comparison-table">\n<tbody><tr>')
    # Target/Phscal images
    wlog.write('<td class="image-comparison-cell">')
    wlog.write(image_source_header('Target', target, peak_target, noise_target, peak_precision=3))
    wlog.write('<div class="image-group">')
    wlog.write(image_block('Image', img_target, img_target_zoom, 'Target image'))
    wlog.write(image_block('Residual', img_target_resid, img_target_resid_zoom, 'Target residual image'))
    wlog.write('</div>')
    wlog.write('</td>')
    wlog.write('<td class="image-comparison-cell">')
    wlog.write(image_source_header('Phase cal', phscal, peak_phscal, noise_phscal, peak_precision=1))
    wlog.write('<div class="image-group">')
    wlog.write(image_block('Image', img_phscal, img_phscal_zoom, 'Phasecal image'))
    wlog.write(image_block('Residual', img_phscal_resid, img_phscal_resid_zoom, 'Phasecal residual image'))
    wlog.write('</div>')
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
            wlog.write(f'<a class="calib-jump-link" href="#{escape(target, quote=True)}">{escape(target)}</a>\n')
    wlog.write('</div>\n')
    wlog.write('</nav>\n')
    wlog.write('</div>')  # close calib-layout
    weblog_foot(wlog)
    wlog.close()
