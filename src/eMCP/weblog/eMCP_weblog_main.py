import os
import logging
from ..utils import eMCP_utils as emutils
from ..utils import eMCP_paths as empaths

logger = logging.getLogger('logger')

# Import each page's content generator (all content logic in their own files)
from .eMCP_weblog_index import weblog_index
from .eMCP_weblog_obssum import weblog_obssum
from .eMCP_weblog_plots import weblog_plots
from .eMCP_weblog_calib import weblog_calibration
from .eMCP_weblog_images import weblog_images
from .eMCP_weblog_download import weblog_download
from .eMCP_weblog_flagstats import weblog_flagstats
from .eMCP_weblog_pipelineinfo import weblog_pipelineinfo


def makedir(pathdir):
    """Make directory if it doesn't exist"""
    if not os.path.exists(pathdir):
        os.makedirs(pathdir)


def start_weblog(eMCP, silent=False):
    """Main function to create weblog, only delegates to subpage scripts."""
    if not isinstance(eMCP, dict):
        logger.warning('Invalid eMCP dictionary provided, weblog not created')
        return

    # Setup directories
    weblog_dir = empaths.WEBLOG_DIR
    info_dir = empaths.INFO_DIR
    calib_dir = empaths.CALIB_DIR
    plots_dir = empaths.PLOTS_DIR
    logs_dir = empaths.LOGS_DIR
    images_dir = empaths.IMAGES_DIR
    flagstats_dir = empaths.FLAGSTATS_DIR

    # Create directories if they don't exist
    for directory in [weblog_dir, info_dir, calib_dir, plots_dir, images_dir, flagstats_dir]:
        makedir(directory)

    # Copy static files (CSS, images).
    utils_path = os.path.dirname(os.path.abspath(emutils.__file__))
    os.system(f'cp -p {utils_path}/emerlin-2.gif {weblog_dir}')
    os.system(f'cp -p {utils_path}/eMCP_logo.png {weblog_dir}')
    os.system(f'cp -p {utils_path}/eMCP_modern.css {weblog_dir}')
    # (No JS copied)

    # Always create the pipeline info page
    weblog_pipelineinfo(eMCP)

    # Generate the rest of the weblog pages if MS info is available
    if 'msinfo' in eMCP and isinstance(eMCP['msinfo'], dict):
        msinfo = eMCP['msinfo']
        pages = [
            ('index', weblog_index, (msinfo,)),
            ('observation summary', weblog_obssum, (msinfo,)),
            ('calibration', weblog_calibration, (eMCP,)),
            ('plots', weblog_plots, (weblog_dir, './', './plots/', msinfo)),
            ('flag statistics', weblog_flagstats, (msinfo,)),
            ('images', weblog_images, (eMCP,)),
            ('download', weblog_download, (msinfo,)),
        ]
        for page_name, page_func, page_args in pages:
            try:
                page_func(*page_args)
            except Exception as e:
                logger.warning('Error creating weblog %s page: %s',
                               page_name, e)
        if not silent:
            logger.info('Created weblog in %s', empaths.WEBLOG_DIR)
    else:
        logger.warning('No MS info available, only pipeline info page created.')


if __name__ == "__main__":
    # For testing or direct execution
    import sys
    if len(sys.argv) > 1:
        import yaml
        with open(sys.argv[1], 'r') as f:
            eMCP = yaml.safe_load(f)
        start_weblog(eMCP)
    else:
        print("Usage: python -m eMCP.weblog.eMCP_weblog_main <eMCP_yaml_file>")
