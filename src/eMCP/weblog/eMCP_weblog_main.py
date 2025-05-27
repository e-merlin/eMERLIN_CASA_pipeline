import os
import shutil
import logging
import glob
from ..utils import eMCP_utils as emutils

logger = logging.getLogger('logger')

# Import all modernized weblog components
from .eMCP_weblog_modern import weblog_index, weblog_obssum
from .eMCP_weblog_plots import weblog_plots
from .eMCP_weblog_calib import weblog_calibration
from .eMCP_weblog_images import weblog_images
from .eMCP_weblog_download import weblog_download
from .eMCP_weblog_flagstats import weblog_flagstats

def makedir(pathdir):
    """Make directory if it doesn't exist"""
    if not os.path.exists(pathdir):
        os.makedirs(pathdir)





def weblog_pipelineinfo(eMCP):
    """Create pipeline info page similar to the original"""
    weblog_dir = './weblog/'
    info_dir = './weblog/info/'
    info_link = './info/'
    calib_dir = './weblog/calib/'
    
    # Import weblog functions from the modern module
    from .eMCP_weblog_modern import weblog_header, weblog_foot, write_link_txt
    
    wlog = open(weblog_dir + "pipelineinfo.html", "w")
    # Get run name if available, otherwise use a default
    run_name = eMCP.get('run', eMCP.get('project', 'Pipeline'))
    weblog_header(wlog, 'Pipeline info', run_name)
    
    # Pipeline version
    pipeline_version = eMCP.get('pipeline_version', 'Unknown')
    wlog.write(f'<p>Pipeline version: {pipeline_version}</p>\n')
    
    # Generate pipeline steps table in original format
    # We'll create a function similar to table_steps from the original
    steps_table = generate_steps_table(eMCP)
    wlog.write(steps_table)
    
    # Color legend
    wlog.write('<p>Green = executed<br>\n')
    wlog.write('Red = executed but outdated by a previous step</p>\n')
    
    # Relevant log files section
    wlog.write('<h4>Relevant log files:</h4>\n')
    
    # Pipeline log
    if os.path.isfile(info_dir + 'eMCP.log.txt'):
        write_link_txt(wlog, info_link + 'eMCP.log.txt', 'Pipeline log', text='eMCP.log')
    
    # CASA log
    if os.path.isfile(info_dir + 'casa_eMCP.log.txt'):
        write_link_txt(wlog, info_link + 'casa_eMCP.log.txt', 'CASA log', text='casa_eMCP.log')
    
    # Relevant parameter files section
    wlog.write('<h4>Relevant parameter files:</h4>\n')
    
    # Create eMCP_info.txt if it doesn't exist
    if isinstance(eMCP, dict) and not os.path.isfile(info_dir + 'eMCP_info.txt'):
        try:
            from ..utils import eMCP_utils as emutils
            emutils.prt_dict_tofile(eMCP, tofilename=info_dir + 'eMCP_info.txt', pre='  ')
        except Exception as e:
            logger.warning(f'Error creating eMCP_info.txt: {e}')
    
    # Link to eMCP_info.txt
    if os.path.isfile(info_dir + 'eMCP_info.txt'):
        write_link_txt(wlog, info_link + 'eMCP_info.txt', 'Pipeline info (dict)', text='eMCP_info.txt')
    
    # Create and link to caltables.txt if available
    if os.path.isfile(calib_dir + 'caltables.yaml'):
        try:
            from pickle import load
            from ..utils import eMCP_utils as emutils
            
            with open(calib_dir + 'caltables.yaml', 'r') as f:
                caltables = load(f)
            
            emutils.prt_dict_tofile(caltables, tofilename=info_dir + 'caltables.txt', pre='  ')
            write_link_txt(wlog, info_link + 'caltables.txt', 'Calibration info (dict)', text='caltables.txt')
        except Exception as e:
            logger.warning(f'Error creating caltables.txt: {e}')
    
    weblog_foot(wlog)
    wlog.close()


def generate_steps_table(eMCP):
    """Generate HTML table for pipeline steps using the original color coding logic"""
    import datetime
    
    if 'input_steps' not in eMCP:
        return '<p>No pipeline step information available</p>'
    
    # Get step information
    steps_info = eMCP.get('steps', {})
    
    html = '<table class="table" style="width:80%; margin: 0 auto;">\n'
    html += '<tr>'
    html += '<th style="width: 15%;">Step name</th>'
    html += '<th style="width: 4%;">Code</th>'
    html += '<th style="width: 25%;">Execution ended</th>'
    html += '<th style="width: 15%;">Execution time</th>'
    html += '<th style="width: 40%;">Notes</th>'
    html += '</tr>\n'
    
    # Used to track executed steps for determining outdated steps
    executed_steps = []
    # Skip these steps when determining outdated status (as in original code)
    nosteps = ['plot_data', 'save_flags', 'plot_corrected', 'first_images']
    
    # Colors from the original
    neutral_color = '#EEEEEE'  # Light gray - not executed
    green_color = '#008B45'    # Green - executed and current
    red_color = '#FF6347'      # Red - executed but outdated
    
    # Build the table
    for step in eMCP.get('input_steps', {}).keys():
        row = '<tr>'
        
        # Step name
        row += f'<td>{step}</td>'
        
        # Get step information
        step_info = steps_info.get(step, [0, 0, ''])
        time_s, delta_min, msg = step_info
        
        # Code/Run column with color coding (using the original logic)
        color = neutral_color
        status = ''
        status_code = ''
        
        if step in eMCP.get('input_steps', {}):
            status_code = str(int(eMCP['input_steps'][step]))
        
        # Determine color based on execution status
        if time_s == 0:
            # Not executed
            color = neutral_color
            status = ''
        elif isinstance(time_s, str):
            # Step has been executed, check if it's outdated
            time_fmt = '%Y-%m-%d %H:%M:%S'
            try:
                # Convert to datetime
                time_step = datetime.datetime.strptime(time_s.split('.')[0], '%Y-%m-%d %H:%M:%S')
                
                if executed_steps:
                    # Find the latest previous step execution time
                    latest_times = []
                    for prev_step in executed_steps:
                        prev_time = steps_info.get(prev_step, [0, 0, ''])[0]
                        if isinstance(prev_time, str):
                            try:
                                latest_times.append(datetime.datetime.strptime(prev_time.split('.')[0], '%Y-%m-%d %H:%M:%S'))
                            except:
                                pass
                    
                    if latest_times:
                        latest_prev = max(latest_times)
                        if time_step < latest_prev:
                            # Outdated - executed before a later step
                            color = red_color
                            status = ''
                        else:
                            # Current - executed after all previous steps
                            color = green_color
                            status = 'OK'
                    else:
                        color = green_color
                        status = 'OK'
                else:
                    # First executed step
                    color = green_color
                    status = 'OK'
            except Exception as e:
                # Handle formatting errors
                color = neutral_color
                status = ''
        
        # Add the code/status cell
        row += f'<td align="center" bgcolor="{color}">{status_code}</td>'
        
        # Execution ended timestamp
        if time_s and time_s != 0:
            row += f'<td align="center">{time_s}</td>'
        else:
            row += '<td align="center">-</td>'
        
        # Execution time (formatted)
        if isinstance(delta_min, (int, float)) and delta_min > 0:
            # Format execution time as in original
            if delta_min < 1:
                formatted_time = '<1 min'
            elif delta_min >= 1 and delta_min < 50:
                formatted_time = '{:.0f} min'.format(delta_min)
            elif delta_min >= 50:
                formatted_time = '{:.1f} h'.format(delta_min / 60.0)
            else:
                formatted_time = '-'
            
            row += f'<td align="center">{formatted_time}</td>'
        else:
            row += '<td align="center">-</td>'
        
        # Notes
        row += f'<td>{msg}</td>'
        
        row += '</tr>\n'
        html += row
        
        # Track executed steps for determining outdated status
        if time_s != 0 and step not in nosteps:
            executed_steps.append(step)
    
    html += '</table>\n'
    return html

def start_weblog(eMCP, silent=False):
    """Main function to create modernized weblog with better error handling"""
    # Early exit if eMCP is not a valid dictionary
    if not isinstance(eMCP, dict):
        logger.warning('Invalid eMCP dictionary provided, weblog not created')
        return
    """Main function to create modernized weblog"""
    # Setup directories
    weblog_dir = './weblog/'
    info_dir = './weblog/info/'
    calib_dir = './weblog/calib/'
    plots_dir = './weblog/plots/'
    logs_dir = './logs/'
    images_dir = './weblog/images/'
    flagstats_dir = './weblog/flagstats/'
    
    # Create directories if they don't exist
    for directory in [weblog_dir, info_dir, calib_dir, plots_dir, images_dir, flagstats_dir]:
        makedir(directory)
    
    # Copy static files (CSS, JS, images)
    utils_path = os.path.dirname(os.path.abspath(emutils.__file__))
    
    # Copy legacy files for compatibility
    os.system('cp -p {0}/emerlin-2.gif {1}'.format(utils_path, weblog_dir))
    os.system('cp -p {0}/eMCP.css {1}'.format(utils_path, weblog_dir))
    os.system('cp -p {0}/eMCP_logo.png {1}'.format(utils_path, weblog_dir))
    
    # Copy modern files
    os.system('cp -p {0}/eMCP_modern.css {1}'.format(utils_path, weblog_dir))
    os.system('cp -p {0}/eMCP_modern.js {1}'.format(utils_path, weblog_dir))
    
    # Always create the pipeline info page which doesn't require msinfo
    weblog_pipelineinfo(eMCP)
        
    # Generate the rest of the weblog pages if we have MS info
    if 'msinfo' in eMCP and isinstance(eMCP['msinfo'], dict):
        msinfo = eMCP['msinfo']
        
        try:
            # Create index page
            weblog_index(msinfo)
            
            # Create observation summary page
            weblog_obssum(msinfo)
            
            # Create calibration page if required info is available
            weblog_calibration(eMCP)
            
            # Create plots page
            weblog_plots(weblog_dir, './', './plots/', msinfo)
            
            # Create flag statistics page
            weblog_flagstats(msinfo)
            
            # Create images page
            weblog_images(eMCP)
            
            # Create download page
            weblog_download(msinfo)
            
            if not silent:
                logger.info('Created modernized weblog in ./weblog/')
        except Exception as e:
            logger.warning(f'Error creating some weblog pages: {e}')
            logger.info('Basic weblog created with available information')
    else:
        logger.warning('No MS info available, only pipeline info page created.')

if __name__ == "__main__":
    # For testing or direct execution
    import sys
    if len(sys.argv) > 1:
        from pickle import load
        with open(sys.argv[1], 'rb') as f:
            eMCP = load(f)
        start_weblog(eMCP)
    else:
        print("Usage: python -m eMCP.weblog.eMCP_weblog_main <eMCP_pickle_file>")
