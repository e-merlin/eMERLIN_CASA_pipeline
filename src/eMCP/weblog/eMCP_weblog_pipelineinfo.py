import os

def weblog_pipelineinfo(eMCP):
    """Create pipeline info page (pipelineinfo.html) for the weblog."""
    from .eMCP_weblog_modern import weblog_header, weblog_foot, write_link_txt

    weblog_dir = './weblog/'
    info_dir = './weblog/info/'
    info_link = './info/'
    calib_dir = './weblog/calib/'

    wlog = open(os.path.join(weblog_dir, "pipelineinfo.html"), "w")
    run_name = eMCP['msinfo'].get('run', eMCP['msinfo'].get('project', 'Pipeline'))
    weblog_header(wlog, 'Pipeline info', run_name)

    # Pipeline version
    pipeline_version = eMCP.get('pipeline_version', 'Unknown')
    wlog.write(f'<p>Pipeline version: {pipeline_version}</p>\n')

    # Steps table
    steps_table = generate_steps_table(eMCP)
    wlog.write(steps_table)

    wlog.write('<p>Green = executed<br>\n')
    wlog.write('Red = executed but outdated by a previous step</p>\n')

    # Relevant log files section
    wlog.write('<h4>Relevant log files:</h4>\n')
    if os.path.isfile(info_dir + 'eMCP.log.txt'):
        write_link_txt(wlog, info_link + 'eMCP.log.txt', 'Pipeline log', text='eMCP.log')
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
            pass  # Silent fail

    if os.path.isfile(info_dir + 'eMCP_info.txt'):
        write_link_txt(wlog, info_link + 'eMCP_info.txt', 'Pipeline info (dict)', text='eMCP_info.txt')

    # Create and link to caltables.txt if available
    if os.path.isfile(calib_dir + 'caltables.yaml'):
        try:
            import yaml
            from ..utils import eMCP_utils as emutils
            with open(calib_dir + 'caltables.yaml', 'r') as f:
                caltables = yaml.safe_load(f)
            emutils.prt_dict_tofile(caltables, tofilename=info_dir + 'caltables.txt', pre='  ')
            write_link_txt(wlog, info_link + 'caltables.txt', 'Calibration info (dict)', text='caltables.txt')
        except Exception as e:
            pass  # Silent fail

    weblog_foot(wlog)
    wlog.close()

def generate_steps_table(eMCP):
    """Generate HTML table for pipeline steps with color coding."""
    import datetime

    if 'input_steps' not in eMCP:
        return '<p>No pipeline step information available</p>'

    steps_info = eMCP.get('steps', {})

    html = '<table class="table" style="width:80%; margin: 0 auto;">\n'
    html += '<tr>'
    html += '<th style="width: 15%;">Step name</th>'
    html += '<th style="width: 4%;">Code</th>'
    html += '<th style="width: 25%;">Execution ended</th>'
    html += '<th style="width: 15%;">Execution time</th>'
    html += '<th style="width: 40%;">Notes</th>'
    html += '</tr>\n'

    executed_steps = []
    nosteps = ['plot_data', 'save_flags', 'plot_corrected', 'first_images']
    neutral_color = '#EEEEEE'
    green_color = '#008B45'
    red_color = '#FF6347'

    for step in eMCP.get('input_steps', {}).keys():
        row = '<tr>'
        row += f'<td>{step}</td>'
        step_info = steps_info.get(step, [0, 0, ''])
        time_s, delta_min, msg = step_info
        color = neutral_color
        status_code = ''
        if step in eMCP.get('input_steps', {}):
            status_code = str(int(eMCP['input_steps'][step]))

        if time_s == 0:
            color = neutral_color
        elif isinstance(time_s, str):
            time_fmt = '%Y-%m-%d %H:%M:%S'
            try:
                time_step = datetime.datetime.strptime(time_s.split('.')[0], '%Y-%m-%d %H:%M:%S')
                if executed_steps:
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
                            color = red_color
                        else:
                            color = green_color
                    else:
                        color = green_color
                else:
                    color = green_color
            except Exception:
                color = neutral_color

        row += f'<td align="center" bgcolor="{color}">{status_code}</td>'

        if time_s and time_s != 0:
            row += f'<td align="center">{time_s}</td>'
        else:
            row += '<td align="center">-</td>'

        if isinstance(delta_min, (int, float)) and delta_min > 0:
            if delta_min < 1:
                formatted_time = '<1 min'
            elif delta_min < 50:
                formatted_time = '{:.0f} min'.format(delta_min)
            elif delta_min >= 50:
                formatted_time = '{:.1f} h'.format(delta_min / 60.0)
            else:
                formatted_time = '-'
            row += f'<td align="center">{formatted_time}</td>'
        else:
            row += '<td align="center">-</td>'

        row += f'<td>{msg}</td>'
        row += '</tr>\n'
        html += row

        if time_s != 0 and step not in nosteps:
            executed_steps.append(step)

    html += '</table>\n'
    return html

