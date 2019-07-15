import os
import numpy as np
import glob
import pickle
import collections
import datetime

# CASA imports
from taskinit import *
from tasks import *
import casadef

import logging
logger = logging.getLogger('logger')

weblog_dir = './weblog/'
info_dir   = './weblog/info/'
calib_dir  = './weblog/calib/'
plots_dir  = './weblog/plots/'
logs_dir   = './logs/'
images_dir = './weblog/images/'
flagstats_dir = './weblog/flagstats/'

weblog_link= './'
info_link  = './info/'
calib_link = './calib/'
plots_link = './plots/'
images_link= './images/'

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def prt_dict_tofile(d, tofilename=None, addfile='', pre=' '):
    if tofilename != None:
        f = open(tofilename, 'wb')
    else:
        f = addfile
    subdict = []
    for key in d.keys():
        if type(d[key]) in [collections.OrderedDict, dict]:
            subdict.append(key)
        else:
            f.write('{0:20s}: {1}\n'.format(pre+key, d[key]))
    if subdict != []:
        for key_inner in subdict:
            f.write('{}\n'.format(pre+key_inner))
            prt_dict_tofile(d[key_inner], addfile=f, pre=pre+pre)


def weblog_button(weblog_link, name, link):
    line = ('<input type="button" value="{0}" '
            'onclick="window.location.href=\''
            '{1}{2}.html\'">\n').format(name, weblog_link, link)
    return line

def weblog_header(wlog, section, project):
    wlog.write('<!DOCTYPE html>\n')
    wlog.write('<html lang="eng">\n')
    wlog.write('<link rel="stylesheet" type="text/css" href="{0}eMCP.css"/>\n'.format(weblog_link))
    wlog.write('<head>\n')
    wlog.write('<title>{} - eMCP</title>\n'.format(project))
    wlog.write('</head>\n')
    wlog.write('<body>\n')
    wlog.write('<h4 id="top">e-MERLIN Pipeline Web Log<br>{}</h4>\n'.format(project))
    wlog.write('<form>\n')
    wlog.write(weblog_button(weblog_link, 'Home', 'index'))
    wlog.write(weblog_button(weblog_link, 'Observation summary',
                             'obs_summary'))
    wlog.write(weblog_button(weblog_link, 'Pipeline info', 'pipelineinfo'))
    wlog.write(weblog_button(weblog_link, 'Calibration', 'calibration'))
    wlog.write(weblog_button(weblog_link, 'Plots', 'plots'))
    wlog.write(weblog_button(weblog_link, 'Flag statistics', 'flagstats'))
    wlog.write(weblog_button(weblog_link, 'Images', 'images'))
    wlog.write(weblog_button(weblog_link, 'Download data', 'download'))
    wlog.write('</form>\n')
    #wlog.write('<br>\n')
    wlog.write('<h2 id="{0}">{0}</h2>\n'.format(section))
    wlog.write('<hr>\n')


def weblog_foot(wlog):
    wlog.write('<br><br>\n<hr>\n')
    wlog.write('<a href="http://www.e-merlin.ac.uk/" target="_blank"><img ' \
               'src="{0}emerlin-2.gif"></a><br>\n'.format(weblog_link))
    wlog.write('<small>')
    links = collections.OrderedDict()
    links['e-MERLIN Pipeline Github'] = 'https://github.com/e-merlin/eMERLIN_CASA_pipeline'
    links['e-MERLIN']  = 'http://www.e-merlin.ac.uk/'
    links['User support and data reduction for e-MERLIN'] = 'http://www.e-merlin.ac.uk/data_red/'
    links['CASA'] = 'https://casa.nrao.edu/'
    for n, l in links.items():
        wlog.write('- <a href="{0}" target="_blank">{1}</a><br>\n'.format(l, n))
    #wlog.write('<hr>\n')
    wlog.write('</small>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext, text='txt'):
    wlog.write('{0}: <a href="{1}" target="_blank">{2}</a><br>\n'.format(intext, infile, text))

def weblog_index(msinfo):
    wlog = open("./weblog/index.html","w")
    weblog_header(wlog, 'Home', msinfo['run'])
    #------------------------------------------
    wlog.write('<table class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:40%">\n')
    wlog.write('<tr><td>Project </td>   <td> {}</td>\n'.format(msinfo['project']))
    wlog.write('<tr><td>Run </td>   <td> {}</td>\n'.format(msinfo['run']))
    wlog.write('<tr><td>MS file </td>   <td> {}</td>\n'.format(msinfo['msfile']))
    #wlog.write('<tr><td>Date observed </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].date()))
    wlog.write('<tr><td>Start </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>End </td>    <td> {0} </td>\n'.format(msinfo['t_end'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>Band </td>   <td> {0} </td>\n'.format(msinfo['band']))
    wlog.write('<tr><td>Antennas </td>   <td> {0} </td>\n'.format(', '.join(msinfo['antennas'])))
    wlog.write('<tr><td>Number of sources </td>   <td> {0} </td>\n'.format(len(msinfo['sources']['mssources'].split(','))))
    wlog.write('<tr><td>Integration time </td>             <td> {0}s</td>\n'.format(msinfo['int_time']))
    wlog.write('<tr><td>Frequency </td>  <td> {0:5.2f} - {1:5.2f} GHz </td>\n'.format(
        msinfo['freq_ini'], msinfo['freq_end']))
    wlog.write('<tr><td>Num. spw </td>             <td> {0} </td>\n'.format(msinfo['num_spw']))
    wlog.write('<tr><td><abbr title="Number of channels per spw after  averaging. The raw data may provide up to 512 channels">Channels/spw<sup>*</sup></abbr></td><td> {0} </td>\n'.format(msinfo['nchan']))
    wlog.write('<tr><td>Channel width </td><td> {0:6.2f} MHz</td>\n'.format(msinfo['chan_res']*1000))
    wlog.write('<tr><td>spw bandwidth </td>   <td> {0:6.0f} MHz</td>\n'.format(
    msinfo['chan_res']*1000*msinfo['nchan']))
    wlog.write('<tr><td>Total bandwidth </td><td> {0:6.0f} MHz</td>\n'.format(
    msinfo['chan_res']*1000*msinfo['nchan']*msinfo['num_spw']))
    #wlog.write('<tr><td><abbr title="Integration time after averaging. The raw data may provide higher time resolution">Integration time<sup>*</sup></abbr></td><td> {0:3.1f} s</td>\n'.format(int_time))
    wlog.write('<tr><td>Polarizations </td>  <td> {0}  </td>\n'.format(msinfo['polarizations']))
    #wlog.write('<tr><td>Stokes </td>         <td> '+', '.join(str(x) for x in stokes)+'</td>\n')
    #wlog.write('<tr><td>Observer </td>       <td> '+uvdata.header['observer']+'</td>\n')
    #wlog.write('<tr><td>Telescope </td>     <td> '+uvdata.header['telescop']+'</td>\n')
    wlog.write('</table><br>\n')
    notes_file = './{0}.notes.txt'.format(msinfo['run'])
    if os.path.isfile(notes_file):
        write_link_txt(wlog, notes_file, 'Notes and observing comments')
#    wlog.write('<br><small>(*) The channel and integration time  displayed are derived from the averaged data that has been processed in the pipeline. Higher time or frequency resolution may be available and can be obtained if required for the science extraction or advanced calibration techniques.  If this is required please contact the e-MERLIN science support team.</small>')
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def weblog_obssum(msinfo):
    ###### Observation summary page
    wlog = open(weblog_dir + "obs_summary.html","w")
    weblog_header(wlog, 'Observation summary', msinfo['run'])
    #------------------------------------------
    wlog.write('<h3>Summary:</h3>\n')
    listobs_file = info_link + msinfo['msfile'] + '.listobs.txt'
    write_link_txt(wlog, listobs_file, 'Summary of the observation (listobs)')
    wlog.write('<h3>Sources:</h3>\n')
    sepfile = info_link + 'source_separations.txt'
    if msinfo['sources']['targets'] != '':
        wlog.write('<table class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:30%">\n')
        wlog.write('<tr><td><b>Target</b> </td><td><b>Phase cal</b></td><td><b>Separation [deg]</b></td>\n')
        separations = msinfo['separations']
        for s1, s2 in zip(msinfo['sources']['targets'].split(','),
                      msinfo['sources']['phscals'].split(',')):
            if s1 not in msinfo['sources']['mssources'].split(','):
                logger.warning('{} not in MS'.format(s1))
                separation = 'Not in MS'
            elif s2 not in msinfo['sources']['mssources'].split(','):
                logger.warning('{} not in MS'.format(s2))
                separation = 'Not in MS'
            else:
                try:
                    separation = '{0:5.2f}'.format(separations[s1+'-'+s2])
                except:
                    try:
                        separation = '{0:5.2f}'.format(separations[s2+'-'+s1])
                    except:
                        separation = '0.0'
            wlog.write('<tr><td>{0} </td><td> {1}</td><td>{2}</td>\n'.format(s1, s2, separation))
        wlog.write('</table>\n')
    wlog.write(('Source pairs and separations: '
               '<a href="{0}" target="_blank">txt</a><br><br>\n').format(sepfile))
    wlog.write('<table  class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:30%">\n')
    wlog.write('Sources in MS:')
    wlog.write('<tr><td><b>Source</b> </td><td> <b>Intent</b></td>\n')
    for source in msinfo['sources']['mssources'].split(','):
        wlog.write('<tr><td>{0} </td><td> {1}</td>\n'.format(source, msinfo['sources']['source_intent'][source]))
    wlog.write('</table><br>\n')
    missing_sources = ', '.join([s for s in
                                 msinfo['sources']['allsources'].split(',') if s not in
                       msinfo['sources']['mssources']])
    if missing_sources != '':
        wlog.write('*Sources specified in the inputs file but not in MS: {0}\n'.format(missing_sources))
    else:
        wlog.write('*All sources specified in the inputs file are in the MS.\n')

    wlog.write('<h3>Antennas:</h3>\n')
    wlog.write('<pre>\n')
    for a in msinfo['antennas']:
        wlog.write("{0}\n".format(a))
    wlog.write('</pre>\n')
    wlog.write('Reference antenna: {}'.format(msinfo['refant']))
    wlog.write('<h3>Source elevation:</h3>\n')
    try:
        elev_plot = glob.glob(weblog_dir+'plots/plots_observation/{0}_elevation.png'.format(msinfo['msfilename']))
        wlog.write('<a href = ".{0}"><img style="max-width:600px" src=".{0}"></a><br>\n'.format(elev_plot[0]))
    except:
        pass
    wlog.write('<br><h3>UV coverage:</h3>\n')
    try:
        uvcov = np.sort(glob.glob(weblog_dir+'plots/plots_observation/{0}_uvcov_*.png'.format(msinfo['msfilename'])))
        for u in uvcov:
            wlog.write('<a href = ".{0}"><img style="max-width:700px" src=".{0}"></a><br>\n'.format(u))
            wlog.write('<hr>\n')
    except:
        pass
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()

def create_pnghtml_baselines(plots_path, source, subtitle, msinfo, datacolumn):
    page_path = weblog_dir +plots_path+'_'+source+".html"
    wlog = open(page_path,"w")
    weblog_header(wlog, source, msinfo['run'])
    wlog.write('<h3>{0}</h3>\n'.format(subtitle))
    #------------------------------------------

    plot_file = weblog_dir+'plots/'+plots_path+'/{0}_4plot_{1}_{2}'.format(msinfo['msfilename'],
                                                      source,
                                                      datacolumn)
    wlog.write('<table bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:40%">\n')
    wlog.write('<tr><td>{0}</td><td>{1}</td><td>{2}</td><td>{3}</td>\n'.format('Amp vs Time',
                                                                               'Phase vs  Time',
                                                                               'Amp vs Freq',
                                                                               'Phase vs Freq'))
    link0 = '<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(plot_file+'0.png')
    link1 = '<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(plot_file+'1.png')
    link2 = '<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(plot_file+'2.png')
    link3 = '<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(plot_file+'3.png')
    wlog.write('<tr>'+link0+link1+link2+link3+'</tr>')
    wlog.write('</table></td>\n')
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()
    return page_path

def plots_data(msinfo, wlog):
    wlog.write('<h3>Uncalibrated visibilities</h3>\n')
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines('plots_data', source,
                                             'Uncalibrated amplitude and phase against time and frequency.',
                                             msinfo, 'data')
        wlog.write('{0} <a href=".{1}" target="_blank">plots</a><br>\n'.format(source, page_path))

def plots_corrected(msinfo, wlog):
    wlog.write('<h3>Calibrated visibilities</h3>\n')
    for source in msinfo['sources']['mssources'].split(','):
        page_path = create_pnghtml_baselines('plots_corrected', source,
                                             'Calibrated amplitude and phase against time and frequency.',
                                             msinfo, 'corrected')
        wlog.write('{0} <a href=".{1}" target="_blank">plots</a><br>\n'.format(source, page_path))

#def plots_caltables(msinfo, wlog):
#    all_plots = np.sort(glob.glob('./weblog/plots/caltables/*png'))
#    for p in all_plots:
#        wlog.write('<a href=".{1}" target="_blank">{0}</a><br>\n'.format(os.path.basename(p), p))

def plots_uvplt(msinfo, wlog):
    wlog.write('<h3>Calibrated UVplots</h3>\n')
    all_plots = np.sort(glob.glob('./weblog/plots/plots_uvplt/*_uvplt_*png'))
    for p in all_plots:
        source_name = os.path.splitext(p)[0].split('_')[-1]
        wlog.write('<h4>{0}</h4>\n'.format(source_name))
        wlog.write('<table cellspacing = "0" cellpadding = "4px" style="width:40%">\n')
        wlog.write('<tr><td  valign="top"><a href = ".{0}"><img style="max-width:960px" src=".{0}"></a></td>\n'.format(p))
        if source_name in msinfo['sources']['calsources'].split(','):
            p_model = './weblog/plots/plots_uvplt/{0}_uvpltmodel_{1}.png'.format(msinfo['msfilename'],source_name)
            wlog.write('<td><a href = ".{0}"><img style="max-width:960px" src=".{0}"></a></td>\n'.format(p_model))
        wlog.write('</tr></table><br><br>\n<hr>\n')

def write_caltable(caltable, wlog):
    wlog.write('<table class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:40%">\n')
    #wlog.write('<tr><td>{} </td>   <td>{}</td>\n'.format('Parameter', 'Value'))
    for name, value in caltable.items():
        if type(value) == str:
            value = value.replace(',', ', ')
        wlog.write('<tr><td>{0}</td><td>{1}</td>\n'.format(name, value))
    wlog.write('</table></td>\n')

def table_colors(steps, step, prev_steps):
    t_fmt = '%Y-%m-%d %H:%M:%S'
    neutral = '#EEEEEE'
    red = '#FF6347'
    green = '#008B45'
    outdated = {True: red, False: green}
    if steps[step][0] == 0:
        color = neutral
    elif type(steps[step][0]) is str:
        time_step = datetime.datetime.strptime(steps[step][0], t_fmt)
        if prev_steps != []:
            latest_prev = max([datetime.datetime.strptime(steps[stepi][0], t_fmt)
                               for stepi in prev_steps])
            color = outdated[time_step < latest_prev]
        else:
            color = green
    else:
        color = '#000000'
    return color

def elapsed_time(delta_min):
    if delta_min > 0.0 and delta_min < 1:
        delta_t = '<1 min'
    elif delta_min >= 1 and delta_min < 50:
        delta_t = '{0:2.0f} min'.format(delta_min)
    elif delta_min >= 50:
        delta_t = '{0:4.1f} h'.format(delta_min/60.)
    else:
        delta_t = '-'
    return delta_t


def table_steps(eMCP):
    table_txt = '<br><h4>Execution summary</h4>\n'
    table_txt += ('<table class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:80%">\n')
    table_txt += ("<tr><th style='width: 15%;'>Step</th>" \
                  "<th style='width: 4%;'>Code</th>" \
                  "<th style='width: 25%;'>Execution ended</th>" \
                  "<th style='width: 15%;'>Execution time</th>" \
                  "<th style='width: 70%;'>Notes</th></tr>\n")
    prev_steps = []
    nosteps = ['plot_data', 'save_flags', 'plot_corrected', 'first_images']
    for step, step_info in eMCP['steps'].items():
        time_s, delta_min, msg = step_info
        delta_t = elapsed_time(delta_min)
        if step == 'start_pipeline':
            delta_t = '-'
        color = table_colors(eMCP['steps'], step, prev_steps)
        table_txt += ('<tr><td>{0}</td>'.format(step))
        table_txt += ('<td bgcolor={0}></td>\n'.format(color))
        table_txt += ('<td align="center">{0}</td>\n'.format(time_s))
        table_txt += ('<td align="center">{0}</td>\n'.format(delta_t))
        table_txt += ('<td>{0}</td></tr>\n'.format(msg))
        if (eMCP['steps'][step][0] != 0) and (step not in nosteps):
            prev_steps.append(step)
    table_txt += ('</table>\n')
    return table_txt

def weblog_pipelineinfo(eMCP):
    msinfo = eMCP['msinfo']
    # Summary of pipeline execution page
    wlog = open(weblog_dir + "pipelineinfo.html","w")
    weblog_header(wlog, 'Pipeline info', msinfo['run'])
    #------------------------------------------         
    wlog.write('CASA version: {}\n<br>'.format(eMCP['casa_version']))
    wlog.write('Pipeline version: {}\n<br>'.format(eMCP['pipeline_version']))
    wlog.write(table_steps(eMCP))
    wlog.write('Green = executed<br>')
    wlog.write('Red = executed but outdated by a previous step<br>')
    wlog.write('<br><h4>Relevant log files:</h4>\n')
    write_link_txt(wlog, info_link + 'eMCP.log.txt', 'Pipeline log', text='eMCP.log')
    write_link_txt(wlog, info_link + 'casa_eMCP.log.txt', 'CASA log', text='casa_eMCP.log')
    # eMCP_info dictionary as text
    prt_dict_tofile(eMCP, tofilename=info_dir + 'eMCP_info.txt', pre='  ')
    eMCP_txtfile = info_link + 'eMCP_info.txt'
    wlog.write('<br><h4>Relevant parameter files:</h4>\n')
    write_link_txt(wlog, eMCP_txtfile, 'Pipeline info (dict)', text='eMCP_info.txt')
    # caltables dictionary as text
    try:
        caltables = load_obj(calib_dir+'caltables.pkl')
        prt_dict_tofile(caltables, tofilename=info_dir + 'caltables.txt', pre='  ')
        caltables_txtfile = info_link + 'caltables.txt'
        write_link_txt(wlog, caltables_txtfile, 'Calibration info (dict)', text='caltables.txt')
    except:
        pass
    #------------------------------------------                                                                       
    weblog_foot(wlog)
    wlog.close()

def write_fluxscale(wlog, msinfo):
    p = glob.glob(calib_dir + '{0}_fluxscale.png'.format(msinfo['msfilename']))[0]
    wlog.write('<td><a href = ".{0}"><img style="max-width:700px" src=".{0}"></a></td>\n'.format(p))

    with open(calib_dir + 'allcal_ap.G1_fluxes.txt') as f:
        lines = f.readlines()
    wlog.write('\n<br><br>\n')
    wlog.write('CASA fluxscale output (not corrected by eMfactor):')
    wlog.write('\n<pre>')
    for line in lines:
        wlog.write(line)
    wlog.write('\n</pre>\n<br>')

def write_applycal_table(wlog, eMCP, ap_dict='applycal_dict'):
    applycal_dict = eMCP[ap_dict]
    wlog.write('<table class="table1" bgcolor="#eeeeee" border="3px"'\
               'cellspacing = "0" cellpadding = "4px" style="width:60%">\n')
    wlog.write("<tr><th>Source</th>" \
          "<th>Table</th>" \
          "<th>Gainfield</th>" \
          "<th>spwmap</th>" \
          "<th>Interpolation</th></tr>\n")
    for source, tables in applycal_dict.items():
        for i, table in enumerate(tables.keys()):
            if i != 0:
                source = ''
            wlog.write('<tr>'\
                       '<td>{0}</td>'\
                       '<td>{1}</td>'\
                       '<td>{2}</td>'\
                       '<td>{3}</td>'\
                       '<td>{4}</td>\n'.format(source, table,
                                               *tables[table]))
    wlog.write('</table><br><br>\n')

def weblog_calibration(eMCP):
    msinfo = eMCP['msinfo']
    ###### Calibration page
    wlog = open(weblog_dir + "calibration.html","w")
    weblog_header(wlog, 'Calibration', msinfo['run'])
    #------------------------------------------
    wlog.write('<a href="#applycal">{0}</a><br>\n'.format('Applycal summary'))
    all_calsteps = ['bpcal_d.K0','bpcal_p.G0','bpcal_ap.G0','bpcal.BP0',
                    'allcal_d.K1','allcal_p.G1','allcal_ap.G1',
                    'fluxscale',
                    'bpcal.BP2',
                    'allcal_p.G3','allcal_ap.G3',
                    'phscal_p_scan.G3','phscal_ap_scan.G3'
                    ]
    if eMCP['is_mixed_mode']:
        all_calsteps.append('narrow_p_offset.G3')
        all_calsteps.append('narrow_bpcal.BP2')
    if os.path.isfile('./weblog/calib/caltables.pkl'):
        caltables = load_obj('./weblog/calib/caltables.pkl')
        for calstep in all_calsteps:
            logger.debug('calstep {}'.format(calstep))
            try:
                if calstep == 'fluxscale' and os.path.isfile(calib_dir + 'allcal_ap.G1_fluxes.txt'):
                    wlog.write('<h4>{}</h4>\n'.format('fluxscale'))
                    write_fluxscale(wlog, msinfo)
                    wlog.write('\n<hr>\n')
                else:
                    wlog.write('<h4>{}</h4>\n'.format(caltables[calstep]['name']))
                    wlog.write('<table cellspacing = "20" cellpadding = "4px" style="width:40%">\n<tr> <td valign="top">\n')
                    write_caltable(caltables[calstep], wlog)
                    all_plots = np.sort(glob.glob('./weblog/plots/caltables/*{}*.png'.format(calstep)))
                    for p in all_plots:
                        wlog.write('<td><a href = ".{0}"><img style="max-width:700px" src=".{0}"></a></td>\n'.format(p))
                    wlog.write('</table><br><br>\n<hr>\n')
            except:
                pass
    try:
        wlog.write('<h2 id="applycal">{}</h2>\n'.format('Applycal'))
        wlog.write('List of tables and apply parameters used to correct each source:')
        write_applycal_table(wlog, eMCP, ap_dict='applycal_dict')

        if eMCP['is_mixed_mode']:
            wlog.write('Tables applied to narrow band data:')
            write_applycal_table(wlog, eMCP, ap_dict='applycal_dict_sp')
    except:
        pass
    wlog.write('<a href="#top"> <small>(up)</small></a><br>\n')
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def weblog_flagstats(msinfo):
    ###### Flag statistics page
    wlog = open(weblog_dir + "flagstats.html","w")
    weblog_header(wlog, 'Flag statistics', msinfo['run'])
    #------------------------------------------
    wlog.write('High-resolution plots in step title.<br><br>')
    flagstats_steps = ['import','aoflagger',
                       'apriori',
                       'flag_manual',
                       'restore_flags',
                       'flag_manual_avg',
                       'initial_bpcal',
                       'initial_gaincal',
                       'applycal_all',
                       'flag_target']
    prev_perc_flagged = 0.0
    for step in flagstats_steps:
        try:
            flag_stats = load_obj('./weblog/plots/plots_flagstats/flagstats_{}.pkl'.format(step))
            perc_flagged = flag_stats['flagged']/flag_stats['total']*100.
            diff_flagged = perc_flagged-prev_perc_flagged
            prev_perc_flagged = perc_flagged
            all_plots = np.sort(glob.glob('./weblog/plots/plots_flagstats/*{}.png'.format(step)))
            if len(all_plots) > 0:
                wlog.write('<a href=".{0}" target="_blank">{1}</a> (Total: '
                           '{2:3.1f}%. Increase: {3:3.1f}%)<br>\n'.format(all_plots[1],
                                                                step,
                                                                perc_flagged,
                                                                diff_flagged))
                wlog.write('<a href = ".{0}"><img style="max-width:1200px" src=".{0}"></a><br>\n'.format(all_plots[0]))
        except:
            pass
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def weblog_plots(msinfo):
    ###### Plots page
    wlog = open(weblog_dir + "plots.html","w")
    weblog_header(wlog, 'Plots', msinfo['run'])
    #------------------------------------------
    if (os.path.isdir('./weblog/plots/plots_data/')) and \
       (os.listdir('./weblog/plots/plots_data/')):
        plots_data(msinfo, wlog)
    if (os.path.isdir('./weblog/plots/plots_corrected/')) and \
       (os.listdir('./weblog/plots/plots_corrected/')):
        plots_corrected(msinfo, wlog)
    if (os.path.isdir('./weblog/plots/plots_uvplt/')) and \
       (os.listdir('./weblog/plots/plots_uvplt/')):
        plots_uvplt(msinfo, wlog)
#    if (os.path.isdir('./weblog/plots/plots_flagstats/')) and \
#    (os.listdir('./weblog/plots/plots_flagstats/')):
#        plots_flagstats(msinfo, wlog)
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()

def check_any_image(eMCP):
    msinfo = eMCP['msinfo']
    images = False
    for i in range(len(msinfo['sources']['targets'].split(','))):
        img_dir = weblog_dir + 'images/{0}/'.format(msinfo['sources']['targets'].split(',')[i])
        if (os.path.isdir(img_dir)) and (os.listdir(img_dir)):
            images = True
    return images

def weblog_images(eMCP):
    msinfo = eMCP['msinfo']
    ###### Images page
    wlog = open(weblog_dir + "images.html","w")
    weblog_header(wlog, 'Crude images', msinfo['run'])
    #------------------------------------------
    if check_any_image(eMCP):
        wlog.write('Jump to target:<br>')
        for target in msinfo['sources']['targets'].split(','):
            wlog.write('<a href="#{0}">{0}</a><br>\n'.format(target))
        note_text = 'Note: These are crude images of targets and phase '\
        'calibrators. The images are produced automatically with tclean '\
        'in CASA using auto-thresh mode, which generates cleaning boxes '\
        'without human intervention. These images should not be used for production.'
        wlog.write('<br><div style="max-width:800px; word-wrap:break-word;">{0}</div>'.format(note_text))
    else:
        wlog.write('No images found')

    for i in range(len(msinfo['sources']['targets'].split(','))):
        img_dir = weblog_dir + 'images/{0}/'.format(msinfo['sources']['targets'].split(',')[i])
        if (os.path.isdir(img_dir)) and (os.listdir(img_dir)):
            show_image(eMCP, wlog, i)
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()

def show_image(eMCP, wlog, i):
    msinfo = eMCP['msinfo']
    nterms = eMCP['defaults']['first_images']['nterms']
    num = 0
    if nterms == 1:
        ext = ''
    else:
        ext = '.tt0'
    target = msinfo['sources']['targets'].split(',')[i]
    phscal = msinfo['sources']['phscals'].split(',')[i]
    try:
        peak_target, noise_target, scaling_target = eMCP['img_stats'][target]
        peak_phscal, noise_phscal, scaling_phscal = eMCP['img_stats'][phscal]
    except:
        peak_target, noise_target, scaling_target = 0.0, 0.0, 0.0
        peak_phscal, noise_phscal, scaling_phscal = 0.0, 0.0, 0.0
    img_target = weblog_dir+'images/{0}/{1}_{0}_img{2:02d}'.format(target, msinfo['msfilename'], num)
    img_phscal = weblog_dir+'images/{0}/{1}_{0}_img{2:02d}'.format(phscal, msinfo['msfilename'], num)
    wlog.write('<hr>\n')
    wlog.write('<h3 id="{0}">{0} <a href="#top"> <small>(up)</small></a></h3>\n'.format(target))
#    wlog.write('<a href="#images"> (up)</a><br>\n')
    wlog.write('<table class="table1" bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:70%">\n')
    wlog.write('<tr><td><b>{0}</b> (Target image) ' \
               'Peak: {1:3.3f} mJy ' \
               '(scaling: {2:4.1f})</td>' \
               '<td>(Target residual) rms: {3:5.3f} mJy</td>' \
               '<td>zoom image</td><td>zoom residual</td></tr>\n'.format(target,
                                                                 peak_target*1000.,
                                                                 scaling_target,
                                                                 noise_target*1000.))
    wlog.write('<tr>')
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(img_target+'.image'+ext+'.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>\n'.format(img_target+'.residual'+ext+'.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(img_target+'.image'+ext+'_zoom.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>\n'.format(img_target+'.residual'+ext+'_zoom.png'))
    wlog.write('</tr>')
    wlog.write('<tr><td><b>{0}</b> (Phasecal image) ' \
               'Peak: {1:5.1f} mJy ' \
               '(scaling: {2:4.1f})</td>' \
               '<td>(Phasecal residual) rms: {3:5.3f} mJy</td>' \
               '<td>zoom image</td><td>zoom residual</td></tr>\n'.format(phscal,
                                                               peak_phscal*1000.,
                                                               scaling_phscal,
                                                               noise_phscal*1000.))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(img_phscal+'.image'+ext+'.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>\n'.format(img_phscal+'.residual'+ext+'.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>'.format(img_phscal+'.image'+ext+'_zoom.png'))
    wlog.write('<td><a href = ".{0}"><img style="max-width:600px" src=".{0}"></a></td>\n'.format(img_phscal+'.residual'+ext+'_zoom.png'))
    wlog.write('</table><br>\n')

def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass


def weblog_download(msinfo):
    ###### Images page
    wlog = open(weblog_dir + "download.html","w")
    weblog_header(wlog, 'Download data', msinfo['run'])
    #------------------------------------------
    wlog.write('This tar file contains the MS and all the plots in the weblog:<br>')
    filepath = '../{}.tar'.format(msinfo['run'])
    wlog.write('{0}: <a href="../{1}" target="_blank">tar</a><br>\n'.format(msinfo['run'], filepath))
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def start_weblog(eMCP, silent=False):
    msinfo = eMCP['msinfo']
    logger.info('Updating weblog')
    ###  Start weblog  ###
    weblog_index(msinfo)
    weblog_obssum(msinfo)
    weblog_pipelineinfo(eMCP)
    weblog_calibration(eMCP)
    weblog_plots(msinfo)
    weblog_flagstats(msinfo)
    weblog_images(eMCP)
    weblog_download(msinfo)



