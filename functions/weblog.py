import os
import numpy as np
import glob
from taskinit import *
from tasks import *

import logging
logger = logging.getLogger('logger')


def weblog_header(wlog, section, project):
    wlog.write('<html>\n')
    wlog.write('<head>\n')
    wlog.write('<title>{} - e-MERLIN Pipeline Web Log</title>\n'.format(project))
    wlog.write('</head>\n')
    wlog.write('<body>\n')
    wlog.write('<h4>e-MERLIN Pipeline Web Log<br>{}</h4>\n'.format(project))
    wlog.write('<form>\n')
    wlog.write('<input type="button" value="Home"  onclick="window.location.href=\'./index.html\'">\n')
    wlog.write('<input type="button" value="Observation summary"  onclick="window.location.href=\'./obs_summary.html\'">\n')
    wlog.write('<input type="button" value="Plots"  onclick="window.location.href=\'./plots.html\'">\n')
    wlog.write('<input type="button" value="Files"  onclick="window.location.href=\'./files.html\'">\n')
    wlog.write('</form>\n')
    #wlog.write('<br>\n')
    wlog.write('<h2>{0}</h2>\n'.format(section))
    wlog.write('<hr>\n')


def weblog_foot(wlog):
    wlog.write('<hr><br>\n')
    wlog.write('<small>')
    wlog.write('<a href="https://github.com/e-merlin/CASA_eMERLIN_pipeline">e-MERLIN Pipeline Github</a><br>\n')
    wlog.write('<a href="http://www.e-merlin.ac.uk/">e-MERLIN</a><br>\n')
    wlog.write('<a href="http://www.e-merlin.ac.uk/data_red/">User support and data reduction for e-MERLIN</a><br>\n')
    wlog.write('<a href="https://casa.nrao.edu/">CASA</a></small><br>\n')
    wlog.write('<hr>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext):
    wlog.write('{0}: <a href="{1}" target="_blank">txt</a><br>\n'.format(intext, infile))

def weblog_index(msinfo):
    wlog = open("./weblog/index.html","w")
    weblog_header(wlog, 'Home', msinfo['run'])
    #------------------------------------------
    wlog.write('<table bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:40%">\n')
    wlog.write('<tr><td>Project </td>   <td> {}</td>\n'.format(msinfo['project']))
    wlog.write('<tr><td>Run </td>   <td> {}</td>\n'.format(msinfo['run']))
    #wlog.write('<tr><td>Date observed </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].date()))
    wlog.write('<tr><td>Start </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>End </td>    <td> {0} </td>\n'.format(msinfo['t_end'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>Band </td>   <td> {0} </td>\n'.format(msinfo['band']))
    wlog.write('<tr><td>Antennas </td>   <td> {0} </td>\n'.format(', '.join(msinfo['antennas'])))
    wlog.write('<tr><td>Number of sources </td>   <td> {0} </td>\n'.format(len(msinfo['sources']['mssources'].split(','))))
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
    write_link_txt(wlog, './{0}.notes.txt'.format(msinfo['run']), 'Notes and observing comments')
#    wlog.write('<br><small>(*) The channel and integration time  displayed are derived from the averaged data that has been processed in the pipeline. Higher time or frequency resolution may be available and can be obtained if required for the science extraction or advanced calibration techniques.  If this is required please contact the e-MERLIN science support team.</small>')

    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def weblog_obssum(msinfo):
    ###### Observation summary page
    wlog = open("./weblog/obs_summary.html","w")
    weblog_header(wlog, 'Observation summary', msinfo['run'])
    #------------------------------------------
    wlog.write('<h3>Summary:</h3>\n')
    write_link_txt(wlog, '../{}'.format(msinfo['msfile']+'.listobs'), 'Summary of the observation (listobs)')
    wlog.write('<h3>Sources:</h3>\n')
    wlog.write('<table bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:30%">\n')
    wlog.write('<tr><td><b>Source in MS</b> </td><td> <b>Intent</b></td>\n')
    for source in msinfo['sources']['mssources'].split(','):
        wlog.write('<tr><td>{0} </td><td> {1}</td>\n'.format(source, msinfo['sources']['source_intent'][source]))
    wlog.write('</table><br>\n')
    missing_sources = ', '.join([s for s in
                                 msinfo['sources']['allsources'].split(',') if s not in
                       msinfo['sources']['mssources']])
    if missing_sources != '':
        wlog.write('Sources in inputs file but not in MS: {0}'.format(missing_sources))
    else:
        wlog.write('All sources in inputs file are in the MS.')
    wlog.write('<h3>Antennas:</h3>\n')
    wlog.write('<pre>\n')
    for a in msinfo['antennas']:
        wlog.write("{0}\n".format(a))
    wlog.write('</pre>\n')
    wlog.write('<h3>Source elevation:</h3>\n')
    try:
        elev_plot = glob.glob('./plots/plots_observation/{0}_elevation.png'.format(msinfo['msfilename']))
        wlog.write('<a href = ".{0}"><img style="max-width:700px" src=".{0}"></a><br>\n'.format(elev_plot[0]))
    except:
        pass
    wlog.write('<br><h3>UV coverage:</h3>\n')
    try:
        uvcov = np.sort(glob.glob('./plots/plots_observation/{0}_uvcov_*.png'.format(msinfo['msfilename'])))
        for u in uvcov:
            wlog.write('<a href = ".{0}"><img style="max-width:700px" src=".{0}"></a><br>\n'.format(u))
            wlog.write('<hr>\n')
    except:
        pass
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()

def create_pnghtml_baselines(plots_path, source, subtitle, msinfo):
    page_path = "./weblog/"+plots_path+'_'+source+".html"
    wlog = open(page_path,"w")
    weblog_header(wlog, source, msinfo['run'])
    wlog.write('<h3>{0}</h3>\n'.format(subtitle))
    #------------------------------------------
    for bsl in msinfo['baselines']:
        try:
            b1, b2 = bsl.split('&')
            png = np.sort(glob.glob('./plots/{0}/{1}_4plot_{2}_{3}-{4}_*.png'.format(plots_path,
                                                            msinfo['msfilename'],
                                                            source,
                                                            b1, b2)))[0]
            wlog.write('<h5>{0}-{1}</h5>\n'.format(b1,b2))
            wlog.write('<a href = ".{0}"><img style="max-width:960px" src=".{0}"></a><br>\n'.format(png))
            wlog.write('<hr>')
        except:
            pass
    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()
    return page_path

def plots_data(msinfo, wlog):
    for source in msinfo['sources']['allsources'].split(','):
        page_path = create_pnghtml_baselines('plots_data', source, 'Uncalibrated amplitude and phase against time and frequency.', msinfo)
        wlog.write('{0} <a href=".{1}" target="_blank">plots</a><br>\n'.format(source, page_path))


def plots_corrected(msinfo, wlog):
    for source in msinfo['sources']['allsources'].split(','):
        page_path = create_pnghtml_baselines('plots_corrected', source, 'Calibrated amplitude and phase against time and frequency.', msinfo)
        wlog.write('{0} <a href=".{1}" target="_blank">plots</a><br>\n'.format(source, page_path))

def plots_caltables(msinfo, wlog):
    all_plots = np.sort(glob.glob('./plots/caltables/*png'))
    for p in all_plots:
        wlog.write('<a href=".{1}" target="_blank">{0}</a><br>\n'.format(os.path.basename(p), p))

def plots_uvplt(msinfo, wlog):
    all_plots = np.sort(glob.glob('./plots/plots_uvplt/*png'))
    for p in all_plots:
        source_name = os.path.splitext(p)[0].split('_')[-1]
        wlog.write('<h4>{0}</h4>\n'.format(source_name))
        wlog.write('<a href = ".{0}"><img style="max-width:960px" src=".{0}"></a><br>\n'.format(p))
        wlog.write('<hr>\n')



def weblog_plots(msinfo):
    ###### Plots page
    wlog = open("./weblog/plots.html","w")
    weblog_header(wlog, 'Plots', msinfo['run'])
    #------------------------------------------
    wlog.write('<h3>Uncalibrated visibilities</h3>\n')
    if (os.path.isdir('./plots/plots_data/')) and (os.listdir('./plots/plots_data/')):
        plots_data(msinfo, wlog)

    wlog.write('<h3>Calibrated visibilities</h3>\n')
    if (os.path.isdir('./plots/plots_corrected/')) and (os.listdir('./plots/plots_corrected/')):
        plots_corrected(msinfo, wlog)

    wlog.write('<h3>Calibration tables</h3>\n')
    if (os.path.isdir('./plots/caltables/')) and (os.listdir('./plots/caltables/')):
        plots_caltables(msinfo, wlog)

    wlog.write('<h3>Calibrated UVplots</h3>\n')
    if (os.path.isdir('./plots/plots_uvplt/')) and (os.listdir('./plots/plots_uvplt/')):
        plots_uvplt(msinfo, wlog)


    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def makedir(pathdir):
    try:
        os.mkdir(pathdir)
        logger.info('Create directory: {}'.format(pathdir))
    except:
        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass



def start_weblog(msinfo):
    logger.info('Start weblog')
    ###  Start weblog  ###
    try: os.mkdir('weblog')
    except: pass

    weblog_index(msinfo)
    weblog_obssum(msinfo)
    weblog_plots(msinfo)
    logger.info('End weblog')




