import os


def weblog_header(wlog, section):
    wlog.write('<html>\n')
    wlog.write('<head>\n')
    wlog.write('<title>e-MERLIN Pipeline Web Log</title>\n')
    wlog.write('</head>\n')
    wlog.write('<body>\n')
    wlog.write('<h5>e-MERLIN Pipeline Web Log</h5>\n')
    wlog.write('<form>\n')
    wlog.write('<input type="button" value="Home"  onclick="window.location.href=\'./index.html\'">\n')
    wlog.write('<input type="button" value="Observation summary"  onclick="window.location.href=\'./obs_summary.html\'">\n')
    wlog.write('<input type="button" value="Plots"  onclick="window.location.href=\'./plots.html\'">\n')
    wlog.write('<input type="button" value="Files"  onclick="window.location.href=\'./files.html\'">\n')
    wlog.write('</form>\n')
    #wlog.write('<br>\n')
    wlog.write('<h1>{0}</h1>\n'.format(section))
    wlog.write('<hr>\n')


def weblog_foot(wlog):
    wlog.write('<hr><br>\n')
    wlog.write('<small><a href="http://www.e-merlin.ac.uk/observe/pipeline/">e-MERLIN Pipeline</a><br>\n')
    wlog.write('<a href="https://github.com/mkargo/pipeline">e-MERLIN Pipeline Github</a><br>\n')
    wlog.write('<a href="http://www.e-merlin.ac.uk/">e-MERLIN</a><br>\n')
    wlog.write('<a href="http://www.e-merlin.ac.uk/data_red/">User support and data reduction for e-MERLIN</a><br>\n')
    wlog.write('<a href="http://www.jive.nl/jivewiki/doku.php?id=parseltongue:parseltongue">Parseltongue</a></small><br>\n')
    wlog.write('<hr>\n')
    wlog.write('</body>\n')
    wlog.write('</html>\n')

def write_link_txt(wlog, infile, intext):
    wlog.write('{0}: <a href="{1}" target="_blank">txt</a><br>\n'.format(intext, infile))

def weblog_index(msinfo):
    wlog = open("weblog/index.html","w")
    weblog_header(wlog, 'Home')
    #------------------------------------------
    wlog.write('<table bgcolor="#eeeeee" border="3px" cellspacing = "0" cellpadding = "4px" style="width:30%">\n')
    wlog.write('<tr><td>Project </td>   <td> {}</td>\n'.format(msinfo['project']))
    wlog.write('<tr><td>Run </td>   <td> {}</td>\n'.format(msinfo['run']))
    #wlog.write('<tr><td>Date observed </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].date()))
    wlog.write('<tr><td>Start </td>  <td> {0} </td>\n'.format(msinfo['t_ini'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>End </td>    <td> {0} </td>\n'.format(msinfo['t_end'].strftime("%Y-%m-%d %H:%M")))
    wlog.write('<tr><td>Band </td>   <td> {0} </td>\n'.format(msinfo['band']))
    wlog.write('<tr><td>Antennas </td>   <td> {0} </td>\n'.format(', '.join(msinfo['antennas'])))
    wlog.write('<tr><td>Number of sources </td>   <td> {0} </td>\n'.format(len(msinfo['mssources'].split(','))))
    wlog.write('<tr><td>Frequency </td>  <td> {0:5.2f} - {1:5.2f} GHz </td>\n'.format(
        msinfo['freq_ini'], msinfo['freq_end']))
    wlog.write('<tr><td>Num. spw </td>             <td> {0} </td>\n'.format(msinfo['num_spw']))
    wlog.write('<tr><td><abbr title="Number of channels per spw after averaging. The raw data may provide more channels">Channels/spw<sup>*</sup></abbr></td><td> {0} </td>\n'.format(msinfo['nchan']))
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
    write_link_txt(wlog, './{0}.notes.txt', 'Notes and observing comments'.format(msinfo['project']))
    wlog.write('<br><small>(*) The channel and integration time  displayed are derived from the averaged data that has been processed in the pipeline. Higher time or frequency resolution may be available and can be obtained if required for the science extraction or advanced calibration techniques.  If this is required please contact the e-MERLIN science support team.</small>')

    #------------------------------------------
    weblog_foot(wlog)
    wlog.close()


def makedir(pathdir):
    try:
        os.mkdir(pathdir)
#        logger.info('Create directory: {}'.format(pathdir))
    except:
#        logger.debug('Cannot create directory: {}'.format(pathdir))
        pass



def start_weblog(msinfo):
    ###  Start weblog  ###
    try: os.mkdir('weblog')
    except: pass

    weblog_index(msinfo)


def make_uvplot(msfile, msinfo, plot_dir):
    uvplt_dir = plot_dir+'/uvplt/'
    makedir(uvplt_dir)
    for f in msinfo['mssources'].split(','):
        print(f)
        outname = '{0}{1}_{2}_{3}.png'.format(uvplt_dir, msinfo['run'], 'uvplt', f)
        plotuv(vis=msfile, field=f, figfile=outname, symb='.')



