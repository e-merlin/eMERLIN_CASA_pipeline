import numpy as np
import sys, os
from callibrary import applycaltocallib

def write_wsclean_command(msfile, imagename, imsize=512, cellsize=None, field='0',
                          gain=0.05, mgain=0.85, robust=0.5, niter=10000,
                          automask=5, autothreshold=3, mask='auto', casamask='',
                          multiscale=None, multiscale_gain=0.1,
                          multiscale_scale_bias=0.8,
                          multiscale_threshold_bias=0.9,
                          joinchannels=None, channelsout=3,
                          column='CORRECTED_DATA', **kwargs):
    if cellsize == None:
        cellsize = find_cellsize(msfile)
    command = []
    command.append('time wsclean')
    command.append('-field {0}'.format(field))
    command.append('-datacolumn {0}'.format(column))
    command.append('-scale {0}arcsec'.format(cellsize))
    command.append('-size {0} {0}'.format(int(imsize*1.1)))
    command.append('-trim {0} {0}'.format(imsize))
    command.append('-gain {0}'.format(gain))
    command.append('-mgain {0}'.format(mgain))
    command.append('-weight briggs {0}'.format(robust))
    if mask == 'auto':
        command.append('-auto-mask {0}'.format(automask))
    elif mask == 'casa':
        command.append('-casamask {0}'.format(casamask))
    if multiscale:
        command.append('-multiscale') # -multiscale-gain 0.1 -multiscale-scale-bias 0.3
        command.append('-multiscale-gain {0}'.format(multiscale_gain))
        command.append('-multiscale-scale-bias {0}'.format(multiscale_scale_bias))
#        command.append('-multiscale-threshold-bias {0}'.format(multiscale_threshold_bias))
    if joinchannels:
        command.append('-joinchannels -channelsout {0}'.format(channelsout))
    command.append('-niter {0}'.format(niter))
    command.append('-auto-threshold {0}'.format(autothreshold))
    command.append('-name {0}'.format(imagename))
    command.append('{0}'.format(msfile))
    outcommand = '  '.join(command)
    print(outcommand)
    return outcommand


def run_wsclean(num, msfile, **kwargs):
    imagename = './selfcalimages/img.{0:02d}'.format(num)
    wsclean_command = write_wsclean_command(msfile, imagename, **kwargs)
    os.system(wsclean_command)
    print('Finished command: \n'+wsclean_command)

def check_band(msfile):
    # Take first frequency in the MS
    ms.open(msfile)
    freq = ms.getdata(['axis_info'])['axis_info']['freq_axis']['chan_freq'][0][0]/1e9
    ms.close()
    band = ''
    if (freq > 1.2) and (freq < 1.7):
        band = 'L'
    if (freq > 4) and (freq < 8):
        band = 'C'
    if (freq > 22) and (freq < 24):
        band = 'K'
    return band

def find_cellsize(msfile):
    band = check_band(msfile)
    if band == 'L':
        cellsize = 0.02
    elif band == 'C':
        cellsize = 0.006
    elif band == 'K':
        cellsize = 0.002
    return cellsize

def run_gaincal(msfile, caltable, calmode, solint, combine, callib, refant):
    print('Starting gaincal')
    gaincal(msfile, caltable, gaintype = 'G', calmode =calmode, solint=solint, uvrange='',
              antenna='', combine=combine, spw='', refant=refant,
              docallib = True, callib=callib)

def write_spwmap(msfile, combine):
    if combine == '':
        spwmap = []
    elif combine == 'spw':
        num_spw = len(vishead(msfile, mode = 'list', listitems = ['spw_name'])['spw_name'][0])
        spwmap = [0]*num_spw
    return spwmap

def run_flagdata(msfile):
    timedevscale = 5
    freqdevscale = 5
    flagdata(vis=msfile, mode='rflag', field='',
         antenna='', scan='',spw='', correlation='',
         ntime='30min', combinescans=True, datacolumn='corrected',
         timedevscale=timedevscale, freqdevscale=freqdevscale,
         action='apply', display='', flagbackup=True)


def selfcal_step(num, msfile, calmode, solint, combine='', **kwargs):
    callib = './selfcal/callib.txt'
    caltable = './selfcal/gcal.{0:02d}'.format(num)
    # Calibration
    refant = kwargs.pop('refant', '')
    run_gaincal(msfile, caltable, calmode, solint, combine, callib, refant)
    applycaltocallib(filename=callib, append=True, gaintable=caltable,
                     spwmap=write_spwmap(msfile, combine), calwt=True)
    # Apply calibration
    applycal(vis=msfile, docallib=True,callib=callib)
    # Flagging
    if calmode =='ap':
        run_flagdata(msfile)
    # Imaging
    imagename = './selfcalimages/img.{0:02d}'.format(num)
    run_wsclean(num, msfile, **kwargs)

def run_selfcal(msfile, kwargs={}):
    os.system('mkdir selfcalimages')
    os.system('mkdir selfcal')
    if not os.path.isfile('./selfcal/callib.txt'):
        os.system('touch ./selfcal/callib.txt')
    else:
        print('./selfcal/callib.txt found! Will use previous calibrations. Remove file to avoid that.')

    # Initial Image
    run_wsclean(0, msfile, column='DATA', automask=10, autothreshold=8, **kwargs)
    # Selfcalibration cycles:
    selfcal_step(1, msfile=msfile, calmode='p',  solint='16s', automask=8, autothreshold=5,  **kwargs)
    selfcal_step(2, msfile=msfile, calmode='p',  solint='16s', automask=6, autothreshold=5, **kwargs)
    selfcal_step(3, msfile=msfile, calmode='p',  solint='1s', combine='spw', automask=6, autothreshold=4, **kwargs)
    selfcal_step(4, msfile=msfile, calmode='ap', solint='inf', automask=5, **kwargs)
    selfcal_step(5, msfile=msfile, calmode='p',  solint='16s', **kwargs)
    selfcal_step(6, msfile=msfile, calmode='ap', solint='32s', **kwargs)
    selfcal_step(7, msfile=msfile, calmode='p',  solint='8s',  **kwargs)
    selfcal_step(8, msfile=msfile, calmode='ap', solint='16s', **kwargs)
    selfcal_step(9, msfile=msfile, calmode='p',  solint='4s',  **kwargs)


msfile='/mirror2/scratch/jmoldon/emerlin/CY3229/CY3229_L_20160117/splits/tar_1703+2615'
flagmanager(vis=msfile, mode='restore', versionname='applycal_1')

#global_args = {'imsize':1024,
#               'mask':'casa',
#               'casamask':'/mirror2/scratch/jmoldon/emerlin/CY3229/CY3229_L_20160117/logs/kk.mask',
#               'refant':'Da'}
global_args = {'imsize':1024, 'refant':'Da'}

run_selfcal(msfile, kwargs=global_args)

