import os
import glob
import logging
import pickle

logger = logging.getLogger('logger')

def write_fluxscale(wlog, msinfo):
    """Display flux calibration results in a modern format"""
    calfluxes_file = './weblog/calib/calfluxes.yaml'
    if os.path.isfile(calfluxes_file):
        try:
            with open(calfluxes_file, 'rb') as fp:
                calfluxes = pickle.load(fp)
            
            wlog.write('<div class="subsection">\n')
            wlog.write('  <h3>Flux calibration results</h3>\n')
            
            # Create a modern, responsive table for the flux calibration results
            wlog.write('  <div class="table-responsive">\n')
            wlog.write('    <table class="table sortable">\n')
            wlog.write('      <thead>\n')
            wlog.write('        <tr>\n')
            wlog.write('          <th>Field</th>\n')
            wlog.write('          <th>Spw</th>\n')
            wlog.write('          <th>Flux Density (Jy)</th>\n')
            wlog.write('          <th>Spectral Index</th>\n')
            wlog.write('        </tr>\n')
            wlog.write('      </thead>\n')
            wlog.write('      <tbody>\n')
            
            # Add rows for each calibrator
            for field in calfluxes:
                if field != msinfo['sources']['fluxcal']:
                    field_id = msinfo['sources']['mssources'].split(',').index(field)
                    for spw in calfluxes[field]:
                        flux = calfluxes[field][spw]['fitFluxd']
                        try:
                            spix = calfluxes[field][spw]['spidx'][1]
                            spix_formatted = f"{spix:.2f}"
                        except:
                            spix_formatted = 'N/A'
                        
                        wlog.write('        <tr>\n')
                        wlog.write(f'          <td>{field}</td>\n')
                        wlog.write(f'          <td>{spw}</td>\n')
                        wlog.write(f'          <td>{flux:.3f}</td>\n')
                        wlog.write(f'          <td>{spix_formatted}</td>\n')
                        wlog.write('        </tr>\n')
            
            wlog.write('      </tbody>\n')
            wlog.write('    </table>\n')
            wlog.write('  </div>\n')
            wlog.write('</div>\n')
        except Exception as e:
            logger.warning('Error reading calfluxes file: {}'.format(e))
            wlog.write('<p>Error reading flux calibration information.</p>\n')
    else:
        wlog.write('<p>No flux calibration information available.</p>\n')

def write_applycal_table(wlog, eMCP, ap_dict='applycal_dict'):
    """Display applycal information in a modern format"""
    if ap_dict in eMCP:
        wlog.write('<div class="subsection">\n')
        wlog.write('  <h3>Apply calibration</h3>\n')
        
        wlog.write('  <div class="table-responsive">\n')
        wlog.write('    <table class="table">\n')
        wlog.write('      <thead>\n')
        wlog.write('        <tr>\n')
        wlog.write('          <th>Field</th>\n')
        wlog.write('          <th>G</th>\n')
        wlog.write('          <th>T</th>\n')
        wlog.write('          <th>A</th>\n')
        wlog.write('          <th>F</th>\n')
        wlog.write('          <th>B</th>\n')
        wlog.write('        </tr>\n')
        wlog.write('      </thead>\n')
        wlog.write('      <tbody>\n')
        
        for field in eMCP[ap_dict]:
            gaintable_txt = eMCP[ap_dict][field]['gaintable']
            interp_txt = eMCP[ap_dict][field]['interp']
            
            # Create abbreviations for each calibration type
            cal_types = {
                'G': False,  # Gain calibration
                'T': False,  # Delay calibration
                'A': False,  # Amplitude calibration
                'F': False,  # Bandpass calibration
                'B': False   # Phase calibration
            }
            
            # Check which calibration types are applied
            if any('g.gcal' in s for s in gaintable_txt):
                cal_types['G'] = True
            if any('K0.gcal' in s for s in gaintable_txt):
                cal_types['T'] = True
            if any('Ga.gcal' in s for s in gaintable_txt):
                cal_types['A'] = True
            if any('B0.gcal' in s for s in gaintable_txt):
                cal_types['F'] = True
            if any('B1.gcal' in s for s in gaintable_txt):
                cal_types['B'] = True
            
            # Write row with styled cells
            wlog.write('        <tr>\n')
            wlog.write(f'          <td>{field}</td>\n')
            
            for cal_type, applied in cal_types.items():
                if applied:
                    wlog.write(f'          <td class="cal-applied">✓</td>\n')
                else:
                    wlog.write(f'          <td class="cal-not-applied">✗</td>\n')
            
            wlog.write('        </tr>\n')
        
        wlog.write('      </tbody>\n')
        wlog.write('    </table>\n')
        wlog.write('  </div>\n')
        
        # Legend for the table
        wlog.write('  <div class="cal-legend">\n')
        wlog.write('    <p><strong>Legend:</strong> G = Gain, T = Delay, A = Amplitude, F = Bandpass, B = Phase</p>\n')
        wlog.write('  </div>\n')
        
        wlog.write('</div>\n')

def weblog_calibration(eMCP):
    """Create a focused calibration page showing only required content"""
    weblog_dir = './weblog/'
    
    # Import weblog functions from the modern module
    from .eMCP_weblog_modern import weblog_header, weblog_foot
    
    wlog = open(weblog_dir + "calibration.html", "w")
    
    # Get run name safely
    run_name = eMCP.get('run', eMCP.get('project', 'Pipeline'))
    weblog_header(wlog, 'Calibration', run_name)
    
    # Create calibration tables section
    wlog.write('<div class="section">\n')
    wlog.write('  <h2>Calibration Tables</h2>\n')
    all_calsteps = [
        'bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G0', 'bpcal.BP0', 'allcal_d.K1',
        'allcal_p.G1', 'allcal_ap.G1', 'fluxscale', 'bpcal.BP2', 'allcal_p.G3',
        'allcal_ap.G3', 'phscal_p_scan.G3', 'phscal_ap_scan.G3'
    ]
    if eMCP['is_mixed_mode']:
        all_calsteps.append('narrow_p_offset.G3')
        all_calsteps.append('narrow_bpcal.BP2')
    
    # Check for calibration plots
    all_plots = glob.glob('./weblog/plots/caltables/*png')
    if len(all_plots) > 0:
        wlog.write('  <div class="plot-grid">\n')
        
        # Sort plots by calibration type and antenna
        sorted_plots = sorted(all_plots)
        
        for p in sorted_plots:
            plot_name = os.path.basename(p).replace('.png', '').replace('_', ' ')
            rel_path = 'plots/caltables/' + os.path.basename(p)
            
            wlog.write('    <div class="plot-card">\n')
            wlog.write(f'      <img src="{rel_path}" class="zoomable-img" alt="{plot_name}">\n')
            wlog.write(f'      <div class="plot-caption">{plot_name}</div>\n')
            wlog.write('    </div>\n')
        
        wlog.write('  </div>\n')
    else:
        wlog.write('  <p class="centered">No calibration plots available yet.</p>\n')
    
    wlog.write('</div>\n')
    
    # Flux calibration section - only if available
    if 'msinfo' in eMCP:
        write_fluxscale(wlog, eMCP['msinfo'])
    
    # ApplyCal section - only if available
    write_applycal_table(wlog, eMCP)
    
    weblog_foot(wlog)
    wlog.close()
