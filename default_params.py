import collections

defaults = collections.OrderedDict()

defaults['global'] = collections.OrderedDict()
defaults['global']['refantmode']                      = 'flex'

defaults['import_eM'] = collections.OrderedDict()
defaults['import_eM']['constobsid']                   = True
defaults['import_eM']['scanreindexgap_s']             = 15.0
defaults['import_eM']['antenna']                      = ''
defaults['import_eM']['timeaverage']                  = True
defaults['import_eM']['timebin']                      = '4s'
defaults['import_eM']['chanaverage']                  = False
defaults['import_eM']['chanbin']                      = 1
defaults['import_eM']['usewtspectrum']                = False
defaults['import_eM']['run_hanning']                  = 'auto'
defaults['import_eM']['ms2mms']                       = False

defaults['aoflagger'] = collections.OrderedDict()
defaults['aoflagger']['run']                          = 'auto'
defaults['aoflagger']['fields']                       = 'all'
defaults['aoflagger']['separate_bands']               = False

defaults['flag_apriori'] = collections.OrderedDict()
defaults['flag_apriori']['do_quack']                  = True
defaults['flag_apriori']['Lo_dropout']                = ''
defaults['flag_apriori']['Lo_datacolumn']             = 'data'
defaults['flag_apriori']['Lo_useflags']               = True
defaults['flag_apriori']['Lo_spws']                   = ['3']
defaults['flag_apriori']['Lo_threshold']              = 0.5
defaults['flag_apriori']['Lo_min_scans']              = 4

defaults['average'] = collections.OrderedDict()
defaults['average']['field']                          = ''
defaults['average']['chanbin']                        = 4
defaults['average']['datacolumn']                     = 'data'
defaults['average']['timerange']                      = ''
defaults['average']['scan']                           = ''
defaults['average']['antenna']                        = ''
defaults['average']['shift_phasecenter']              = False

defaults['plot_data'] = collections.OrderedDict()
defaults['plot_data']['num_proc']                     = 1

defaults['init_models'] = collections.OrderedDict()
defaults['init_models']['calibrator_models']          = 'calibrator_models/'

defaults['bandpass'] = collections.OrderedDict()
defaults['bandpass']['delay_tablename']               = 'bpcal_d.K0'
defaults['bandpass']['delay_solint']                  = '180s'
defaults['bandpass']['delay_combine']                 = 'spw'
defaults['bandpass']['delay_prev_cal']                = []
defaults['bandpass']['delay_interp']                  = 'linear'
defaults['bandpass']['delay_spw']                     = ['*','innerchan']
defaults['bandpass']['delay_minblperant']             = 3
defaults['bandpass']['delay_minsnr']                  = 2
defaults['bandpass']['phase_tablename']               = 'bpcal_p.G0'
defaults['bandpass']['phase_solint']                  = 'int'
defaults['bandpass']['phase_prev_cal']                = ['bpcal_d.K0']
defaults['bandpass']['phase_interp']                  = 'linear'
defaults['bandpass']['phase_combine']                 = ''
defaults['bandpass']['phase_spw']                     = ['*','innerchan']
defaults['bandpass']['phase_minblperant']             = 3
defaults['bandpass']['phase_minsnr']                  = 2
defaults['bandpass']['ap_tablename']                  = 'bpcal_ap.G1'
defaults['bandpass']['ap_solint']                     = '32s'
defaults['bandpass']['ap_prev_cal']                   = ['bpcal_d.K0', 'bpcal_p.G0']
defaults['bandpass']['ap_interp']                     = 'linear'
defaults['bandpass']['ap_combine']                    = ''
defaults['bandpass']['ap_spw']                        = ['*','innerchan']
defaults['bandpass']['ap_minblperant']                = 3
defaults['bandpass']['ap_minsnr']                     = 2
defaults['bandpass']['bp_tablename']                  = 'bpcal.B0'
defaults['bandpass']['bp_solint']                     = 'inf'
defaults['bandpass']['bp_combine']                    = 'field,scan'
defaults['bandpass']['bp_interp']                     = 'nearest,cubicflag'
defaults['bandpass']['bp_spw']                        = ['*','']
defaults['bandpass']['bp_uvrange']                    = ''
defaults['bandpass']['bp_fillgaps']                   = 8
defaults['bandpass']['bp_solnorm']                    = True
defaults['bandpass']['bp_prev_cal']                   = ['bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G1']
defaults['bandpass']['apply_calibrators']             = ['bpcal.B0']
defaults['bandpass']['apply_targets']                 = []
defaults['bandpass']['run_flag']                      = True
defaults['bandpass']['tfcrop'] = collections.OrderedDict()
defaults['bandpass']['tfcrop']['mode']                = 'tfcrop'
defaults['bandpass']['tfcrop']['sources']             = 'bpcal'
defaults['bandpass']['tfcrop']['antenna']             = ''
defaults['bandpass']['tfcrop']['scan']                = ''
defaults['bandpass']['tfcrop']['spw']                 = ''
defaults['bandpass']['tfcrop']['correlation']         = ''
defaults['bandpass']['tfcrop']['ntime']               = ''
defaults['bandpass']['tfcrop']['combinescans']        = False
defaults['bandpass']['tfcrop']['datacolumn']          = 'corrected'
defaults['bandpass']['tfcrop']['winsize']             = 3
defaults['bandpass']['tfcrop']['timecutoff']          = 4.5
defaults['bandpass']['tfcrop']['freqcutoff']          = 4.5
defaults['bandpass']['tfcrop']['maxnpieces']          = 7
defaults['bandpass']['tfcrop']['uwstats']             = 'none'
defaults['bandpass']['tfcrop']['halfwin']             = 1
defaults['bandpass']['tfcrop']['extendflags']         = True
defaults['bandpass']['tfcrop']['action']              = 'apply'
defaults['bandpass']['tfcrop']['display']             = ''
defaults['bandpass']['tfcrop']['flagbackup']          = True
defaults['bp_apply_mid'] = collections.OrderedDict()
defaults['bp_apply_mid']['apply_calibrators']         = ['bpcal_d.K0','bpcal_p.G0','bpcal_ap.G1','bpcal.B0']
defaults['bp_apply_mid']['apply_targets']             = []

defaults['initial_gaincal'] = collections.OrderedDict()
defaults['initial_gaincal']['delay'] = collections.OrderedDict()
defaults['initial_gaincal']['delay']['use_fringefit'] = False
defaults['initial_gaincal']['delay']['tablename']     = 'delay.K1'
defaults['initial_gaincal']['delay']['solint']        = '180s'
defaults['initial_gaincal']['delay']['combine']       = 'spw'
defaults['initial_gaincal']['delay']['prev_cal']      = ['bpcal.B0']
defaults['initial_gaincal']['delay']['interp']        = 'linear'
defaults['initial_gaincal']['delay']['spw']           = ['*','innerchan']
defaults['initial_gaincal']['delay']['zerorates']     = True
defaults['initial_gaincal']['delay']['minblperant']   = 3
defaults['initial_gaincal']['delay']['minsnr']        = 2
defaults['initial_gaincal']['p_tablename']            = 'allcal_p.G0'
defaults['initial_gaincal']['p_prev_cal']             = ['bpcal.B0','delay.K1']
defaults['initial_gaincal']['p_solint']               = 'int'
defaults['initial_gaincal']['p_spw']                  = ['*','innerchan']
defaults['initial_gaincal']['p_combine']              = ''
defaults['initial_gaincal']['p_interp']               = 'linear'
defaults['initial_gaincal']['p_minblperant']          = 3
defaults['initial_gaincal']['p_minsnr']               = 2
defaults['initial_gaincal']['ap_tablename']           = 'allcal_ap.G1'
defaults['initial_gaincal']['ap_prev_cal']            = ['bpcal.B0','delay.K1','allcal_p.G0']
defaults['initial_gaincal']['ap_solint']              = '32s'
defaults['initial_gaincal']['ap_spw']                 = ['*','innerchan']
defaults['initial_gaincal']['ap_combine']             = ''
defaults['initial_gaincal']['ap_interp']              = 'linear'
defaults['initial_gaincal']['ap_minblperant']         = 3
defaults['initial_gaincal']['ap_minsnr']              = 2
defaults['initial_gaincal']['apply_calibrators']      = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_ap.G1']
defaults['initial_gaincal']['apply_targets']          = []
defaults['initial_gaincal']['flagmode']               = 'tfcrop'
defaults['initial_gaincal']['tfcrop'] = collections.OrderedDict()
defaults['initial_gaincal']['tfcrop']['mode']         = 'tfcrop'
defaults['initial_gaincal']['tfcrop']['sources']      = 'calsources'
defaults['initial_gaincal']['tfcrop']['antenna']      = ''
defaults['initial_gaincal']['tfcrop']['scan']         = ''
defaults['initial_gaincal']['tfcrop']['spw']          = ''
defaults['initial_gaincal']['tfcrop']['correlation']  = ''
defaults['initial_gaincal']['tfcrop']['ntime']        = ''
defaults['initial_gaincal']['tfcrop']['combinescans'] = False
defaults['initial_gaincal']['tfcrop']['datacolumn']   = 'residual'
defaults['initial_gaincal']['tfcrop']['winsize']      = 3
defaults['initial_gaincal']['tfcrop']['timecutoff']   = 5.0
defaults['initial_gaincal']['tfcrop']['freqcutoff']   = 5.0
defaults['initial_gaincal']['tfcrop']['maxnpieces']   = 7
defaults['initial_gaincal']['tfcrop']['uwstats']      = 'none'
defaults['initial_gaincal']['tfcrop']['halfwin']      = 1
defaults['initial_gaincal']['tfcrop']['extendflags']  = True
defaults['initial_gaincal']['tfcrop']['action']       = 'apply'
defaults['initial_gaincal']['tfcrop']['display']      = ''
defaults['initial_gaincal']['tfcrop']['flagbackup']   = True
defaults['initial_gaincal']['rflag'] = collections.OrderedDict()
defaults['initial_gaincal']['rflag']['mode']          = 'rflag'
defaults['initial_gaincal']['rflag']['sources']       = 'calsources'
defaults['initial_gaincal']['rflag']['antenna']       = ''
defaults['initial_gaincal']['rflag']['scan']          = ''
defaults['initial_gaincal']['rflag']['spw']           = ''
defaults['initial_gaincal']['rflag']['correlation']   = ''
defaults['initial_gaincal']['rflag']['ntime']         = ''
defaults['initial_gaincal']['rflag']['combinescans']  = False
defaults['initial_gaincal']['rflag']['datacolumn']    = 'residual'
defaults['initial_gaincal']['rflag']['timedevscale']  = 4.5
defaults['initial_gaincal']['rflag']['freqdevscale']  = 4.5
defaults['initial_gaincal']['rflag']['extendflags']   = True
defaults['initial_gaincal']['rflag']['action']        = 'apply'
defaults['initial_gaincal']['rflag']['display']       = ''
defaults['initial_gaincal']['rflag']['flagbackup']    = True

defaults['fluxscale'] = collections.OrderedDict()
defaults['fluxscale']['tablename']                    = 'allcal_ap.G1_fluxscaled'
defaults['fluxscale']['ampcal_table']                 = 'allcal_ap.G1'
defaults['fluxscale']['apply_calibrators']            = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled']
defaults['fluxscale']['apply_targets']                = []

defaults['bandpass_sp'] = collections.OrderedDict()
defaults['bandpass_sp']['bp_tablename']               = 'bpcal_sp.B1'
defaults['bandpass_sp']['bp_prev_cal']                = ['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled']
defaults['bandpass_sp']['bp_solint']                  = 'inf'
defaults['bandpass_sp']['bp_spw']                     = ['*','']
defaults['bandpass_sp']['bp_combine']                 = 'field,scan'
defaults['bandpass_sp']['bp_interp']                  = 'nearest,cubicflag'
defaults['bandpass_sp']['bp_uvrange']                 = ''
defaults['bandpass_sp']['bp_fillgaps']                = 8
defaults['bandpass_sp']['bp_solnorm']                 = False
defaults['bandpass_sp']['apply_calibrators']          = ['delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1']
defaults['bandpass_sp']['apply_targets']              = []

defaults['gain_amp_sp'] = collections.OrderedDict()
defaults['gain_amp_sp']['ap_tablename']               = 'allcal_ap.G3'
defaults['gain_amp_sp']['ap_prev_cal']                = ['delay.K1','allcal_p.G0','bpcal_sp.B1']
defaults['gain_amp_sp']['ap_solint']                  = '32s'
defaults['gain_amp_sp']['ap_spw']                     = ['*','innerchan']
defaults['gain_amp_sp']['ap_combine']                 = ''
defaults['gain_amp_sp']['ap_interp']                  = 'linear'
defaults['gain_amp_sp']['ap_minblperant']             = 3
defaults['gain_amp_sp']['ap_minsnr']                  = 2
defaults['gain_amp_sp']['p_scan_tablename']           = 'phscal_p_scan.G3'
defaults['gain_amp_sp']['p_scan_prev_cal']            = ['delay.K1','bpcal_sp.B1']
defaults['gain_amp_sp']['p_scan_solint']              = 'inf'
defaults['gain_amp_sp']['p_scan_spw']                 = ['*','innerchan']
defaults['gain_amp_sp']['p_scan_combine']             = ''
defaults['gain_amp_sp']['p_scan_interp']              = 'linear'
defaults['gain_amp_sp']['p_scan_minblperant']         = 3
defaults['gain_amp_sp']['p_scan_minsnr']              = 2
defaults['gain_amp_sp']['ap_scan_tablename']          = 'phscal_ap_scan.G3'
defaults['gain_amp_sp']['ap_scan_prev_cal']           = ['delay.K1','bpcal_sp.B1','allcal_p.G0']
defaults['gain_amp_sp']['ap_scan_solint']             = 'inf'
defaults['gain_amp_sp']['ap_scan_spw']                = ['*','innerchan']
defaults['gain_amp_sp']['ap_scan_combine']            = ''
defaults['gain_amp_sp']['ap_scan_interp']             = 'linear'
defaults['gain_amp_sp']['ap_scan_minblperant']        = 3
defaults['gain_amp_sp']['ap_scan_minsnr']             = 2
defaults['gain_amp_sp']['apply_calibrators']          = ['delay.K1','allcal_p.G0','bpcal_sp.B1','allcal_ap.G3']
defaults['gain_amp_sp']['apply_targets']              = []

defaults['applycal_all'] = collections.OrderedDict()
defaults['applycal_all']['apply_calibrators']         = ['delay.K1','allcal_p.G0','bpcal_sp.B1','allcal_ap.G3']
defaults['applycal_all']['apply_targets']             = ['delay.K1','bpcal_sp.B1','phscal_p_scan.G3','phscal_ap_scan.G3']

defaults['flag_target'] = collections.OrderedDict()
defaults['flag_target']['mode_to_run']                = 'tfcrop'
defaults['flag_target']['rflag'] = collections.OrderedDict()
defaults['flag_target']['rflag']['mode']              = 'rflag'
defaults['flag_target']['rflag']['sources']           = 'targets'
defaults['flag_target']['rflag']['antenna']           = ''
defaults['flag_target']['rflag']['scan']              = ''
defaults['flag_target']['rflag']['spw']               = ''
defaults['flag_target']['rflag']['correlation']       = ''
defaults['flag_target']['rflag']['ntime']             = 'scan'
defaults['flag_target']['rflag']['combinescans']      = False
defaults['flag_target']['rflag']['datacolumn']        = 'corrected'
defaults['flag_target']['rflag']['timedevscale']      = 4.5
defaults['flag_target']['rflag']['freqdevscale']      = 4.5
defaults['flag_target']['rflag']['extendflags']       = True
defaults['flag_target']['rflag']['action']            = 'apply'
defaults['flag_target']['rflag']['display']           = ''
defaults['flag_target']['rflag']['flagbackup']        = True
defaults['flag_target']['tfcrop'] = collections.OrderedDict()
defaults['flag_target']['tfcrop']['mode']             = 'tfcrop'
defaults['flag_target']['tfcrop']['sources']          = 'targets'
defaults['flag_target']['tfcrop']['antenna']          = ''
defaults['flag_target']['tfcrop']['scan']             = ''
defaults['flag_target']['tfcrop']['spw']              = ''
defaults['flag_target']['tfcrop']['correlation']      = ''
defaults['flag_target']['tfcrop']['ntime']            = ''
defaults['flag_target']['tfcrop']['combinescans']     = False
defaults['flag_target']['tfcrop']['datacolumn']       = 'corrected'
defaults['flag_target']['tfcrop']['winsize']          = 3
defaults['flag_target']['tfcrop']['timecutoff']       = 4.5
defaults['flag_target']['tfcrop']['freqcutoff']       = 4.5
defaults['flag_target']['tfcrop']['maxnpieces']       = 7
defaults['flag_target']['tfcrop']['uwstats']          = 'none'
defaults['flag_target']['tfcrop']['halfwin']          = 1
defaults['flag_target']['tfcrop']['extendflags']      = True
defaults['flag_target']['tfcrop']['action']           = 'apply'
defaults['flag_target']['tfcrop']['display']          = ''
defaults['flag_target']['tfcrop']['flagbackup']       = True


defaults['first_images'] = collections.OrderedDict()
defaults['first_images']['run_statwt']                = True
defaults['first_images']['statwt_timebin']            = '32s'
defaults['first_images']['imsize']                    = 1024
defaults['first_images']['niter']                     = 80
defaults['first_images']['deconvolver']               = 'hogbom'
defaults['first_images']['nterms']                    = 1
#defaults['first_images']['deconvolver']              = 'mtmfs'
#defaults['first_images']['nterms']                   = 2
defaults['first_images']['weighting']                 = 'briggs'
defaults['first_images']['robust']                    = 0.5
defaults['first_images']['gain']                      = 0.1
defaults['first_images']['nsigma']                    = 5.0
defaults['first_images']['sidelobethreshold']         = 1.0
defaults['first_images']['noisethreshold']            = 8.0
defaults['first_images']['lownoisethreshold']         = 1.5
defaults['first_images']['minbeamfrac']               = 0.2
defaults['first_images']['growiterations']            = 25
defaults['first_images']['level0']                    = 3.
