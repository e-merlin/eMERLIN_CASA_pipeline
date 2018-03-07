import collections

defaults = collections.OrderedDict()

defaults['hanning'] = collections.OrderedDict()
defaults['hanning']['run_hanning']           = True
defaults['hanning']['deloriginal']           = True

defaults['ms2mms'] = collections.OrderedDict()
defaults['ms2mms']['run_ms2mms']             = False
defaults['ms2mms']['mode']                   = 'parallel'

defaults['aoflagger'] = collections.OrderedDict()
defaults['aoflagger']['fields']              = 'all'

defaults['flag_apriori'] = collections.OrderedDict()
defaults['flag_apriori']['do_quack']         = True

defaults['average'] = collections.OrderedDict()
defaults['average']['width']                 = 4

defaults['init_models'] = collections.OrderedDict()
defaults['init_models']['calibrator_models'] = 'calibrator_models/'

defaults['bandpass'] = collections.OrderedDict()
defaults['bandpass']['delay_tablename']      = 'bpcal_d.K0'
defaults['bandpass']['delay_solint']         = '180s'
defaults['bandpass']['delay_combine']        = 'spw'
defaults['bandpass']['delay_prev_cal']       = []
defaults['bandpass']['delay_interp']         = 'linear'
defaults['bandpass']['delay_spw']            = ['*','innerchan']
defaults['bandpass']['delay_minblperant']    = 3
defaults['bandpass']['delay_minsnr']         = 2
defaults['bandpass']['phase_tablename']      = 'bpcal_p.G0'
defaults['bandpass']['phase_solint']         = 'int'
defaults['bandpass']['phase_prev_cal']       = ['bpcal_d.K0']
defaults['bandpass']['phase_interp']         = 'linear'
defaults['bandpass']['phase_combine']        = ''
defaults['bandpass']['phase_spw']            = ['*','innerchan']
defaults['bandpass']['phase_minblperant']    = 3
defaults['bandpass']['phase_minsnr']         = 2
defaults['bandpass']['ap_tablename']         = 'bpcal_ap.G1'
defaults['bandpass']['ap_solint']            = '32s'
defaults['bandpass']['ap_prev_cal']          = ['bpcal_d.K0', 'bpcal_p.G0']
defaults['bandpass']['ap_interp']            = 'linear'
defaults['bandpass']['ap_combine']           = ''
defaults['bandpass']['ap_spw']               = ['*','innerchan']
defaults['bandpass']['ap_minblperant']       = 3
defaults['bandpass']['ap_minsnr']            = 2
defaults['bandpass']['bp_tablename']         = 'bpcal.B0'
defaults['bandpass']['bp_solint']            = 'inf'
defaults['bandpass']['bp_combine']           = 'field,scan'
defaults['bandpass']['bp_interp']            = 'nearest,linear'
defaults['bandpass']['bp_spw']               = ['*','']
defaults['bandpass']['bp_uvrange']           = ''
defaults['bandpass']['bp_solnorm']           = True
defaults['bandpass']['bp_prev_cal']          = ['bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G1']
defaults['bandpass']['apply_calibrators']    = ['bpcal.B0']
defaults['bandpass']['apply_targets']        = ['bpcal.B0']

defaults['flag_tfcrop'] = collections.OrderedDict()
defaults['flag_tfcrop']['mode']              = 'tfcrop'
defaults['flag_tfcrop']['antenna']           = ''
defaults['flag_tfcrop']['scan']              = ''
defaults['flag_tfcrop']['spw']               = ''
defaults['flag_tfcrop']['correlation']       = 'ABS_ALL'
defaults['flag_tfcrop']['ntime']             = '90min'
defaults['flag_tfcrop']['combinescans']      = True
defaults['flag_tfcrop']['datacolumn']        = 'corrected'
defaults['flag_tfcrop']['winsize']           = 3
defaults['flag_tfcrop']['timecutoff']        = 4.0
defaults['flag_tfcrop']['freqcutoff']        = 3.0
defaults['flag_tfcrop']['maxnpieces']        = 1
defaults['flag_tfcrop']['usewindowstats']    = 'sum'
defaults['flag_tfcrop']['halfwin']           = 3
defaults['flag_tfcrop']['extendflags']       = True
defaults['flag_tfcrop']['action']            = 'apply'
defaults['flag_tfcrop']['display']           = ''
defaults['flag_tfcrop']['flagbackup']        = True

defaults['delay'] = collections.OrderedDict()
defaults['delay']['use_fringefit']           = False
defaults['delay']['tablename']               = 'delay.K1'
defaults['delay']['solint']                  = '300s'
defaults['delay']['combine']                 = 'spw'
defaults['delay']['prev_cal']                = ['bpcal.B0']
defaults['delay']['interp']                  = 'linear'
defaults['delay']['spw']                     = ['*','innerchan']
defaults['delay']['zerorates']               = True
defaults['delay']['minblperant']             = 3
defaults['delay']['minsnr']                  = 3
defaults['delay']['apply_calibrators']       = ['bpcal.B0','delay.K1']
defaults['delay']['apply_targets']           = ['bpcal.B0','delay.K1']

defaults['gain_p_ap'] = collections.OrderedDict()
defaults['gain_p_ap']['p_tablename']         = 'allcal_p.G0'
defaults['gain_p_ap']['p_prev_cal']          = ['bpcal.B0','delay.K1']
defaults['gain_p_ap']['p_solint']            = '16s'
defaults['gain_p_ap']['p_spw']               = ['*','innerchan']
defaults['gain_p_ap']['p_combine']           = ''
defaults['gain_p_ap']['p_interp']            = 'linear'
defaults['gain_p_ap']['p_minblperant']       = 3
defaults['gain_p_ap']['p_minsnr']            = 3
defaults['gain_p_ap']['ji_tablename']        = 'allcal_p_jitter.G0'
defaults['gain_p_ap']['ji_prev_cal']         = ['bpcal.B0','delay.K1','allcal_p.G0']
defaults['gain_p_ap']['ji_solint']           = 'int'
defaults['gain_p_ap']['ji_spw']              = ['*','innerchan']
defaults['gain_p_ap']['ji_combine']          = 'spw'
defaults['gain_p_ap']['ji_interp']           = 'linear'
defaults['gain_p_ap']['ji_minblperant']      = 3
defaults['gain_p_ap']['ji_minsnr']           = 3
defaults['gain_p_ap']['ap_tablename']        = 'allcal_ap.G1'
defaults['gain_p_ap']['ap_prev_cal']         = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_p_jitter.G0']
defaults['gain_p_ap']['ap_solint']           = '32s'
defaults['gain_p_ap']['ap_spw']              = ['*','innerchan']
defaults['gain_p_ap']['ap_combine']          = ''
defaults['gain_p_ap']['ap_interp']           = 'linear'
defaults['gain_p_ap']['ap_minblperant']      = 3
defaults['gain_p_ap']['ap_minsnr']           = 3
defaults['gain_p_ap']['apply_calibrators']   = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1']
defaults['gain_p_ap']['apply_targets']       = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_ap.G1']

defaults['flux'] = collections.OrderedDict()
defaults['flux']['tablename']           = 'allcall_ap.G1_fluxscaled'
defaults['flux']['ampcal_table']        = 'allcal_ap.G1'
defaults['flux']['apply_calibrators']   = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled']
defaults['flux']['apply_targets']       = ['bpcal.B0','delay.K1','allcal_p.G0','allcal_ap.G1_fluxscaled']


defaults['flag_rflag'] = collections.OrderedDict()
defaults['flag_rflag']['mode']               = 'rflag'
defaults['flag_rflag']['antenna']            = ''
defaults['flag_rflag']['scan']               = ''
defaults['flag_rflag']['spw']                = ''
defaults['flag_rflag']['correlation']        = ''
defaults['flag_rflag']['ntime']              = '90min'
defaults['flag_rflag']['combinescans']       = True
defaults['flag_rflag']['datacolumn']         = 'corrected'
defaults['flag_rflag']['timedevscale']       = 5.0
defaults['flag_rflag']['freqdevscale']       = 5.0
defaults['flag_rflag']['action']             = 'apply'
defaults['flag_rflag']['display']            = ''
defaults['flag_rflag']['flagbackup']         = True

