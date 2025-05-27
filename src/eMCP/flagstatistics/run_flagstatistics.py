"""
Flag statistics module for eMERLIN CASA Pipeline.

Collects and saves flag statistics from measurement sets.
"""

import os
from casatasks import flagdata
from ..utils.eMCP_utils import save_obj

weblog_dir = './weblog/'


def run_flagstats(msfile, step):
    """
    Calculate flag statistics and save results.
    
    Args:
        msfile: Path to measurement set
        step: Pipeline step identifier for filename
    """
    flag_stats = flagdata(vis=msfile,
                          mode='summary',
                          action='calculate',
                          display='none',
                          antenna='*&*',
                          flagbackup=False)
    outfile = os.path.join(
        weblog_dir, 'plots/plots_flagstats/flagstats_{}.yaml'.format(step))
    save_obj(flag_stats, outfile)
