"""
Flag statistics module for eMERLIN CASA Pipeline.

Collects and saves flag statistics from measurement sets.
"""

import os
import pickle
from casatasks import flagdata

weblog_dir = './weblog/'


def save_obj(obj, name):
    """
    Save a Python object to disk using pickle.
    
    Args:
        obj: Object to save
        name: Output file path
    """
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


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
        weblog_dir, 'plots/plots_flagstats/flagstats_{}.pkl'.format(step))
    save_obj(flag_stats, outfile)
