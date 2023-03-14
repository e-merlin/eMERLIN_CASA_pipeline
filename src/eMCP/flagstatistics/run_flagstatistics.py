import os
import pickle
from casatasks import flagdata

weblog_dir = './weblog/'


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def run_flagstats(msfile, step):
    flag_stats = flagdata(vis=msfile,
                          mode='summary',
                          action='calculate',
                          display='none',
                          antenna='*&*',
                          flagbackup=False)
    outfile = os.path.join(
        weblog_dir, 'plots/plots_flagstats/flagstats_{}.pkl'.format(step))
    save_obj(flag_stats, outfile)
