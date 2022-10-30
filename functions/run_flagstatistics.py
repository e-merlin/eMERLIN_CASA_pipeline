import os
import pickle
import argparse


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-msfile', dest='msfile', help='msfile path')
    parser.add_argument('-step', dest='step', help='step')
    args = parser.parse_args()
    return args


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


weblog_dir = './weblog/'
args = get_args()
run_flagstats(args.msfile, args.step)
