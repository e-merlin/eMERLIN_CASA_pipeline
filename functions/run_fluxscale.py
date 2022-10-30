import os
import pickle
import argparse


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vis', dest='vis')
    parser.add_argument('-reference', dest='reference')
    parser.add_argument('-transfer', dest='transfer')
    parser.add_argument('-antenna', dest='antenna')
    parser.add_argument('-caltable', dest='caltable')
    parser.add_argument('-fluxtable', dest='fluxtable')
    parser.add_argument('-listfile', dest='listfile')
    args = parser.parse_args()
    return args


def run_fluxscale(args):
    calfluxes = fluxscale(vis=args.vis,
                          reference=args.reference,
                          transfer=args.transfer,
                          antenna=args.antenna,
                          caltable=args.caltable,
                          fluxtable=args.fluxtable,
                          listfile=args.listfile)
    outfile = os.path.join(weblog_dir, 'calib/calfluxes.pkl')
    save_obj(calfluxes, outfile)


weblog_dir = './weblog/'
args = get_args()
run_fluxscale(args)
