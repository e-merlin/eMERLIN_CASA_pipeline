import os
import pickle
from casatasks import fluxscale

weblog_dir = './weblog/'


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)


def run_fluxscale(vis, reference, transfer, antenna, caltable, fluxtable,
                  listfile):
    calfluxes = fluxscale(vis=vis,
                          reference=reference,
                          transfer=transfer,
                          antenna=antenna,
                          caltable=caltable,
                          fluxtable=fluxtable,
                          listfile=listfile)
    outfile = os.path.join(weblog_dir, 'calib/calfluxes.pkl')
    save_obj(calfluxes, outfile)
