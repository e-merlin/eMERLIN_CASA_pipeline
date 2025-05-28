import os
import numpy as np
import yaml
from casatasks import fluxscale

weblog_dir = './weblog/'


def save_obj(obj, name):
    with open(name, 'w') as f:
        yaml.dump(obj, f, default_flow_style=False)

def convert_ndarray(obj):
    # This is to avoid numpy arrays in the yaml file
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: convert_ndarray(v) for k, v in obj.items()}
    return obj

def run_fluxscale(vis, reference, transfer, antenna, caltable, fluxtable,
                  listfile):
    calfluxes = fluxscale(vis=vis,
                          reference=reference,
                          transfer=transfer,
                          antenna=antenna,
                          caltable=caltable,
                          fluxtable=fluxtable,
                          listfile=listfile)

    outfile = os.path.join(weblog_dir, 'calib/calfluxes.yaml')
    calfluxes_converted = convert_ndarray(calfluxes)
    save_obj(calfluxes_converted, outfile)
