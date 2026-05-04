import os
import numpy as np
from casatasks import fluxscale
from ..utils import eMCP_paths as empaths
from ..utils.eMCP_utils import save_obj

weblog_dir = empaths.WEBLOG_DIR

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
