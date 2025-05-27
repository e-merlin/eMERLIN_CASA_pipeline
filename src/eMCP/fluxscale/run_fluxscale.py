import os
import yaml
from casatasks import fluxscale

weblog_dir = './weblog/'


def save_obj(obj, name):
    with open(name, 'w') as f:
        yaml.dump(obj, f, default_flow_style=False)


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
    save_obj(calfluxes, outfile)
