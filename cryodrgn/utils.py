import os
import numpy as np

from pyworkflow.utils.process import runJob
from cryodrgn import Plugin


def generateVolumes(zValues, weights, config, outdir, apix=1, flip=False,
                    downsample=None, invert=False):
    """
    Function to call cryodrgn eval_vol and generate new volumes
    """
    program = 'eval_vol'
    args = _getEvalVolArgs(zValues, weights, config, outdir, apix, flip,
                           downsample, invert)
    runJob(None, Plugin.getProgram(program, gpus='0'), ' '.join(args),
           env=Plugin.getEnviron())


def _getEvalVolArgs(zvalues, weights, config, outdir, apix, flip,
                    downsample, invert):
    os.makedirs(outdir, exist_ok=True)
    np.savetxt(f'{outdir}/zfile.txt', zvalues)
    zfilePath = os.path.abspath(os.path.join(outdir, 'zfile.txt'))

    return [
        weights,
        f"--config {config}",
        f"--zfile {zfilePath}",
        f"-o {outdir}",
        f"--Apix {apix}",
        "--flip" if flip else "",
        f"-d {downsample}" if downsample is not None else "",
        "--invert" if invert else ""
    ]
