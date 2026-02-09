import os
import numpy as np

from pyworkflow.utils.process import runJob
from cryodrgn import Plugin
from cryodrgn.constants import *


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

def getAnnotateSpaceArguments(particles, gpu_id=None):
    server_functions_path = os.path.join(os.path.dirname(cryodrgn.__file__), "utils", "annotate_space_server.py")
    args = (f"--config {particles.getFlexInfo().getAttr(CONFIG)} --load {particles.getFlexInfo().getAttr(WEIGHTS)}"
            f"--server_functions_path {server_functions_path} --env_name {opusdsd.Plugin.getCryoDrgnEnvActivation().split(' ')[-1]}")

    if gpu_id is not None:
        args += f" --gpu_id {gpu_id}"

    return args