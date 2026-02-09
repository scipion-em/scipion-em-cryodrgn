# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Eduardo García (eduardo.garcia@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] Unidad de  Biocomputacion, Centro Nacional de Biotecnologia, CSIC (CNB-CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import numpy as np

from pyworkflow.utils.process import runJob
import cryodrgn
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
            f" --server_functions_path {server_functions_path} --env_name {cryodrgn.Plugin.getCryoDrgnEnvActivation().split(' ')[-1]}")

    if gpu_id is not None:
        args += f" --gpu_id {gpu_id}"

    return args