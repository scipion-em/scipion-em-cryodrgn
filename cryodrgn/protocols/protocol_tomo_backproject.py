# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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
from enum import Enum
from emtable import Table

from pyworkflow.constants import BETA
from pyworkflow.utils import Message, createAbsLink, findRootFrom
import pyworkflow.protocol.params as params
from pyworkflow.object import Float
from pwem.objects import SetOfFSCs, FSC
from tomo.objects import AverageSubTomogram

from cryodrgn import Plugin
from cryodrgn.constants import FSC_COLUMNS, V3_4_0
from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class outputs(Enum):
    subtomogramAverage = AverageSubTomogram
    FSCs = SetOfFSCs


class CryoDrgnProtBackProject(CryoDrgnProtBase):
    """ Backprojection of tilt-series particles from Warp/M.

    Can be used to confirm our inputs were correctly
    parsed using traditional homogeneous reconstruction.
    """
    _label = 'reconstruct tilt-series particles'
    _devStatus = BETA
    _possibleOutputs = outputs

    @classmethod
    def isDisabled(cls):
        return not Plugin.versionGE(V3_4_0)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.pixelSize = Float()
        self.imgPath = None  # no need to save between steps

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addHidden(params.GPU_LIST, params.StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")
        form.addParam('starFile', params.FileParam, filter="*.star",
                      label='Star file from Warp',
                      help="cryoDRGN-ET expects 2D particle tilt series "
                           "images in a .star file exported from Windows "
                           "Warp/M. 2D particle tilt series images without "
                           "CTF premultiplication can be extracted in RELION5 "
                           "using the --no_ctf option however, RELION5 and "
                           "Linux Warp/M star files are not currently supported.")
        form.addParam('dosePerTilt', params.FloatParam,
                      default=3.0,
                      label='Dose per tilt image')
        form.addParam('ntilts', params.IntParam,
                      default=10,
                      label="Number of tilts per particle",
                      help="Number of tilts per particle to backproject "
                           "(default: 10)")

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.reconstructStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Create the input star, poses and ctf pkl files as expected by cryoDRGN. """
        inputFn = self.starFile.get().strip()
        starFn = self._getFileName('input_parts')
        if os.path.exists(inputFn):
            createAbsLink(os.path.abspath(inputFn), starFn)

            # find pixel size and box size
            row = Table(fileName=starFn)[0]
            self.pixelSize.set(row.rlnDetectorPixelSize * 10000 / row.rlnMagnification)
            imgFn = row.rlnImageName.split("@")[-1]
            self.imgPath = findRootFrom(inputFn, imgFn)
            boxSize = self._getBoxSize(os.path.join(self.imgPath, imgFn))

            extraArgs = [
                f"-D {boxSize}",
                f"--Apix {self.pixelSize.get()}",
            ]
            self._runProgram('parse_pose_star',
                             self._getParsePosesArgs() + extraArgs)
            self._runProgram('parse_ctf_star',
                             self._getParseCtfArgs() + extraArgs)

    def reconstructStep(self):
        self._runProgram('backproject_voxel',
                         self._getBackprjArgs())

    def createOutputStep(self):
        volume = AverageSubTomogram()
        volume.setFileName(self._getFileName('subtomo_avg'))
        volume.setHalfMaps([self._getFileName('subtomo_avg_half1'),
                            self._getFileName('subtomo_avg_half2')])
        volume.setSamplingRate(self.pixelSize.get())

        # Output FSC
        setOfFSC = self.genFSCs(self._getFileName('subtomo_avg_fsc'),
                                FSC_COLUMNS)

        self._defineOutputs(**{outputs.subtomogramAverage.name: volume,
                               outputs.FSCs.name: setOfFSC})

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []

        return errors

    def _citations(self):
        return ['Rangan2023']

    # --------------------------- UTILS functions -----------------------------
    def _getBackprjArgs(self):
        args = [
            self._getFileName('input_parts'),
            f"--poses {self._getFileName('input_poses')}",
            f"--ctf {self._getFileName('input_ctfs')}",
            "--tilt",
            f"--ntilts {self.ntilts}",
            f"--dose-per-tilt {self.dosePerTilt.get()}",
            f"-o {self.getOutputDir()}",
            f"--datadir {self.imgPath}"
        ]

        return args

    def genFSCs(self, fscFile, fscColumns):
        fscSet = self._createSetOfFSCs()
        data = np.loadtxt(fscFile, dtype=float, skiprows=1)
        resolution_inv = (data[:, 0] / self.pixelSize.get()).tolist()
        for columnIndex, columnName in enumerate(fscColumns):
            columnValues = data[:, columnIndex+1].tolist()
            fsc = FSC(objLabel=columnName)
            fsc.setData(resolution_inv, columnValues)
            fscSet.append(fsc)

        fscSet.write()
        return fscSet

    @staticmethod
    def _getBoxSize(fn):
        from pwem.emlib.image import ImageHandler
        ih = ImageHandler()
        x, *_ = ih.getDimensions(fn)
        return x
