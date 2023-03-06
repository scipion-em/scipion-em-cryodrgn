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

import pyworkflow.utils as pwutils
from pyworkflow.constants import NEW
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj

from .. import Plugin
from ..constants import *
from .protocol_base import CryoDrgnProtBase


class CryoDrgnProtAbinitio(CryoDrgnProtBase):
    """
    Protocol to run ab-initio reconstruction with cryoDRGN2 neural network.
    """
    _label = 'training ab initio'
    _devStatus = NEW

    @classmethod
    def isDisabled(cls):
        return not Plugin.versionGE(V2_1_0)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        CryoDrgnProtBase._defineParams(self, form)
        form.getParam('zDim').default = pwobj.Integer(1)
        form.getParam('numEpochs').default = pwobj.Integer(30)

    def _defineAdvancedParams(self, form):
        form.addParam('protType', params.EnumParam,
                      choices=['homogeneous', 'heterogeneous'],
                      default=AB_INITIO_HETERO,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Ab initio type')

        form.addSection(label='Advanced')
        form.addParam('searchRange', params.IntParam, default=10,
                      label='Translational search range (px)')
        form.addParam('psFreq', params.IntParam, default=5,
                      label='Update poses every N epochs')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn abinit_homo -h* for help.")

    # --------------------------- STEPS functions -----------------------------
    def runTrainingStep(self):
        # Create output folder
        pwutils.cleanPath(self.getOutputDir())
        pwutils.makePath(self.getOutputDir())

        # Call cryoDRGN with the appropriate parameters
        protType = self.protType.get()
        program = "homo" if protType == AB_INITIO_HOMO else "het"
        self._runProgram("abinit_" + program, self._getTrainingArgs(protType))

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["Training ab initio for %d epochs." % self.numEpochs]

        return summary

    def _validate(self):
        errors = CryoDrgnProtBase._validateBase(self)

        if self.zDim > 1 and self.protType.get() == AB_INITIO_HOMO:
            errors.append("Latent variable must be 1 for "
                          "homogeneous reconstruction")
        if self.zDim == 1 and self.protType.get() == AB_INITIO_HETERO:
            errors.append("Latent variable must be >1 for "
                          "heterogeneous reconstruction")

        return errors

    def _citations(self):
        return ['Zhong2021b']

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self, protType=AB_INITIO_HOMO):
        parts = self.inputParticles.get()

        args = [
            parts.filename.get(),
            '--ctf %s' % parts.ctfs,
            '-o %s ' % self.getOutputDir(),
            '-n %d' % self.numEpochs,
            '--t-extent %d' % self.searchRange,
            '--ps-freq %d' % self.psFreq,
        ]

        if protType == AB_INITIO_HETERO:
            args.extend([
                '--zdim %d' % self.zDim,
                '--preprocessed',
                '--max-threads %d' % self.numberOfThreads,
            ])

            if len(self.getGpuList()) > 1:
                args.append('--multigpu')

        if self.extraParams.hasValue():
            args.append(self.extraParams.get())

        return args
