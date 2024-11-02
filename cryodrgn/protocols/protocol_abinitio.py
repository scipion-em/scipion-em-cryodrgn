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
from enum import Enum

from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
from pwem.objects import SetOfParticlesFlex, Volume

from cryodrgn.constants import AB_INITIO_HOMO, AB_INITIO_HETERO
from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class outputs(Enum):
    Particles = SetOfParticlesFlex
    Volumes = Volume


class CryoDrgnProtAbinitio(CryoDrgnProtBase):
    """
    Protocol to run ab-initio reconstruction with cryoDRGN neural network.
    """
    _label = 'training ab initio'
    _devStatus = PROD
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        CryoDrgnProtBase._defineParams(self, form)
        form.getParam('numEpochs').default = pwobj.Integer(30)
        form.getParam('zDim').default = pwobj.Integer(8)

    def _defineAdvancedParams(self, form):
        form.addParam('protType', params.EnumParam,
                      condition='not doContinue',
                      choices=['homogeneous', 'heterogeneous'],
                      default=AB_INITIO_HETERO,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Ab initio type')

        form.addSection(label='Advanced')
        form.addParam('searchRange', params.IntParam, default=10,
                      condition='not doContinue',
                      label='Translational search range (px)')
        form.addParam('psFreq', params.IntParam, default=5,
                      condition='not doContinue',
                      label='Update poses every N epochs')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn abinit_het -h* for help.")

    # --------------------------- STEPS functions -----------------------------
    def runTrainingStep(self):
        run = self._getRun()
        protType = run.protType.get()
        program = "homo" if protType == AB_INITIO_HOMO else "het"
        self._runProgram("abinit_" + program, self._getTrainingArgs(protType))

    def createOutputStep(self):
        """ Creating a set of particles with z_values. """
        run = self._getRun()
        protType = run.protType.get()
        if protType == AB_INITIO_HETERO:
            CryoDrgnProtBase.createOutputStep(self)
        else:
            # Creating output volume
            inputSet = self._getInputParticles()
            vol = Volume()
            vol.setFileName(self.getOutputDir("reconstruct.mrc"))
            vol.setSamplingRate(inputSet.getSamplingRate())

            self._defineOutputs(**{outputs.Volumes.name: vol})
            self._defineSourceRelation(self._getInputParticles(pointer=True), vol)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        run = self._getRun()
        protType = run.getEnumText("protType")
        summary = [f"Training ab initio ({protType}) for {self.numEpochs} epochs."]

        return summary

    def _validate(self):
        errors = super()._validate()

        if not self.doContinue:
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
        run = self._getRun()

        args = [
            self._getFileName('input_parts'),
            f"--ctf {self._getFileName('input_ctfs')}",
            f"-o {self.getOutputDir()}",
            f"-n {self.numEpochs}",
            f"--t-extent {run.searchRange}",
            f"--ps-freq {run.psFreq}",
            "--load latest" if self.doContinue else "",
            f"--datadir {self._getExtraPath('input')}"
        ]

        if protType == AB_INITIO_HETERO:
            args.extend([
                f"--zdim {run.zDim}",
                f"--max-threads {self.numberOfThreads}",
            ])

            if len(self.getGpuList()) > 1:  # only for hetero
                args.append('--multigpu')

        if self.extraParams.hasValue():
            args.append(self.extraParams.get())

        return args

    def _getRun(self):
        return self.continueRun.get() if self.doContinue else self
