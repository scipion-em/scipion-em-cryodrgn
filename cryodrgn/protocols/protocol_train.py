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

import pyworkflow.object as pwobj
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params

from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class CryoDrgnProtTrain(CryoDrgnProtBase):
    """ Protocol to train cryoDRGN neural network. """
    _label = 'training VAE'
    _devStatus = PROD
    _possibleOutputs = CryoDrgnProtBase._possibleOutputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        CryoDrgnProtBase._defineParams(self, form)
        form.getParam('numberOfThreads').default = pwobj.Integer(1)

    def _defineAdvancedParams(self, form):
        form.addSection(label='Advanced')
        group = form.addGroup('Encoder', condition='not doContinue')
        group.addParam('qLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('qDim', params.IntParam, default=1024,
                       label='Number of nodes in hidden layers')

        group = form.addGroup('Decoder', condition='not doContinue')
        group.addParam('pLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('pDim', params.IntParam, default=1024,
                       label='Number of nodes in hidden layers')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn train_vae -h* for help.")

    # --------------------------- STEPS functions -----------------------------
    def runTrainingStep(self):
        self._runProgram('train_vae', self._getTrainingArgs())

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = [f"Training VAE for {self.numEpochs} epochs."]

        return summary

    def _validate(self):
        errors = super()._validate()

        if not self._inputHasAlign():
            errors.append("Input particles have no alignment information!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self):
        run = self.continueRun.get() if self.doContinue else self

        args = [
            self._getFileName('input_parts'),
            f"--poses {self._getFileName('input_poses')}",
            f"--ctf {self._getFileName('input_ctfs')}",
            f"--zdim {run.zDim}",
            f"-o {self.getOutputDir()}",
            f"-n {self.numEpochs}",
            f"--max-threads {self.numberOfThreads}",
            f"--enc-layers {run.qLayers}",
            f"--enc-dim {run.qDim}",
            f"--dec-layers {run.pLayers}",
            f"--dec-dim {run.pDim}",
            "--load latest" if self.doContinue else "",
            f"--datadir {self._getExtraPath('input')}"
        ]

        if run.doWindow:
            args.append(f"--window-r {run.winSize}")

        if not run.doInvert:  # neg. stain only
            args.append('--uninvert-data')

        if run.lazyLoad:
            args.append("--lazy")

        if len(self.getGpuList()) > 1:
            args.append('--multigpu')

        if self._getInputParticles().getXDim() % 8 != 0:
            args.append("--no-amp")

        if self.extraParams.hasValue():
            args.append(self.extraParams.get())

        return args
