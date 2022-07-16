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

import pyworkflow.utils as pwutils
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles

from cryodrgn import Plugin
from cryodrgn.constants import *


class CryoDrgnProtTrain(ProtProcessParticles):
    """
    Protocol to train cryoDRGN neural network.

    Find more information at https://github.com/zhonge/cryodrgn
    """
    _label = 'training'
    _devStatus = PROD

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden('doInvert', params.BooleanParam, default=True)
        form.addSection(label='Input')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass="CryoDrgnParticles",
                      label='CryoDrgn particles')

        form.addParam('zDim', params.IntParam, default=8,
                      validators=[params.Positive],
                      label='Dimension of latent variable',
                      help='It is recommended to first train on lower '
                           'resolution images (e.g. D=128) with '
                           '--zdim 1 and with --zdim 10 using the '
                           'default architecture (fast).')
        form.addParam('numEpochs', params.IntParam, default=20,
                      label='Number of epochs',
                      help='The number of epochs refers to the number '
                           'of full passes through the dataset for '
                           'training, and should be modified depending '
                           'on the number of particles in the dataset. '
                           'For a 100k particle dataset, the above '
                           'settings required ~6 min per epoch for D=128 '
                           'images + default architecture, ~12 min/epoch '
                           'for D=128 images + large architecture, and ~47 '
                           'min per epoch for D=256 images + large architecture.')

        form.addSection(label='Advanced')
        group = form.addGroup('Encoder')
        group.addParam('qLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('qDim', params.IntParam, default=1024,
                       label='Number of nodes in hidden layers')

        group = form.addGroup('Decoder')
        group.addParam('pLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('pDim', params.IntParam, default=1024,
                       label='Number of nodes in hidden layers')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn train_vae -h* for help.")

        form.addSection(label='Analysis')

        form.addParam('viewEpoch', params.EnumParam,
                      choices=['last', 'selection'], default=EPOCH_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Epoch to analyze")

        form.addParam('epochNum', params.IntParam,
                      condition='viewEpoch==%d' % EPOCH_SELECTION,
                      label="Epoch number")

        form.addParam('ksamples', params.IntParam, default=20,
                      label='Number of K-means samples to generate',
                      help="*cryodrgn analyze* uses the k-means clustering "
                           "algorithm to partition the latent space into "
                           "regions (by default k=20 regions), and generate a "
                           "density map from the center of each of these "
                           "regions. The goal is to provide a tractable number "
                           "of representative density maps to visually inspect. ")

        form.addParallelSection(threads=16, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runTrainingStep')
        if self.viewEpoch == EPOCH_LAST:
            epoch = self.numEpochs.get() - 1
        else:
            epoch = self.epochNum.get()
        self._insertFunctionStep('runAnalysisStep', epoch)

    # --------------------------- STEPS functions -----------------------------
    def runTrainingStep(self):
        # Create output folder
        pwutils.cleanPath(self.getOutputDir())
        pwutils.makePath(self.getOutputDir())

        # Call cryoDRGN with the appropriate parameters.
        self._runProgram('train_vae', self._getTrainingArgs())

    def runAnalysisStep(self, epoch):
        """ Run analysis step.
        Args:
            epoch: epoch number to be analyzed.
        """
        self._runProgram('analyze', self._getAnalyzeArgs(epoch))

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["Training VAE for %d epochs." % self.numEpochs]

        return summary

    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self):
        parts = self.inputParticles.get()

        args = [
            parts.filename.get(),
            '--poses %s' % parts.poses,
            '--ctf %s' % parts.ctfs,
            '--zdim %d' % self.zDim,
            '-o %s ' % self.getOutputDir(),
            '-n %d' % self.numEpochs,
            '--preprocessed',
            '--max-threads %d ' % self.numberOfThreads,
            '--enc-layers %d' % self.qLayers,
            '--enc-dim %d' % self.qDim,
            '--dec-layers %d' % self.pLayers,
            '--dec-dim %d' % self.pDim
        ]

        if not Plugin.versionGE(V1_0_0):
            args.append('--relion31')

        if len(self.getGpuList()) > 1:
            args.append('--multigpu')

        if self.extraParams.hasValue():
            args.append(self.extraParams.get())

        return args

    def _getAnalyzeArgs(self, epoch):
        return [
            self.getOutputDir(),
            str(epoch),
            '--Apix %0.3f' % self.inputParticles.get().getSamplingRate(),
            '--ksample %d' % self.ksamples,
        ]

    def _runProgram(self, program, args):
        gpus = ','.join(str(i) for i in self.getGpuList())
        self.runJob(Plugin.getProgram(program, gpus), ' '.join(args))

    def getOutputDir(self):
        return self._getPath('output')

    def getEpochZFile(self, epoch):
        return os.path.join(self.getOutputDir(), 'z.%d.pkl' % epoch)

    def getLastEpoch(self):
        epoch = -1

        while os.path.exists(self.getEpochZFile(epoch + 1)):
            epoch += 1

        return epoch
