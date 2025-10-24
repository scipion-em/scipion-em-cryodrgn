# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Eduardo García Delgado (eduardo.garcia@cnb.csic.es) [2]
# *
# * [1] MRC Laboratory of Molecular Biology (MRC-LMB)
# * [2] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import pickle
from pwem.constants import ALIGN_PROJ, ALIGN_NONE
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.plugin import Domain
from pyworkflow.constants import PROD
from pwem.protocols import ProtProcessParticles, ProtFlexBase
from .. import Plugin
from ..constants import *

convert = Domain.importFromPlugin('relion.convert', doRaise=True)

class CryoDrgnProtTrain(ProtProcessParticles, ProtFlexBase):
    """ Protocol to train cryoDRGN neural network. """

    _label = 'training VAE'
    _devStatus = PROD

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
            'input_parts': self._getExtraPath('input_particles.star'),
            'input_poses': self._getExtraPath('poses.pkl'),
            'input_ctfs': self._getExtraPath('ctf.pkl'),
            'z': self._getOutputDir('z.%(epoch)d.pkl'),
            'z_final': self._getOutputDir('z.pkl'),
            'weights': self._getOutputDir('weights.%(epoch)d.pkl'),
            'weights_final': self._getOutputDir('weights.pkl'),
            'config': self._getOutputDir('config.yaml')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('doContinue', params.BooleanParam, default=False,
                      label="Continue previous run?",
                      help="Training will resume from the latest epoch.")

        form.addParam('continueRun', params.PointerParam,
                      condition='doContinue', important=True,
                      pointerClass='CryoDrgnProtTrain, CryoDrgnProtAbinitio',
                      label="Previous run", allowsNull=True)

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      label="Input particles", important=True,
                      help='Select a set of particles from a consensus C1 '
                           '3D refinement.')

        form.addSection(label='Training')
        form.addParam('zDim', params.IntParam, default=8,
                      condition='not doContinue',
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

        group = form.addGroup('Encoder', condition='not doContinue', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('qLayers', params.IntParam, default=3,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Number of hidden layers')
        group.addParam('qDim', params.IntParam, default=1024,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Number of nodes in hidden layers')

        group = form.addGroup('Decoder', condition='not doContinue', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('pLayers', params.IntParam, default=3,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Number of hidden layers')
        group.addParam('pDim', params.IntParam, default=1024,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Number of nodes in hidden layers')

        form.addSection('Network parameters')
        form.addParam('batchSize', params.IntParam, default=8,
                      condition='not doContinue',
                      label="Batch size",
                      help="Minibatch size for training")

        form.addParam('learningRate', params.FloatParam, default=0.0001,
                      condition='not doContinue',
                      label="Learning rate",
                      help="Learning rate in Adam optimizer")

        form.addParam('weightDecay', params.FloatParam, default=0.0,
                      condition='not doContinue',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Weight decay",
                      help="Weight decay for Adam optimizer")

        form.addParam('doInvert', params.BooleanParam, default=True,
                      condition='not doContinue',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Invert",
                      help="Invert input particles?")

        form.addParam('doWindow', params.BooleanParam, default=False,
                      condition='not doContinue',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Circular mask",
                      help="Apply for circular mask?")

        form.addParam('winSize', params.FloatParam, default=0.85,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='doWindow and not doContinue',
                      label="Window size",
                      help="Circular windowing mask inner radius")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()

        if self.doContinue:
            self._insertFunctionStep(self.continueStep, needsGPU=False)
        else:
            self._insertFunctionStep(self.convertInputStep, needsGPU=False)

        self._insertFunctionStep(self.runTrainingStep, needsGPU=True)

    # --------------------------- STEPS functions -----------------------------

    def continueStep(self):
        """ Copy previous run outputs. """
        prevRun = self.continueRun.get()
        pwutils.cleanPath(self._getExtraPath())
        pwutils.copyTree(prevRun._getExtraPath(), self._getExtraPath())

    def convertInputStep(self):
        """ Create the input star, poses and ctf pkl files as expected by cryoDRGN. """
        imgSet = self._getInputParticles()
        alignType = ALIGN_PROJ if self._inputHasAlign() else ALIGN_NONE
        convert.writeSetOfParticles(imgSet,
                                    self._getExtraPath('input_particles.star'),
                                    outputDir=self._getExtraPath(),
                                    alignType=alignType)

        self._runProgram('parse_pose_star', self._getParsePosesArgs())
        self._runProgram('parse_ctf_star', self._getParseCtfArgs())

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

        if not self.doContinue:
            if self.zDim == 1:
                errors.append("Latent variable must be >1 for "
                              "heterogeneous reconstruction")
        else:
            if not self.continueRun.hasValue():
                errors.append("Select the input run to continue from!")

            prevEpochs = self.continueRun.get().numEpochs.get()
            if self.numEpochs <= prevEpochs:
                errors.append(f"Number of epochs must be larger than {prevEpochs} "
                              "that are already completed!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self):
        run = self.continueRun.get() if self.doContinue else self

        args = [
            self._getFileName('input_parts'),
            f"--poses {self._getFileName('input_poses')}",
            f"--ctf {self._getFileName('input_ctfs')}",
            f"--zdim {run.zDim}",
            f"-o {self._getOutputDir()}",
            f"-n {self.numEpochs}",
            f"--lr {self.learningRate}",
            f"--wd {self.weightDecay}",
            f"--batch-size {self.batchSize}",
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

        if len(self.getGpuList()) > 1:
            args.append('--multigpu')

        if self._getInputParticles().getXDim() % 8 != 0:
            args.append("--no-amp")

        return args

    def _getParsePosesArgs(self):
        args = [
            self._getFileName('input_parts'),
            f"-o {self._getFileName('input_poses')}"
        ]

        return args

    def _getParseCtfArgs(self):
        args = [
            self._getFileName('input_parts'),
            f"-o {self._getFileName('input_ctfs')}",
            "--ps 0"  # required due to cryodrgn parsing bug
        ]

        return args

    def _runProgram(self, program, args):
        gpus = ','.join(str(i) for i in self.getGpuList())
        self.runJob(Plugin.getProgram(program, gpus), ' '.join(args))

    def _getParticlesZvalues(self):
        """
        Read from z.pkl file the particles z_values
        :return: a numpy array with the particles z_values
        """
        zEpochFile = self._getFileName("z_final")
        with open(zEpochFile, 'rb') as f:
            zValues = pickle.load(f)

        return zValues

    def _setZValues(self, item, row=None):
        item.getFlexInfo().setProgName(CRYODRGN)
        item.setZFlex(list(row))

    def _getOutputDir(self, *paths):
        return self._getExtraPath("output", *paths)

    def _getInputParticles(self, pointer=False):
        if self.doContinue and self.continueRun.hasValue():
            parts = self.continueRun.get().inputParticles
        else:
            parts = self.inputParticles

        return parts if pointer else parts.get()

    def _inputHasAlign(self):
        return self._getInputParticles().hasAlignmentProj()
