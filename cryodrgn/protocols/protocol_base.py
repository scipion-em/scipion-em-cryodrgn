# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es) [2]
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
import re
from glob import glob
from enum import Enum

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pwem.constants import ALIGN_PROJ, ALIGN_NONE
from pwem.protocols import ProtProcessParticles, ProtFlexBase
from pwem.objects import SetOfParticlesFlex, ParticleFlex

from cryodrgn import Plugin
from cryodrgn.constants import WEIGHTS, CONFIG, CRYODRGN


convert = Domain.importFromPlugin('relion.convert', doRaise=True)


class outputs(Enum):
    Particles = SetOfParticlesFlex


class CryoDrgnProtBase(ProtProcessParticles, ProtFlexBase):
    _label = None
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
            'input_parts': self._getExtraPath('input_particles.star'),
            'input_poses': self._getExtraPath('poses.pkl'),
            'input_ctfs': self._getExtraPath('ctf.pkl'),
            'z': self.getOutputDir('z.%(epoch)d.pkl'),
            'z_final': self.getOutputDir('z.pkl'),
            'weights': self.getOutputDir('weights.%(epoch)d.pkl'),
            'weights_final': self.getOutputDir('weights.pkl'),
            'config': self.getOutputDir('config.yaml')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use multiple GPUs - in that case"
                            " set to i.e. *0 1 2*.")

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

        form.addParam('doInvert', params.BooleanParam, default=True,
                      condition='not doContinue',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Are particles white?")

        form.addParam('doWindow', params.BooleanParam, default=True,
                      condition='not doContinue',
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Apply circular mask?")

        form.addParam('winSize', params.FloatParam, default=0.85,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition='doWindow and not doContinue',
                      label="Window size",
                      help="Circular windowing mask inner radius")

        form.addParam('lazyLoad', params.BooleanParam, default=False,
                      condition='not doContinue',
                      label="Use lazy loading?",
                      help="Lazy loading if full dataset is too large to "
                           "fit in memory.")

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

        self._defineAdvancedParams(form)

        form.addParallelSection(threads=16, mpi=0)

    def _defineAdvancedParams(self, form):
        """ Should be defined in subclasses. """
        raise NotImplementedError

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()

        if self.doContinue:
            self._insertFunctionStep(self.continueStep)
        else:
            self._insertFunctionStep(self.convertInputStep)

        self._insertFunctionStep(self.runTrainingStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create the input star, poses and ctf pkl files as expected by cryoDRGN. """
        imgSet = self._getInputParticles()
        alignType = ALIGN_PROJ if self._inputHasAlign() else ALIGN_NONE
        convert.writeSetOfParticles(imgSet,
                                    self._getExtraPath('input_particles.star'),
                                    outputDir=self._getExtraPath(),
                                    alignType=alignType)

        if self._inputHasAlign() and self.getClassName() != "CryoDrgnProtAbinitio":
            self._runProgram('parse_pose_star', self._getParsePosesArgs())
        self._runProgram('parse_ctf_star', self._getParseCtfArgs())

    def continueStep(self):
        """ Copy previous run outputs. """
        prevRun = self.continueRun.get()

        pwutils.cleanPath(self._getExtraPath())
        self.info("Copying previous run results...")
        pwutils.copyTree(prevRun._getExtraPath(), self._getExtraPath())

    def runTrainingStep(self):
        """ Should be implemented in subclasses. """
        raise NotImplementedError

    def createOutputStep(self):
        """ Creating a set of particles with z_values. """
        inputSet = self._getInputParticles()
        zValues = iter(self._getParticlesZvalues())
        outImgSet = self._createSetOfParticlesFlex(progName=CRYODRGN)
        outImgSet.copyInfo(inputSet)
        outImgSet.setHasCTF(inputSet.hasCTF())
        outImgSet.getFlexInfo().setProgName(CRYODRGN)

        for particle, zValue in zip(inputSet, zValues):
            outParticle = ParticleFlex(progName=CRYODRGN)
            outParticle.copyInfo(particle)
            outParticle.getFlexInfo().setProgName(CRYODRGN)

            outParticle.setZFlex(list(zValue))

            outImgSet.append(outParticle)

        outImgSet.getFlexInfo().setAttr(WEIGHTS, self._getFileName('weights_final'))
        outImgSet.getFlexInfo().setAttr(CONFIG, self._getFileName('config'))

        self._defineOutputs(**{outputs.Particles.name: outImgSet})
        self._defineSourceRelation(self._getInputParticles(pointer=True), outImgSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        if self.doContinue:
            if not self.continueRun.hasValue():
                errors.append("Select the input run to continue from!")

            prevEpochs = self.continueRun.get().numEpochs.get()
            if self.numEpochs <= prevEpochs:
                errors.append(f"Number of epochs must be larger than {prevEpochs} "
                              "that are already completed!")

        return errors

    # --------------------------- UTILS functions -----------------------------
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

    def getOutputDir(self, *paths):
        return self._getExtraPath("output", *paths)

    def _getLastEpoch(self):
        """ Return the last iteration number. """
        epoch = None
        epochRegex = re.compile(r'weights.(\d).pkl')
        files = sorted(glob(self._getFileName("weights", epoch=0).replace('0', '*')))
        if files:
            f = files[-1]
            s = epochRegex.search(f)
            if s:
                epoch = int(s.group(1))  # group 1 is a digit iteration number

        return epoch

    def _canContinue(self):
        return self._getLastEpoch() is not None

    def _getInputParticles(self, pointer=False):
        if self.doContinue and self.continueRun.hasValue():
            parts = self.continueRun.get().inputParticles
        else:
            parts = self.inputParticles

        return parts if pointer else parts.get()

    def _inputHasAlign(self):
        return self._getInputParticles().hasAlignmentProj()
