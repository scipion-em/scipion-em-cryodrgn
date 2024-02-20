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
import numpy as np

from pyworkflow.plugin import Domain
from pyworkflow.constants import PROD
import pyworkflow.protocol.params as params
from pwem.constants import ALIGN_PROJ, ALIGN_NONE
from pwem.protocols import ProtProcessParticles
from pwem.objects import SetOfParticles

from cryodrgn import Plugin

convert = Domain.importFromPlugin('relion.convert', doRaise=True)


class outputs(Enum):
    Particles = SetOfParticles


class CryoDrgnProtPreprocess(ProtProcessParticles):
    """ Protocol to downsample a particle stack. """
    _label = 'preprocess particles'
    _devStatus = PROD
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden('usePreprocess', params.BooleanParam, default=True)
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select a set of particles from a consensus C1 '
                           '3D refinement.')

        form.addParam('doScale', params.BooleanParam, default=True,
                      label='Downsample particles?')

        form.addParam('scaleSize', params.IntParam, default=128,
                      condition='doScale',
                      validators=[params.Positive],
                      label='New box size (px)',
                      help='New box size in pixels, must be even.')

        form.addParam('chunk', params.IntParam, default=0,
                      label='Split in chunks',
                      help='Chunk size (in # of images) to split '
                           'particle stack when saving.')

        form.addParallelSection(threads=16, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.runDownSampleStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a star file as expected by cryoDRGN."""
        imgSet = self._getInputParticles()
        # Create links to binary files and write the relion .star file
        alignType = ALIGN_PROJ if self._inputHasAlign() else ALIGN_NONE
        convert.writeSetOfParticles(imgSet,
                                    self._getTmpPath('input_particles.star'),
                                    outputDir=self._getTmpPath(),
                                    alignType=alignType)

    def runDownSampleStep(self):
        self._runProgram('downsample', self._getArgs())

    def createOutputStep(self):
        inputSet = self._getInputParticles()
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(inputSet)

        newSampling = self._getSamplingRate()
        imgSet.setSamplingRate(newSampling)

        itemIter = self._getOutputFn(inputSet.getSize(), self.chunk.get())
        imgSet.copyItems(inputSet,
                         itemDataIterator=itemIter,
                         updateItemCallback=self._updateLocation)
        self._defineOutputs(**{outputs.Particles.name: imgSet})
        self._defineTransformRelation(self.inputParticles, imgSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []

        particles = self._getInputParticles()

        if self.doScale and self.scaleSize > particles.getXDim():
            errors.append("You cannot upscale particles!")

        if self._getBoxSize() % 2 != 0:
            errors.append("Box size must be even!")

        return errors

    def _warnings(self):
        warnings = []

        if not self._inputHasAlign():
            warnings.append("Input particles have no alignment, you will only "
                            "be able to use the output for *ab initio* training!")

        if self._getBoxSize() % 8 != 0:
            warnings.append("CryoDRGN mixed-precision (AMP) training will "
                            "require box size divisible by 8. Alternatively, "
                            "you will have to provide --no-amp option.")

        return warnings

    # --------------------------- UTILS functions -----------------------------
    def _getArgs(self):
        newBox = self._getBoxSize()
        args = [
            self._getTmpPath('input_particles.star'),
            f"-o {self._getExtraPath('particles.%d.mrcs' % newBox)}",
            f"--datadir {self._getTmpPath('input')}",
            f"-D {newBox}",
            f"--max-threads {self.numberOfThreads}"
        ]

        if self.chunk > 0:
            args.append(f"--chunk {self.chunk}")

        return args

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _getBoxSize(self):
        if self.doScale:
            return self.scaleSize.get()
        else:
            return self._getInputParticles().getXDim()

    def _getSamplingRate(self):
        inputSet = self._getInputParticles()
        oldSampling = inputSet.getSamplingRate()
        scaleFactor = self._getScaleFactor()

        return oldSampling * scaleFactor

    def _getScaleFactor(self):
        return self._getInputParticles().getXDim() / self._getBoxSize()

    def _inputHasAlign(self):
        return self._getInputParticles().hasAlignmentProj()

    def _runProgram(self, program, args):
        self.runJob(Plugin.getProgram(program), ' '.join(args))

    def _getOutputFn(self, totalSize, chunkSize):
        newBox = self._getBoxSize()
        if chunkSize == 0:
            indexes = np.arange(totalSize)
            fnames = np.full(totalSize, f"particles.{newBox}.mrcs")
        else:
            q, mod = divmod(totalSize, chunkSize)
            chunks = q * [chunkSize] + [mod]
            indexes = np.concatenate([np.arange(i) for i in chunks])
            fnames = np.concatenate([np.full(i, f"particles.{newBox}.{n}.mrcs") for n, i in enumerate(chunks)])

        for index, fn in zip(indexes, fnames):
            yield index+1, self._getExtraPath(fn)

    def _updateLocation(self, item, row):
        """ Update the output item location.
        :item: output item
        :row: new (index, fn) output location
        """
        item.setLocation(row)
        invFactor = 1 / self._getScaleFactor()

        if invFactor != 1.0:
            if item.hasCoordinate():
                item.scaleCoordinate(invFactor)
            if item.hasTransform():
                item.getTransform().scaleShifts(invFactor)
