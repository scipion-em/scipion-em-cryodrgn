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
    """
    Preprocesses cryo-EM particle stacks for cryoDRGN analysis by
    optionally downsampling particles and preparing them for neural
    network training workflows. The protocol standardizes particle
    dimensions and sampling properties while preserving the metadata
    required for downstream heterogeneous reconstruction and latent
    space analysis.

    AI Generated:

    CryoDRGN Particle Preprocessing (CryoDrgnProtPreprocess) — User Manual
        Overview

        The CryoDRGN Particle Preprocessing protocol prepares cryo-EM
        particle stacks for efficient use in cryoDRGN neural network
        workflows. Its primary purpose is to generate particle datasets
        with dimensions and sampling properties suitable for deep learning
        training while maintaining consistency with the original
        experimental information.

        In practical cryo-EM workflows, preprocessing is often one of the
        first steps before heterogeneous reconstruction or latent space
        analysis. Large particle box sizes can significantly increase GPU
        memory usage and computational cost during neural network training.
        Downsampling particles allows users to reduce these requirements
        while preserving the structural features necessary for studying
        conformational variability.

        Biological Motivation and Practical Context

        CryoDRGN workflows are commonly used to analyze continuous
        heterogeneity, flexible molecular motions, and multiple structural
        states within cryo-EM datasets. In these analyses, computational
        efficiency becomes particularly important because neural network
        training may require many epochs and large particle populations.

        Downsampling provides a balance between structural detail and
        computational feasibility. Lower-resolution representations are
        often sufficient during exploratory analysis or initial latent
        space training, especially when the goal is to identify major
        conformational trends rather than high-resolution features.

        Biological users frequently begin with reduced box sizes to
        accelerate experimentation and parameter optimization before
        training more computationally demanding models on higher-resolution
        particles.

        Input Requirements and Data Consistency

        The protocol accepts a set of cryo-EM particles that typically
        originate from consensus refinements or other reconstruction
        workflows. Ideally, particles should already be reasonably aligned
        and centered prior to preprocessing.

        Alignment information is particularly important for standard
        cryoDRGN heterogeneous reconstruction workflows because particle
        orientations guide the relationship between experimental images and
        reconstructed conformational states. When alignment information is
        absent, the resulting dataset may still be suitable for ab initio
        training, although downstream interpretation becomes more limited.

        The preprocessing stage also preserves important metadata such as
        coordinates and geometric relationships. This ensures that the
        processed particles remain compatible with later refinement and
        reconstruction steps.

        Downsampling Strategy

        The central operation of this protocol is particle downsampling.
        Reducing the box size decreases memory usage, accelerates disk
        access, and substantially shortens neural network training time.

        From a biological perspective, the chosen box size should remain
        large enough to preserve the structural features relevant to the
        intended analysis. Excessive downsampling may remove subtle
        conformational signals or blur flexible domains that are important
        for interpreting molecular dynamics.

        A moderate reduction in particle dimensions is often a good
        starting point for exploratory latent space analysis. Once
        meaningful conformational organization has been identified, users
        may decide to retrain using larger particle dimensions for improved
        structural detail.

        Computational Considerations

        Neural network training efficiency strongly depends on particle
        dimensions. Smaller box sizes allow larger batch sizes and more
        stable GPU memory usage, which can significantly accelerate
        experimentation.

        The protocol also supports splitting particle stacks into smaller
        chunks. This is useful when working with extremely large datasets
        or limited storage environments because it improves manageability
        and reduces the size of individual output files.

        Certain neural network optimizations benefit from box sizes that
        are divisible by specific values. Choosing dimensions compatible
        with mixed-precision training can substantially improve training
        performance on modern GPU hardware.

        Interpretation of the Output

        The output consists of a processed particle stack with updated
        sampling information and geometry adjusted to match the new image
        dimensions. The resulting particles remain associated with their
        original metadata and can be directly used in downstream cryoDRGN
        workflows.

        Biologically, the processed dataset should be interpreted as a
        computationally optimized representation of the original particle
        population rather than a fundamentally altered dataset. Structural
        relationships and conformational variability are preserved within
        the limits imposed by the selected resolution and box size.

        Practical Recommendations

        In routine practice, many users begin with moderate downsampling
        during exploratory heterogeneity analysis. This allows rapid
        testing of training parameters and latent space organization before
        investing computational resources into larger-scale models.

        Care should be taken to avoid reducing particle dimensions below
        the scale required to visualize the biological variability of
        interest. Flexible domains, ligand-binding regions, or subtle
        conformational transitions may become difficult to resolve if the
        box size is excessively reduced.

        Users should also verify that the selected box size remains even
        and computationally compatible with the intended cryoDRGN training
        configuration.

        Final Perspective

        Particle preprocessing is a foundational preparation step for
        neural network-based cryo-EM heterogeneity analysis. By balancing
        computational efficiency with preservation of biologically relevant
        structural information, the CryoDRGN Particle Preprocessing
        protocol enables practical and scalable exploration of molecular
        flexibility within large cryo-EM datasets.

        Thoughtful selection of particle dimensions and preprocessing
        strategy can strongly influence both the efficiency of training and
        the interpretability of downstream conformational analysis.
    """
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
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.runDownSampleStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

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
