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
    """
    Trains a cryoDRGN variational autoencoder model to learn continuous
    structural heterogeneity directly from cryo-EM particle images. The
    protocol enables reconstruction of flexible molecular landscapes by
    embedding particles into a latent space that captures conformational
    variability across the dataset.

    AI Generated:

    CryoDRGN Training (CryoDrgnProtTrain) — User Manual
        Overview

        The CryoDRGN Training protocol is designed to analyze continuous
        conformational heterogeneity in cryo-electron microscopy datasets
        using deep learning approaches based on variational autoencoders.
        Instead of describing molecular variability through a small number
        of discrete classes, the protocol models structural differences as
        a continuous latent landscape. This allows researchers to study
        gradual motions, domain rearrangements, and flexible transitions
        that are often difficult to capture using conventional classification
        strategies.

        In practical cryo-EM workflows, this protocol is especially useful
        for studying dynamic macromolecular assemblies, flexible enzymes,
        membrane proteins, ribosomes, and complexes undergoing functional
        transitions. By learning a compact latent representation of particle
        variability, the protocol enables downstream visualization and
        interpretation of conformational continua rather than isolated
        structural snapshots.

        Inputs and Biological Context

        The protocol requires a set of particles associated with a consensus
        reconstruction and containing valid projection alignment information.
        The quality of these alignments is biologically important because
        accurate orientation estimates strongly influence the ability of the
        neural network to separate structural variability from alignment
        uncertainty.

        In most biological applications, the input dataset should correspond
        to a single biochemical composition or assembly state. Mixing
        particles from unrelated complexes, severe contaminants, or strongly
        heterogeneous stoichiometries may produce latent representations
        that are difficult to interpret biologically.

        The protocol also supports continuation from a previous training run.
        This capability is particularly valuable for large cryo-EM datasets
        where extended optimization may be required to adequately capture
        complex conformational landscapes. Continuing a previous run allows
        refinement of the latent representation without restarting the entire
        learning process.

        Latent Space Representation

        One of the most important parameters is the dimensionality of the
        latent space. This latent representation defines how structural
        variability is encoded internally. Smaller dimensions typically
        capture only the largest conformational motions, while larger latent
        spaces can describe more subtle and complex variability.

        For exploratory biological analyses, lower-dimensional latent spaces
        are often easier to interpret visually because they tend to reveal
        dominant motions such as hinge movements or domain opening events.
        Higher-dimensional representations may capture richer variability,
        but interpretation becomes progressively more difficult and may
        require additional downstream analysis.

        The latent representation should not be interpreted as a direct
        physical coordinate system. Instead, it reflects a mathematical
        embedding of structural variability inferred from the particle data.
        Nearby points in the latent space generally correspond to similar
        conformations, while distant points represent more distinct states.

        Neural Network Architecture

        The protocol allows advanced control over both encoder and decoder
        architectures. These neural network components determine how particle
        images are compressed into latent coordinates and subsequently used
        to reconstruct structural information.

        Larger and deeper architectures may improve the representation of
        complex structural variability, particularly for high-resolution
        datasets or highly flexible systems. However, increasing network
        complexity also raises computational cost and memory requirements.
        For many biological applications, default architectures provide a
        good balance between performance and stability.

        Users working with extremely heterogeneous systems or very large
        particle datasets may benefit from experimenting with deeper models,
        while smaller datasets may train more reliably with simpler
        architectures.

        Training Strategy and Optimization

        The training process iteratively refines the neural network over
        multiple epochs. Each epoch corresponds to a complete pass through
        the dataset. The appropriate number of epochs depends on dataset
        size, structural complexity, particle quality, and desired level
        of convergence.

        Batch size influences both computational efficiency and optimization
        stability. Larger batches may accelerate training on powerful GPUs,
        while smaller batches are often more memory efficient and stable on
        limited hardware resources.

        The learning rate controls how aggressively the optimization updates
        the neural network parameters. Excessively high learning rates may
        destabilize training, whereas very small values may lead to slow
        convergence. In most biological workflows, moderate default values
        provide reliable behavior.

        Weight decay can optionally regularize the optimization process and
        reduce overfitting. This may become useful when training on smaller
        datasets or when the learned latent space appears excessively noisy.

        Particle Preprocessing and Image Conditioning

        The protocol provides optional preprocessing operations that influence
        how particle images are interpreted during training. One option
        applies image inversion, which is particularly relevant depending on
        whether particles appear with dark or bright contrast relative to
        the background.

        Another important option is the application of a circular mask.
        Masking helps focus the neural network on biologically meaningful
        regions while suppressing background solvent noise. In many cryo-EM
        datasets, especially those containing flexible peripheral regions,
        masking improves training stability and enhances interpretability of
        the latent landscape.

        From a biological perspective, the mask radius should include the
        relevant molecular density while avoiding excessive solvent area.
        Overly tight masking may suppress meaningful flexible motions, while
        excessively loose masking may allow noise to dominate the learning
        process.

        GPU Acceleration and Computational Considerations

        CryoDRGN training is computationally intensive and benefits strongly
        from GPU acceleration. The protocol supports execution on one or
        multiple GPUs, allowing efficient handling of large cryo-EM datasets
        and high-dimensional latent spaces.

        Runtime depends heavily on particle count, image size, network
        architecture, and latent dimensionality. Large datasets with high
        resolution images may require substantial GPU memory and extended
        training times. Users should therefore balance biological ambition
        with available computational resources.

        Outputs and Biological Interpretation

        The primary result of the protocol is a trained latent representation
        describing the conformational variability present in the particle
        dataset. Each particle becomes associated with a coordinate in the
        latent space, enabling downstream visualization, clustering, and
        trajectory analysis.

        Biologically, these latent coordinates can reveal continuous motions,
        conformational pathways, and relationships between structural states.
        Researchers often use these embeddings to identify transition
        pathways, reconstruct representative conformations, or explore energy
        landscapes associated with molecular function.

        The protocol also generates trained model parameters and associated
        metadata required for later reconstruction and visualization steps.
        These outputs serve as the foundation for subsequent exploration of
        heterogeneous structural ensembles.

        Practical Recommendations

        In routine biological practice, it is often advisable to begin with
        moderate image sizes, conservative latent dimensionality, and default
        network architectures. This approach provides rapid initial insight
        into the heterogeneity landscape before investing computational
        resources into larger and more detailed models.

        When the latent space appears noisy or poorly organized, improving
        particle quality, refining consensus alignments, or applying
        appropriate masking often provides larger benefits than increasing
        model complexity. Similarly, biologically meaningful variability is
        usually easier to interpret when the dataset is relatively clean and
        compositionally homogeneous.

        For highly dynamic systems, exploring multiple latent dimensions and
        comparing the resulting landscapes can help distinguish robust
        conformational signals from optimization artifacts.

        Final Perspective

        Continuous heterogeneity analysis represents a major conceptual shift
        in cryo-EM structural biology because it allows molecular flexibility
        to be modeled as a continuum rather than a collection of isolated
        classes. The CryoDRGN Training protocol provides a framework for
        uncovering these complex conformational landscapes directly from
        experimental particle images.

        Successful biological interpretation depends not only on neural
        network optimization, but also on careful dataset preparation,
        reliable particle alignments, appropriate preprocessing, and critical
        interpretation of the resulting latent space. When used thoughtfully,
        the protocol can reveal biologically meaningful motions that are
        inaccessible through traditional discrete reconstruction approaches.
    """

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
