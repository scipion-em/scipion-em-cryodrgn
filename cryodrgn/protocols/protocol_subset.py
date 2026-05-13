# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk) [1]
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es) [2]
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
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles, ProtFlexBase
from pyworkflow.constants import PROD
from cryodrgn.constants import CRYODRGN

class CryoDrgnProtSubset(ProtProcessParticles, ProtFlexBase):
    """
    Creates a subset of cryoDRGN particles based on externally selected
    particle indices. The protocol is intended for extracting biologically
    meaningful conformational populations identified during latent space
    analysis and interactive exploration of cryoDRGN heterogeneity results.

    AI Generated:

    CryoDRGN Particle Subset (CryoDrgnProtSubset) — User Manual
        Overview

        The CryoDRGN Particle Subset protocol is designed to generate a new
        particle set from a previously analyzed cryoDRGN dataset using a
        predefined selection of particle indices. Its main purpose is to
        isolate specific conformational populations, structural states, or
        regions of interest discovered during latent space exploration.

        In practical cryo-EM workflows, this protocol becomes particularly
        useful after training a cryoDRGN model and inspecting the resulting
        latent landscape. Researchers often identify clusters, trajectories,
        or continuous conformational regions that correspond to biologically
        relevant molecular states. This protocol allows those selected
        particles to be separated into a dedicated subset for further
        refinement, reconstruction, classification, or interpretation.

        Biological Motivation and Typical Applications

        Continuous heterogeneity analysis frequently reveals that particle
        datasets contain multiple conformational states connected through
        smooth structural transitions. Instead of treating the dataset as a
        single homogeneous population, users may wish to isolate particles
        corresponding to specific functional states, ligand-binding
        conformations, domain motions, or assembly intermediates.

        This protocol provides a practical bridge between exploratory latent
        space analysis and downstream structural biology workflows. Selected
        subsets can subsequently be refined independently, reconstructed at
        higher resolution, or compared against biochemical hypotheses.

        Typical applications include isolating open and closed states of
        molecular machines, separating flexible domain arrangements,
        identifying rare conformations, or extracting particles along a
        continuous reaction trajectory inferred from cryoDRGN analysis.

        Input Particle Requirements

        The protocol requires particles containing cryoDRGN flexibility
        information. These particles usually originate from cryoDRGN training
        or ab initio heterogeneity analysis workflows and already contain
        latent space annotations describing their conformational placement.

        It is biologically important that the subset selection corresponds to
        the same particle set used during cryoDRGN analysis. Mismatched
        particle indexing between datasets can lead to invalid selections and
        biologically meaningless subsets.

        The particle selection itself is provided through an external file
        containing particle indices. In most workflows, these indices are
        generated interactively during latent space exploration using
        visualization notebooks or custom analysis tools.

        Latent Space Selection Strategies

        The biological meaning of the resulting subset depends entirely on
        how particles are selected from the latent landscape. Different
        selection strategies can emphasize distinct aspects of molecular
        variability.

        Cluster-based selection is commonly used when the latent space
        contains well-separated conformational populations. In these cases,
        the resulting subsets often correspond to discrete structural states
        that can be independently reconstructed and interpreted.

        Trajectory-based selection is useful for studying gradual conformational
        changes. Selecting particles along a continuous path through latent
        space may reveal intermediate states involved in molecular motions or
        functional transitions.

        Density-based selection may also help isolate rare or transient
        conformations that are underrepresented in the original dataset but
        biologically important.

        Interpretation of the Output

        The resulting output is a new particle set containing only the
        selected particles while preserving the associated flexibility
        information and metadata. This allows the subset to remain compatible
        with downstream cryoDRGN analyses as well as conventional cryo-EM
        refinement workflows.

        Biologically, the output subset should be interpreted as a focused
        representation of a particular region of conformational space rather
        than a completely homogeneous population. Depending on the selection
        criteria, residual variability may still remain within the subset.

        The protocol preserves important experimental metadata such as CTF
        information and particle relationships, enabling further refinement
        and reconstruction without loss of contextual information.

        Practical Recommendations

        In routine biological analysis, it is often beneficial to begin with
        broad exploratory selections before progressively refining the subset
        boundaries. Visual inspection of latent distributions and reconstructed
        volumes is strongly recommended to confirm that selected particles
        correspond to meaningful structural variability.

        Overly narrow selections may produce subsets with insufficient
        particle counts for high-resolution refinement, whereas excessively
        broad selections may reintroduce heterogeneity and blur structural
        features. Balancing structural purity and particle number is therefore
        an important practical consideration.

        Users should also verify that the selected indices correspond exactly
        to the intended dataset. Incorrect indexing is one of the most common
        causes of invalid subsets and misleading biological interpretation.

        Final Perspective

        Particle subsetting is a critical step in transforming continuous
        heterogeneity analysis into biologically interpretable structural
        models. By isolating regions of latent space associated with specific
        conformational behaviors, the CryoDRGN Particle Subset protocol helps
        researchers move from abstract latent representations toward concrete
        structural and functional interpretation.

        Careful selection strategy, validation of reconstructed subsets, and
        thoughtful biological interpretation remain essential for extracting
        reliable insights from flexible cryo-EM datasets.
    """

    _label = "particles subset"
    _devStatus = PROD
    doContinue = False

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticlesFlex',
                      label="Input particles with Flex info", important=True,
                      help="Select a set of output particles from CryoDrgn "
                           "training or ab-initio protocol.")

        form.addParam('pklFile', params.FileParam, important=True,
                      filter="*.pkl", default='',
                      label='Choose *.pkl file with particle selection',
                      help="This usually comes from filtering particles using "
                           "the Jupyter notebook.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def createOutputStep(self):
        inputSet = self._getInputParticles()
        outImgSet = self._createSetOfParticlesFlex(progName=CRYODRGN)
        outImgSet.copyInfo(inputSet)
        outImgSet.setHasCTF(inputSet.hasCTF())
        outImgSet.copyItems(inputSet, self._updateItem,
                            itemDataIterator=iter(range(inputSet.getSize())))

        self._defineOutputs(Particles=outImgSet)
        self._defineSourceRelation(self._getInputParticles(pointer=True), outImgSet)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        if self.isFinished():
            summary.append(
                f"Input particles: {self._getInputParticles().getSize()}\n"
                f"Selected particles: {self.Particles.getSize()}")

        return summary

    def _validate(self):
        errors = []

        inputSize = self._getInputParticles().getSize()
        subsetSize = len(self._getParticlesIndices())

        if subsetSize > inputSize:
            errors.append("Subset size is larger than the input set size.")
        if max(self._getParticlesIndices()) > inputSize-1:
            errors.append("Subset has particle indices bigger "
                          "than the input set size. Make sure you are "
                          "selecting matching sets!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getParticlesIndices(self):
        """ Get zero-based indices of particles. """
        with open(self.pklFile.get(), "rb") as f:
            x = pickle.load(f)
        return x

    def _updateItem(self, item, index):
        if index not in self._getParticlesIndices():
            item._appendItem = False
