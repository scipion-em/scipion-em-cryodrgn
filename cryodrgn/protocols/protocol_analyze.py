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

import os
import numpy as np
import pickle
import pyworkflow.utils as pwutils
from pyworkflow.object import *
import pyworkflow.protocol.params as params
from pyworkflow.constants import PROD
import pwem.objects as emobj
from pwem.emlib import *
from pwem.emlib.image import ImageHandler
from pwem.protocols import ProtProcessParticles, ProtFlexBase
from .. import Plugin
from ..constants import *


class CryoDrgnProtAnalyze(ProtProcessParticles, ProtFlexBase):
    """
    Visualizes latent conformational space and generates representative
    3D volumes from CryoDRGN training results. The protocol enables the
    exploration of structural heterogeneity by analyzing latent embeddings,
    reconstructing density maps, and organizing conformational variability
    into interpretable trajectories and clusters.

    AI Generated:

    Analyze Results (CryoDrgnProtAnalyze) - User Manual
        Overview

        The Analyze Results protocol is designed to interpret and explore
        conformational variability learned by CryoDRGN models. After a neural
        network has been trained on cryo-EM particles, the resulting latent
        space contains information describing continuous or discrete structural
        changes within the dataset. This protocol transforms those latent
        representations into biologically meaningful visualizations and
        representative density maps.

        In practical cryo-EM studies, structural heterogeneity is often one of
        the central scientific questions. Flexible assemblies, dynamic enzymes,
        membrane proteins, and large molecular machines frequently populate a
        continuum of conformational states rather than a single rigid structure.
        This protocol allows users to inspect that variability, identify
        dominant conformations, and reconstruct representative maps suitable
        for downstream interpretation or refinement.

        Inputs and General Workflow

        The protocol requires a previously completed CryoDRGN training or
        ab initio reconstruction together with the associated particle set.
        Users may choose to analyze either the final training epoch or a
        specific epoch from the training trajectory. This flexibility is useful
        because latent spaces can evolve significantly during optimization, and
        intermediate epochs may occasionally provide cleaner or more stable
        conformational organization.

        During execution, the protocol reads the latent embeddings associated
        with the particles and performs a series of analyses that may include
        dimensionality reduction, clustering, graph traversal, and volume
        generation. The outputs include updated particle sets carrying latent
        coordinates as metadata together with representative reconstructed
        density maps.

        Biological Interpretation of Latent Space

        The latent space represents conformational variability learned directly
        from the experimental particles. Nearby points generally correspond to
        structurally related conformations, whereas distant points indicate more
        pronounced structural differences. The analysis tools provided here are
        intended to help transform these abstract numerical coordinates into
        interpretable biological transitions.

        In many biological systems, latent dimensions may capture motions such
        as domain opening and closing, ligand-dependent rearrangements, or
        continuous assembly dynamics. However, users should interpret latent
        organization carefully. Distances in latent space do not necessarily
        correspond to physical free-energy differences, and poorly sampled
        regions may generate unreliable reconstructions.

        Principal Component Traversals

        The protocol can generate traversals along principal directions of the
        latent space. These trajectories are useful for visualizing dominant
        modes of structural variability across the dataset. Biologically, such
        traversals often reveal gradual transitions between conformations and
        can help identify hinge motions, coordinated domain rearrangements, or
        continuous flexibility.

        In practice, principal component traversals are especially informative
        when the conformational landscape is smooth and well sampled. When the
        dataset contains multiple disconnected states or strong compositional
        heterogeneity, interpretation becomes more complex and additional
        clustering analyses are often beneficial.

        K-Means Sampling and Representative Volumes

        One of the core objectives of the protocol is to generate a manageable
        collection of representative density maps from the latent space. This
        is achieved through clustering strategies that partition the conformational
        landscape into representative regions.

        The resulting volumes provide biologically meaningful snapshots of the
        molecular ensemble. Users commonly inspect these maps visually to detect
        distinct structural states, compare ligand occupancy, identify flexible
        domains, or select subsets for higher-resolution refinement workflows.

        Choosing the number of representative samples requires biological
        judgment. Too few samples may oversimplify the landscape and hide
        important intermediates, while too many may produce highly redundant
        volumes that complicate interpretation.

        Graph Traversal and Continuous Conformational Pathways

        The graph traversal option attempts to identify continuous pathways
        through latent space while remaining within regions supported by the
        experimental data. This approach is particularly valuable for studying
        smooth conformational transitions rather than isolated structural states.

        From a biological perspective, graph traversal can help visualize
        molecular trajectories connecting different conformations. Examples
        include ribosomal rotations, channel gating motions, or domain
        rearrangements in molecular motors. Because the generated paths remain
        constrained by occupied regions of latent space, the resulting
        trajectories are generally more reliable than naive interpolation
        between distant conformations.

        Nevertheless, users should remain cautious when interpreting these
        pathways as true kinetic or energetic transitions. The trajectories
        represent geometrical continuity within the learned embedding space
        rather than experimentally measured reaction coordinates.

        Conformational Landscape Analysis

        The conformational landscape analysis mode provides a more comprehensive
        framework for studying structural heterogeneity. It combines clustering,
        dimensionality reduction, and map generation to produce a structured
        overview of the conformational organization learned during training.

        This analysis is especially useful for large and heterogeneous datasets
        where visual inspection alone becomes difficult. By assigning particles
        into conformational regions, the protocol facilitates downstream focused
        refinement strategies and enables quantitative comparisons between
        states.

        In many practical workflows, this analysis serves as a bridge between
        exploratory heterogeneous reconstruction and high-resolution refinement
        of biologically relevant substates.

        Masking Strategies

        The protocol supports both automatic and user-provided masking during
        landscape analysis. Masking is biologically important because it defines
        which regions of the reconstruction contribute most strongly to the
        analysis and clustering procedures.

        Automatic masking is convenient for exploratory studies and generally
        performs well for compact particles with moderate flexibility. However,
        custom masks are often preferable for complex assemblies containing
        highly mobile domains, detergent micelles, disordered regions, or large
        solvent regions.

        Biologically meaningful masks should focus on structurally conserved
        regions while excluding highly noisy or irrelevant density. Poor masking
        may distort clustering results or artificially emphasize non-biological
        variability.

        Downsampling and Computational Considerations

        Optional volume downsampling can substantially reduce computational
        requirements during exploratory analyses. This is often useful when
        studying very large complexes or when rapidly screening conformational
        variability before committing to high-resolution reconstruction.

        Downsampling reduces memory consumption and accelerates volume
        generation, although it also decreases structural detail. For final
        interpretation or publication-quality analyses, users generally return
        to the original sampling whenever feasible.

        Outputs and Their Interpretation

        The protocol produces an updated particle set containing latent-space
        coordinates associated with each particle. These coordinates may be used
        in downstream workflows for particle selection, clustering, or focused
        refinement.

        In addition, the protocol generates representative volumes sampled from
        the latent space. These maps are intended to summarize the major
        conformational states present in the dataset and provide interpretable
        structural snapshots for biological analysis.

        Depending on the selected analysis options, additional outputs may
        include graph traversal trajectories, dimensionality reduction
        embeddings, clustering assignments, and conformational landscape
        representations.

        Practical Recommendations

        For exploratory analysis, it is often useful to begin with a moderate
        number of representative samples and inspect the resulting volumes
        visually. If the latent organization appears smooth and continuous,
        principal component traversals and graph traversal analyses can provide
        valuable insight into molecular motions.

        For highly heterogeneous datasets, conformational landscape analysis
        combined with carefully designed masks often produces more interpretable
        results. Users should also verify that generated maps correspond to
        physically meaningful conformations rather than noise-driven artifacts.

        When studying subtle structural rearrangements, maintaining the original
        box size and sampling rate is generally preferable. Conversely,
        downsampling can accelerate exploratory workflows during early stages of
        analysis.

        Final Perspective

        The Analyze Results protocol provides a bridge between neural-network
        latent representations and biologically interpretable structural
        variability. Rather than treating heterogeneity as a nuisance, the
        protocol enables researchers to directly explore conformational continua,
        identify representative molecular states, and visualize dynamic
        transitions embedded within cryo-EM datasets.

        For many modern cryo-EM studies, understanding conformational dynamics
        is as important as obtaining high-resolution structures. Careful
        interpretation of latent-space organization, thoughtful masking, and
        biologically informed selection of representative states are essential
        for extracting reliable insights from heterogeneous datasets.
    """

    _label = "analyze results"
    _devStatus = PROD

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """
        out = lambda p: self._getOutputDir(f'analyze.{self._epoch}', p)

        myDict = {
            'input_mask': self._getExtraPath("input_mask.mrc"),
            'output_vol': out('vol_%(id)03d.mrc'),
            'output_volN': out('kmeans%(ksamples)d/vol_%(id)03d.mrc'),
            'z_values': out('z_values.txt'),
            'z_valuesN': out('kmeans%(ksamples)d/z_values.txt'),
            'kmeans_centers': out('kmeans%(ksamples)d/centers_ind.txt'),
            'graph_path': out('graph_traversal/path.txt'),
            'graph_pathZ': out('graph_traversal/z.path.txt'),
            'graph_vols': out('graph_traversal'),
            'umaps': out('umap.pkl')
        }
        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass="SetOfParticles, SetOfParticlesFlex",
                      label='Input Particles')

        form.addParam('inputProt', params.PointerParam, important=True,
                      pointerClass='CryoDrgnProtTrain, CryoDrgnProtAbinitio',
                      label="Previous run to analyse")

        form.addParam('inputEpoch', params.EnumParam,
                      choices=['last', 'selection'], default=EPOCH_LAST,
                      display=params.EnumParam.DISPLAY_LIST,
                      label="Epoch to analyze")

        form.addParam('epochNum', params.IntParam,
                      condition='inputEpoch==%d' % EPOCH_SELECTION,
                      label="Epoch number")

        form.addSection(label='Analysis')
        form.addParam('doGraphTraversal', params.BooleanParam, default=False,
                      label="Do graph traversal?",
                      help="CryoDRGN's graph traversal algorithm builds a nearest "
                           "neighbor graph between all the latent embeddings, and "
                           "then performs Dijkstra's algorithm to find the shortest "
                           "path on the graph between the anchors nodes. The "
                           "idea is to define a trajectory in latent space while "
                           "remaining on the data manifold since we don't want "
                           "to generate structures from unoccupied regions of "
                           "the latent space.")

        form.addParam('doLandscape', params.BooleanParam, default=False,
                      label="Perform conformational landscape analysis?",
                      help="Runs landscape analysis tool for comprehensive and "
                           "quantitative analysis of a trained cryodrgn model, "
                           "including *1) assigning discrete conformational "
                           "states (and providing their particle lists for "
                           "refinement) and 2) visualizing continuous "
                           "conformational landscapes*. This tool also allows "
                           "the user to focus their analysis on specific regions "
                           "of interest by providing custom masks.")

        form.addParam('numVols', params.IntParam, default=500,
                      condition='doLandscape',
                      label="Number of volumes to generate")

        group = form.addGroup('Masking', condition='doLandscape')
        group.addParam('autoMask', params.BooleanParam, default=True,
                       label="Mask volumes automatically?")

        group.addParam('inputMask', params.PointerParam, important=True,
                       condition="not autoMask",
                       pointerClass='VolumeMask', allowsNull=True,
                       label="Custom mask")

        group.addParam('threshold', params.FloatParam, default=0.,
                       condition="autoMask",
                       label="Threshold for masking",
                       help="Default 0 means a half of max density value",
                       expertLevel=params.LEVEL_ADVANCED)

        group.addParam('dilate', params.IntParam, default=5,
                       condition="autoMask",
                       label="Dilation (px)",
                       help="Dilate initial mask by this amount",
                       expertLevel=params.LEVEL_ADVANCED)

        group = form.addGroup('Clustering', condition='doLandscape')
        group.addParam('linkage', params.EnumParam,
                       choices=['average', 'ward'],
                       default=CLUSTER_WARD,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label="Linkage for agglomerative clustering")

        group.addParam('numClusters', params.IntParam, default=10,
                       label="Number of clusters")

        form.addParam('doDownsample', params.BooleanParam, default=False,
                      label="Downsample volumes?")

        form.addParam('boxSize', params.IntParam, default=128,
                      condition='doDownsample', label="New box size (px)")

        form.addParam('pc', params.IntParam, default=2,
                      label="Number of principal components",
                      help="Number of principal component traversals to generate.",
                      expertLevel=params.LEVEL_ADVANCED)

        form.addParam('ksamples', params.IntParam, default=20,
                      label='Number of K-means samples to generate',
                      help="*cryodrgn analyze* uses the k-means clustering "
                           "algorithm to partition the latent space into "
                           "regions (by default k=20 regions), and generate a "
                           "density map from the center of each of these "
                           "regions. The goal is to provide a tractable number "
                           "of representative density maps to visually inspect.",
                      expertLevel=params.LEVEL_ADVANCED)

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
        inputProt = self._getInputProt()
        inputProt._createFilenameTemplates()

        if self.inputEpoch == EPOCH_LAST:
            self._epoch = self._getLastEpoch()
        else:
            self._epoch = self.epochNum.get() - 1

        self._createFilenameTemplates()

        self.weights = self._getInputProt()._getFileName('weights_final')
        self.config = self._getInputProt()._getFileName('config')

        if self.doLandscape and self.hasMultLatentVars():
            self._insertFunctionStep(self.convertInputStep, needsGPU=False)

        self._insertFunctionStep(self.runAnalysisStep, self._epoch, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        if not self.autoMask:
            maskFn = self.inputMask.get().getFileName()
            if pwutils.getExt(maskFn) == ".mrc":
                pwutils.createLink(maskFn, self._getFileName("input_mask"))
            else:
                ih = ImageHandler()
                ih.convert(maskFn, self._getFileName("input_mask"), DT_FLOAT)

    def runAnalysisStep(self, epoch):
        pwutils.makePath(self._getOutputDir())

        self._runProgram('analyze', self._getAnalyzeArgs(epoch))

        pwutils.makePath(self._getOutputDir(f'landscape.{epoch}'))
        pwutils.copyFile(self._getFileName('umaps'),
                         self._getOutputDir(f'landscape.{epoch}/umap.pkl'))

        if self.doGraphTraversal and self.hasMultLatentVars():
            self._runProgram('graph_traversal', self._getGraphArgs())
            self._runProgram('eval_vol', self._getEvalArgs())

        if self.doLandscape and self.hasMultLatentVars():
            self._runProgram('analyze_landscape', self._getLandscapeArgs(epoch))

    def createOutputStep(self):
        # Creating a set of particles with z_values.
        inSet = self._getInputParticles()
        zIterValues = iter(self._getParticlesZvalues())

        outSet = self._createSetOfParticlesFlex(progName=CRYODRGN)
        outSet.copyInfo(inSet)
        outSet.setHasCTF(inSet.hasCTF())
        outSet.getFlexInfo().setProgName(CRYODRGN)

        for particle, zValue in zip(inSet, zIterValues):
            outParticle = emobj.ParticleFlex(progName=CRYODRGN)
            outParticle.copyInfo(particle)
            outParticle.getFlexInfo().setProgName(CRYODRGN)
            outParticle.setZFlex(list(zValue))
            outSet.append(outParticle)

        outSet.getFlexInfo().setAttr(WEIGHTS, String(self.weights))
        outSet.getFlexInfo().setAttr(CONFIG, String(self.config))
        outSet.getFlexInfo().setAttr(ZDIM, Integer(self._getInputProt().zDim))

        self._defineOutputs(outputParticles=outSet)
        self._defineSourceRelation(inSet, outSet)

        # Create a set of k-means sample volumes with z_values.
        fn = self._getExtraPath('volumes.sqlite')
        samplingRate = self._getOutputSampling()
        files, zValues = self._getVolumes()
        setOfVolumes = self._createVolumeSet(files, zValues, fn, samplingRate)

        self._defineOutputs(outputVolumes=setOfVolumes)
        self._defineSourceRelation(self._getInputProt()._getInputParticles(), setOfVolumes)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _warnings(self):
        warnings = []

        if not self.hasMultLatentVars():
            warnings.append("Input protocol has *zDim=1*, the following "
                            "parameters will be ignored:\n"
                            "\t- Number of principal components\n"
                            "\t- Number of K-means samples to generate\n"
                            "\t- Do graph traversal?\n"
                            "\t- Perform conformational landscape analysis?")

        return warnings

    def _validate(self):
        errors = []
        inputProt = self._getInputProt()

        # ab initio homo is not allowed
        if inputProt.getClassName() == "CryoDrgnProtAbinitio":
            run = inputProt.continueRun.get() if inputProt.doContinue else inputProt
            if run.protType.get() == AB_INITIO_HOMO:
                errors.append("Cannot analyze ab initio homogeneous run!")

        if self.inputEpoch == EPOCH_SELECTION:
            inputProt._createFilenameTemplates()
            ep = self.epochNum.get() - 1
            total = self._getLastEpoch()
            if ep > total:
                errors.append(f"You can analyse only epochs 1-{total + 1}")

        if self.doDownsample:
            origBox = self._getBoxSize()
            newBox = self.boxSize.get()
            if newBox > origBox:
                errors.append("You cannot upscale volumes!")

        if self.doLandscape:
            if self.inputMask.hasValue():
                maskSize = self.inputMask.get().getXDim()
                origBox = self._getBoxSize()
                volSize = self.boxSize if self.doDownsample else origBox
                if maskSize != volSize:
                    errors.append("Mask dimensions do not match the output volumes!")

            if not self.autoMask and not self.inputMask.hasValue():
                errors.append("Please provide an input mask or choose auto-masking!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getAnalyzeArgs(self, epoch):
        args = [
            self._getInputProt()._getExtraPath("output"),
            f"{epoch}",
            f"-o {self._getOutputDir(f'analyze.{epoch}')}",
            f"--Apix {self._getSamplingRate()}",
            f"--device {self.gpuList.get()}",
            f"-d {self.boxSize}" if self.doDownsample else "",
            f"--ksample {self.ksamples}" if self.hasMultLatentVars() else "",
            f"--pc {self.pc}" if self.hasMultLatentVars() else ""
        ]
        return args

    def _getGraphArgs(self):
        args = [
            self._getInputProt()._getFileName('z_final'),
            f"--anchors {self._getFileName('kmeans_centers', ksamples=self.ksamples)}",
            f"--outtxt {self._getFileName('graph_pathZ')}",
            f"--outind {self._getFileName('graph_path')}"
        ]
        return args

    def _getEvalArgs(self):
        args = [
            self._getInputProt()._getFileName('weights_final'),
            f"-c {self._getInputProt()._getFileName('config')}",
            f"--zfile {self._getFileName('graph_pathZ')}",
            f"-o {self._getFileName('graph_vols')}"
        ]
        return args

    def _getLandscapeArgs(self, epoch):
        args = [
            self._getInputProt()._getExtraPath("output"),
            f"{epoch}",
            f"-o {self._getOutputDir(f'landscape.{epoch}')}",
            f"--Apix {self._getSamplingRate()}",
            f"--device {self.gpuList.get()}",
            "--skip-umap",
            f"-N {self.numVols}",
            f"--linkage {self.getEnumText('linkage')}",
            f"-M {self.numClusters}",
            f"-d {self.boxSize if self.doDownsample else self._getBoxSize()}",
            f"--pc-dim {min(self.numVols, 20)}"
        ]

        if self.autoMask:
            args.append(f"--dilate {self.dilate}")

            if not self.threshold < 0.001:  # consider as 0
                args.append(f"--thresh {self.threshold}")

        else:
            args.append(f"--mask {self._getFileName('input_mask')}")
        return args

    def _getParticlesZvalues(self):
        """
        Read from z.npz file the particles z_values
        """
        zfile = self._getInputProt()._getFileName('z_final')
        zValues = pickle.load(open(zfile, "rb"))
        return zValues

    def _getVolumes(self):
        """ Returns a list of volume names and their zValues. """
        vols = []
        if self.hasMultLatentVars():
            fn = 'output_volN'
            num = self.ksamples.get()
            zValue = 'z_valuesN'
            zValues = self._getVolumeZvalues(self._getFileName(zValue,
                                                               ksamples=num))
        else:
            fn = 'output_vol'
            num = 10
            zValue = 'z_values'
            zValues = self._getVolumeZvalues(self._getFileName(zValue))

        for volId in range(num):
            if self.hasMultLatentVars():
                volFn = self._getFileName(fn, ksamples=num, epoch=self._epoch,
                                          id=volId + 1)
            else:
                volFn = self._getFileName(fn, epoch=self._epoch, id=volId + 1)

            if os.path.exists(volFn):
                vols.append(volFn)
            else:
                raise FileNotFoundError("Volume %s does not exists. \n"
                                        "Please select a valid epoch "
                                        "number." % volFn)

        return vols, zValues

    def _getVolumeZvalues(self, zValueFile):
        """
        Read from z_values.txt file the volume z_values
        :return: a list with the volumes z_values
        """
        return np.loadtxt(zValueFile, dtype=float).tolist()

    def _createVolumeSet(self, files, zValues, path, samplingRate,
                         updateItemCallback=None):
        """
        Create a set of volume with the associated z_values
        :param files: list of the volumes path
        :param zValues: list with the volumes z_values
        :param path: output path
        :param samplingRate: volumes sampling rate
        :return: a set of volumes
        """
        pwutils.cleanPath(path)
        volSet = emobj.SetOfVolumes(filename=path)
        volSet.setSamplingRate(samplingRate)
        volSet.setObjComment("k-means sample volumes")
        volId = 0
        if type(zValues[0]) is not list:
            # csvList requires each item as a list
            zValues = [[i] for i in zValues]

        for volFn in files:
            vol = emobj.Volume()
            vol.setFileName(volFn)
            vector = CsvList()
            # We assume that each row "i" of z_values corresponds to each
            # volumes with ID "i"
            volZValues = zValues[volId]
            vector._convertValue(volZValues)
            # Creating a new column in the volumes with the z_value
            setattr(vol, Z_VALUES, vector)
            if updateItemCallback:
                updateItemCallback(vol)
            volSet.append(vol)
            volId += 1

        return volSet

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _runProgram(self, program, args):
        gpus = ','.join(str(i) for i in self.getGpuList())
        self.runJob(Plugin.getProgram(program, gpus), ' '.join(args))

    def _getSamplingRate(self):
        return self._getInputProt()._getInputParticles().getSamplingRate()

    def _getBoxSize(self):
        return self._getInputProt()._getInputParticles().getXDim()

    def _getOutputDir(self, *paths):
        return self._getExtraPath("output", *paths)

    def _getOutputSampling(self):
        if self.doDownsample:
            origBox = self._getBoxSize()
            newBox = self.boxSize.get()
            return origBox / newBox * self._getSamplingRate()
        else:
            return self._getSamplingRate()

    def hasMultLatentVars(self):
        inputProt = self._getInputProt()
        if inputProt.doContinue:
            return inputProt.continueRun.get().zDim.get() > 1
        else:
            return inputProt.zDim.get() > 1

    def _getInputProt(self):
        return self.inputProt.get()

    def _getLastEpoch(self):
        outDir = self._getInputProt()._getExtraPath("output")
        files = [file for file in os.listdir(outDir) if file.startswith("weights")]
        print(os.path.basename(files[0]))
        if len(files) != 0:
            epochs = [os.path.basename(file).split('.')[1] for file in files]
            lastEpoch = max([int(epoch) for epoch in epochs if epoch != "pkl"])
        return lastEpoch