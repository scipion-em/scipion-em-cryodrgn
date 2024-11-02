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

import os
import numpy as np
from enum import Enum

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW
import pyworkflow.object as pwobj
from pwem.protocols import ProtAnalysis3D
from pwem.objects import SetOfVolumes, Volume
from pwem.emlib.image import ImageHandler, DT_FLOAT

from cryodrgn import Plugin
from cryodrgn.constants import (EPOCH_LAST, EPOCH_SELECTION, Z_VALUES,
                                AB_INITIO_HOMO, CLUSTER_WARD, V3_3_2)
from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class outputs(Enum):
    Volumes = SetOfVolumes


class CryoDrgnProtAnalyze(ProtAnalysis3D, CryoDrgnProtBase):
    """ CryoDrgn protocol to visualize latent space and generate volumes. """

    _label = "analyze results"
    _devStatus = NEW
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """
        out = lambda p: self.getOutputDir(f'analyze.{self._epoch}', p)

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
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="GPU may have several cores. Set it to zero"
                            " if you do not know what we are talking about."
                            " First core index is 0, second 1 and so on."
                            " You can use only a single GPU.")
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

        form.addSection('Volume generation')
        form.addParam('doFlip', params.BooleanParam, default=False,
                      label="Flip handedness of output volumes")
        form.addParam('doInvert', params.BooleanParam, default=False,
                      label="Invert contrast of output volumes")
        form.addParam('doDownsample', params.BooleanParam, default=False,
                      label="Downsample volumes?")
        form.addParam('boxSize', params.IntParam, default=128,
                      condition='doDownsample', label="New box size (px)")

        form.addSection(label='Latent space')
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
        form.addParam('pc', params.IntParam, default=2,
                      label="Number of principal components",
                      help="Number of principal component traversals to generate.")
        form.addParam('ksamples', params.IntParam, default=20,
                      label='Number of K-means samples to generate',
                      help="*cryodrgn analyze* uses the k-means clustering "
                           "algorithm to partition the latent space into "
                           "regions (by default k=20 regions), and generate a "
                           "density map from the center of each of these "
                           "regions. The goal is to provide a tractable number "
                           "of representative density maps to visually inspect.")

        form.addSection(label='Landscape analysis')
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
        group = form.addGroup('Masking',
                              condition='doLandscape')
        group.addParam('autoMask', params.BooleanParam, default=True,
                       label="Mask volumes automatically?")
        group.addParam('inputMask', params.PointerParam, important=True,
                       condition="not autoMask",
                       pointerClass='VolumeMask', allowsNull=True,
                       label="Custom mask")
        group.addParam('threshold', params.FloatParam, default=0.,
                       condition="autoMask",
                       label="Threshold for masking",
                       help="Default 0 means a half of max density value")
        group.addParam('dilate', params.IntParam, default=5,
                       condition="autoMask",
                       label="Dilation (px)",
                       help="Dilate initial mask by this amount")

        group = form.addGroup('Clustering',
                              condition='doLandscape')
        group.addParam('linkage', params.EnumParam,
                       choices=['average', 'ward'],
                       default=CLUSTER_WARD,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label="Linkage for agglomerative clustering")
        group.addParam('numClusters', params.IntParam, default=10,
                       label="Number of clusters")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        inputProt = self._getInputProt()
        inputProt._createFilenameTemplates()

        if self.inputEpoch == EPOCH_LAST:
            self._epoch = self._getLastEpoch()
        else:
            self._epoch = self.epochNum.get() - 1

        self._createFilenameTemplates()

        self._insertFunctionStep(self.runAnalysisStep, self._epoch)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def runAnalysisStep(self, epoch):
        pwutils.makePath(self.getOutputDir())
        self._runProgram('analyze', self._getAnalyzeArgs(epoch))

        if self.doGraphTraversal and self.hasMultLatentVars():
            self._runProgram('graph_traversal', self._getGraphArgs())
            self._runProgram('eval_vol', self._getEvalArgs())

        if self.doLandscape and self.hasMultLatentVars():
            self.convertInputs(epoch)
            self._runProgram('analyze_landscape', self._getLandscapeArgs(epoch))

    def createOutputStep(self):
        """ Create a set of k-means sample volumes with z_values. """
        fn = self._getExtraPath('volumes.sqlite')
        samplingRate = self._getOutputSampling()
        files, zValues = self._getVolumes()
        setOfVolumes = self._createVolumeSet(files, zValues, fn, samplingRate)

        self._defineOutputs(**{outputs.Volumes.name: setOfVolumes})
        self._defineSourceRelation(self._getInputProt()._getInputParticles(pointer=True),
                                   setOfVolumes)

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
                errors.append(f"You can analyse only epochs 1-{total+1}")

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
    def convertInputs(self, epoch):
        # Copy analyze.epoch/umap.pkl to landscape.epoch folder
        pwutils.makePath(self.getOutputDir(f'landscape.{epoch}'))
        pwutils.copyFile(self._getFileName('umaps'),
                         self.getOutputDir(f'landscape.{epoch}/umap.pkl'))

        if not self.autoMask:
            # convert mask to mrc
            maskFn = self.inputMask.get().getFileName()
            if pwutils.getExt(maskFn) == ".mrc":
                pwutils.createLink(maskFn, self._getFileName("input_mask"))
            else:
                ih = ImageHandler()
                ih.convert(maskFn, self._getFileName("input_mask"), DT_FLOAT)

    def _getAnalyzeArgs(self, epoch):
        args = [
            self._getInputProt()._getExtraPath("output"),
            f"{epoch}",
            f"-o {self.getOutputDir(f'analyze.{epoch}')}",
            f"--Apix {self._getSamplingRate()}",
            f"--device {self.gpuList.get()}",
            f"-d {self.boxSize}" if self.doDownsample else "",
            "--flip" if self.doFlip else "",
            "--invert" if self.doInvert else "",
            f"--ksample {self.ksamples}" if self.hasMultLatentVars() else "",
            f"--pc {self.pc}" if self.hasMultLatentVars() else ""
        ]

        return args

    def _getGraphArgs(self):
        args = [
            self._getInputProt()._getFileName('z_final'),
            f"--anchors $(cat {self._getFileName('kmeans_centers', ksamples=self.ksamples)})"
        ]

        if Plugin.versionGE(V3_3_2):
            args.extend([
                f"--outtxt {self._getFileName('graph_pathZ')}",
                f"--outind {self._getFileName('graph_path')}"
            ])
        else:
            args.extend([
                f"-o {self._getFileName('graph_path')}",
                f"--out-z {self._getFileName('graph_pathZ')}"
            ])

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
            f"-o {self.getOutputDir(f'landscape.{epoch}')}",
            f"--Apix {self._getSamplingRate()}",
            f"--device {self.gpuList.get()}",
            "--skip-umap",
            f"-N {self.numVols}",
            f"--linkage {self.getEnumText('linkage')}",
            f"-M {self.numClusters}",
            f"-d {self.boxSize if self.doDownsample else self._getBoxSize()}",
            "--flip" if self.doFlip else "",
            f"--pc-dim {min(self.numVols, 20)}"
        ]

        if self.autoMask:
            args.append(f"--dilate {self.dilate}")

            if not self.threshold < 0.001:  # consider as 0
                args.append(f"--thresh {self.threshold}")

        else:
            args.append(f"--mask {self._getFileName('input_mask')}")

        return args

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
                                          id=volId)
            else:
                volFn = self._getFileName(fn, epoch=self._epoch, id=volId)

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
        volSet = SetOfVolumes(filename=path)
        volSet.setSamplingRate(samplingRate)
        volSet.setObjComment("k-means sample volumes")
        volId = 0
        if type(zValues[0]) is not list:
            # csvList requires each item as a list
            zValues = [[i] for i in zValues]

        for volFn in files:
            vol = Volume()
            vol.setFileName(volFn)
            vector = pwobj.CsvList()
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

    def _getSamplingRate(self):
        return self._getInputProt()._getInputParticles().getSamplingRate()

    def _getBoxSize(self):
        return self._getInputProt()._getInputParticles().getXDim()

    def _getOutputSampling(self):
        if self.doDownsample:
            origBox = self._getBoxSize()
            newBox = self.boxSize.get()
            return origBox/newBox * self._getSamplingRate()
        else:
            return self._getSamplingRate()

    def hasMultLatentVars(self):
        inputProt = self._getInputProt()
        if inputProt.doContinue:
            return inputProt.continueRun.get().zDim.get() > 1
        else:
            return inputProt.zDim.get() > 1

    def _getLastEpoch(self):
        return self._getInputProt()._getLastEpoch()

    def _getInputProt(self):
        return self.inputProt.get()
