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
import pyworkflow.object as pwobj
from pwem.protocols import ProtAnalysis3D
from pwem.objects import SetOfVolumes, Volume

from cryodrgn.constants import EPOCH_LAST, EPOCH_SELECTION, Z_VALUES
from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class outputs(Enum):
    Volumes = SetOfVolumes


class CryoDrgnProtAnalyze(ProtAnalysis3D, CryoDrgnProtBase):
    """ CryoDrgn protocol to visualize latent space and generate volumes or
        to perform conformational landscape analysis. """

    _label = "analyze results"
    _possibleOutputs = outputs

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """
        def out(*p):
            return os.path.join(self.getOutputDir(f'analyze.{self._epoch}'), *p)

        myDict = {
            'output_vol': out('vol_%(id)03d.mrc'),
            'output_volN': out('kmeans%(ksamples)d/vol_%(id)03d.mrc'),
            'z_values': out('z_values.txt'),
            'z_valuesN': out('kmeans%(ksamples)d/z_values.txt'),
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

        form.addSection(label='Latent space')
        form.addParam('skipVol', params.BooleanParam, default=False,
                      label="Skip generation of volumes")
        form.addParam('skipUmap', params.BooleanParam, default=False,
                      label="Skip running UMAP")

        group = form.addGroup('Volume generation')
        group.addParam('doFlip', params.BooleanParam, default=False,
                       label="Flip handedness of output volumes")
        group.addParam('doInvert', params.BooleanParam, default=False,
                       label="Invert contrast of output volumes")
        group.addParam('doDownsample', params.BooleanParam, default=False,
                       label="Downsample volumes?")
        group.addParam('boxSize', params.IntParam, default=128,
                       condition='doDownsample', label="New box size (px)")
        group.addParam('pc', params.IntParam, default=2,
                       label="Number of principal components",
                       help="Number of principal component traversals to generate.")
        group.addParam('ksamples', params.IntParam, default=20,
                       label='Number of K-means samples to generate',
                       help="*cryodrgn analyze* uses the k-means clustering "
                            "algorithm to partition the latent space into "
                            "regions (by default k=20 regions), and generate a "
                            "density map from the center of each of these "
                            "regions. The goal is to provide a tractable number "
                            "of representative density maps to visually inspect.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        inputProt = self.inputProt.get()
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

    def createOutputStep(self):
        """ Create a set of volumes with z_values. """
        fn = self._getExtraPath('volumes.sqlite')
        samplingRate = self._getSamplingRate()
        files, zValues = self._getVolumes()
        setOfVolumes = self._createVolumeSet(files, zValues, fn, samplingRate)

        self._defineOutputs(**{outputs.Volumes.name: setOfVolumes})
        self._defineSourceRelation(self.inputProt.get()._getInputParticles(pointer=True),
                                   setOfVolumes)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []

        return summary

    def _validate(self):
        errors = []

        if self.inputEpoch == EPOCH_SELECTION:
            self.inputProt.get()._createFilenameTemplates()
            ep = self.epochNum.get() - 1
            total = self._getLastEpoch()
            if ep > total:
                errors.append(f"You can analyse only epochs 1-{total+1}")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getAnalyzeArgs(self, epoch):
        args = [
            self.inputProt.get()._getExtraPath("output"),
            f"{epoch}",
            f"-o {self.getOutputDir(f'analyze.{epoch}')}",
            f"--Apix {self._getSamplingRate()}",
            "--skip-vol" if self.skipVol else "",
            "--skip-umap" if self.skipVol else "",
            f"-d {self.boxSize}" if self.doDownsample else "",
            "--flip" if self.doFlip else "",
            "--invert" if self.doInvert else "",
            f"--pc {self.pc}",
            f"--ksample {self.ksamples}",
            f"--device {self.gpuList.get()}"
        ]

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
        return self.inputProt.get()._getInputParticles().getSamplingRate()

    def hasMultLatentVars(self):
        inputProt = self.inputProt.get()
        if inputProt.doContinue:
            return inputProt.continueRun.get().zDim.get() > 1
        else:
            return inputProt.zDim.get() > 1

    def _getLastEpoch(self):
        return self.inputProt.get()._getLastEpoch()
