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
import pickle
import numpy as np
import re
from glob import glob

import pyworkflow.utils as pwutils
from pyworkflow.constants import NEW
import pyworkflow.protocol.params as params
import pyworkflow.object as pwobj
from pwem.protocols import ProtProcessParticles
import pwem.objects as emobj

from cryodrgn import Plugin
from cryodrgn.constants import *


class CryoDrgnProtAbinitio(ProtProcessParticles):
    """
    Protocol to run ab-initio reconstruction with cryoDRGN2 neural network.
    """
    _label = 'training ab initio'
    _devStatus = NEW

    @classmethod
    def isDisabled(cls):
        return not Plugin.versionGE(V2_0_0)

    def _createFilenameTemplates(self):
        """ Centralize how files are called within the protocol. """

        def out(*p):
            return os.path.join(self.getOutputDir(f'analyze.{self._epoch}'), *p)

        myDict = {
            'output_vol': out('vol_%(id)03d.mrc'),
            'output_volN': out('kmeans%(ksamples)d/vol_%(id)03d.mrc'),
            'z_values': out('z_values.txt'),
            'z_valuesN': out('kmeans%(ksamples)d/z_values.txt'),
            'weights': self.getOutputDir(f'weights.{self._epoch}.pkl'),
            'config': self.getOutputDir('config.pkl')
        }
        self._updateFilenamesDict(myDict)

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

        form.addParam('protType', params.EnumParam,
                      choices=['homogeneous', 'heterogeneous'],
                      default=AB_INITIO_HOMO,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Ab initio type')

        form.addParam('zDim', params.IntParam, default=1,
                      validators=[params.Positive],
                      label='Dimension of latent variable',
                      help='It is recommended to first train on lower '
                           'resolution images (e.g. D=128) with '
                           '--zdim 1 and with --zdim 10 using the '
                           'default architecture (fast).')
        form.addParam('numEpochs', params.IntParam, default=30,
                      label='Number of epochs',
                      help='The number of epochs refers to the number '
                           'of full passes through the dataset for '
                           'training. For reference, ab initio heterogeneous '
                           'reconstruction of a recent dataset of ~218k '
                           '128x128 particles took 20 hours to train '
                           '(1 V100 GPU).')

        form.addSection(label='Advanced')
        form.addParam('searchRange', params.IntParam, default=10,
                      label='Translational search range (px)')
        form.addParam('psFreq', params.IntParam, default=5,
                      label='Update poses every N epochs')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn abinit_homo -h* for help.")

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
        self._insertFunctionStep(self.runTrainingStep)
        if self.viewEpoch == EPOCH_LAST:
            self._epoch = self.numEpochs.get() - 1
        else:
            self._epoch = self.epochNum.get()
        self._createFilenameTemplates()
        self._insertFunctionStep(self.runAnalysisStep, self._epoch)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def runTrainingStep(self):
        # Create output folder
        pwutils.cleanPath(self.getOutputDir())
        pwutils.makePath(self.getOutputDir())

        # Call cryoDRGN with the appropriate parameters
        protType = self.protType.get()
        program = "homo" if protType == AB_INITIO_HOMO else "het"
        self._runProgram("abinit_" + program, self._getTrainingArgs(protType))

    def runAnalysisStep(self, epoch):
        """ Run analysis step.
        Args:
            epoch: epoch number to be analyzed.
        """
        self._runProgram('analyze', self._getAnalyzeArgs(epoch))

    def createOutputStep(self):
        """ Create the protocol outputs. """
        # Creating a set of particles with z_values
        outImgSet = self._createParticleSet()
        self._defineOutputs(Particles=outImgSet)

        # Creating a set of volumes with z_values
        fn = self._getExtraPath('volumes.sqlite')
        samplingRate = self.inputParticles.get().getSamplingRate()
        files, zValues = self._getVolumes()
        setOfVolumes = self._createVolumeSet(files, zValues, fn, samplingRate)
        self._defineOutputs(Volumes=setOfVolumes)
        self._defineSourceRelation(self.inputParticles.get(), setOfVolumes)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["Training ab initio for %d epochs." % self.numEpochs]

        return summary

    def _validate(self):
        errors = []

        if self.zDim > 1 and self.protType.get() == AB_INITIO_HOMO:
            errors.append("Latent variable must be 1 for "
                          "homogeneous reconstruction")

        if self.viewEpoch == EPOCH_SELECTION:
            ep = self.epochNum.get()
            total = self.numEpochs.get()
            if ep > total:
                errors.append(f"You can analyse only epochs 1-{total}")

        return errors

    def _citations(self):
        return ['Zhong2021b']

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self, protType=AB_INITIO_HOMO):
        parts = self.inputParticles.get()

        args = [
            parts.filename.get(),
            '--ctf %s' % parts.ctfs,
            '-o %s ' % self.getOutputDir(),
            '-n %d' % self.numEpochs,
            '--t-extent %d' % self.searchRange,
            '--ps-freq %d' % self.psFreq,
        ]

        if protType == AB_INITIO_HETERO:
            args.extend([
                '--zdim %d' % self.zDim,
                '--preprocessed',
                '--max-threads %d' % self.numberOfThreads,
            ])

            if len(self.getGpuList()) > 1:
                args.append('--multigpu')

        if self.extraParams.hasValue():
            args.append(self.extraParams.get())

        return args

    def _getAnalyzeArgs(self, epoch):
        return [
            self.getOutputDir(),
            f"{epoch}",
            '--Apix %0.3f' % self.inputParticles.get().getSamplingRate(),
            '--ksample %d' % self.ksamples,
        ]

    def _runProgram(self, program, args):
        gpus = ','.join(str(i) for i in self.getGpuList())
        self.runJob(Plugin.getProgram(program, gpus), ' '.join(args))

    def _getVolumes(self):
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

    def _getParticlesZvalues(self):
        """
        Read from z.pkl file the particles z_values
        :return: a numpy array with the particles z_values
        """
        zEpochFile = self.getEpochZFile(self._epoch)
        with open(zEpochFile, 'rb') as f:
            zValues = pickle.load(f)
        return zValues

    def _createParticleSet(self):
        """
        Create a set of particles with the associated z_values
        :return: a set of particles
        """
        cryoDRGParticles = self.inputParticles.get()
        ImgSet = self.getProject().getProtocol(cryoDRGParticles.getObjParentId()).inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(ImgSet)
        outImgSet.copyItems(ImgSet, updateItemCallback=self._setZValues)

        setattr(outImgSet, WEIGHTS, pwobj.String(self._getFileName('weights')))
        setattr(outImgSet, CONFIG, pwobj.String(self._getFileName('config')))

        return outImgSet

    def _setZValues(self, item, row=None):
        zValues = self._getParticlesZvalues()
        vector = pwobj.CsvList()
        # We assume that each row "i" of z_values corresponds to each
        # particle with ID "i"
        vector._convertValue(list(zValues[item.getObjId()-1]))
        setattr(item, Z_VALUES, vector)

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
        :param zValues: array with the volumes z_values
        :param path: output path
        :param samplingRate: volumes sampling rate
        :return: a set of volumes
        """
        pwutils.cleanPath(path)
        volSet = emobj.SetOfVolumes(filename=path)
        volSet.setSamplingRate(samplingRate)
        volId = 0

        for volFn in files:
            vol = emobj.Volume()
            vol.setFileName(volFn)
            vector = pwobj.CsvList()
            # We assume that each row "i" of z_values corresponds to each
            # volumes with ID "i"
            vector._convertValue(zValues[volId])
            # Creating a new column in the volumes with the z_value
            setattr(vol, Z_VALUES, vector)
            if updateItemCallback:
                updateItemCallback(vol)
            volSet.append(vol)
            volId += 1
        volSet.write()
        volSet.close()

        return volSet

    def getOutputDir(self, *paths):
        return os.path.join(self._getPath('output'), *paths)

    def getEpochZFile(self, epoch):
        return self.getOutputDir('z.%d.pkl' % epoch)

    def getLastEpoch(self):
        """ Return the last iteration number. """
        epoch = None
        self._epochRegex = re.compile('z.(\d).pkl')
        files = sorted(glob(self.getEpochZFile(0).replace('0', '*')))
        if files:
            f = files[-1]
            s = self._epochRegex.search(f)
            if s:
                epoch = int(s.group(1))  # group 1 is a digit iteration number

        return epoch

    def hasMultLatentVars(self):
        return self.zDim > 1