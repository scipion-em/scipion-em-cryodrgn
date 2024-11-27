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

import pyworkflow.protocol.params as params
from pyworkflow.constants import NEW

from cryodrgn.constants import CRYODRGN
from cryodrgn.protocols.protocol_base import CryoDrgnProtBase


class CryoDrgnProtSubset(CryoDrgnProtBase):
    """ CryoDrgn protocol to make a particles subset using a pkl file. """

    _label = "particles subset"
    _devStatus = NEW
    _possibleOutputs = CryoDrgnProtBase._possibleOutputs
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
        self._insertFunctionStep(self.createOutputStep)

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
