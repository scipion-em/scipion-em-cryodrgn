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
from emtools.utils import Timer

from pyworkflow.object import Set
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol

from ..objects import TiltParticle1, SetOfTiltSeriesParticles


class outputs(Enum):
    outputTSP = SetOfTiltSeriesParticles


class CryoDrgnExtractTSParticles(EMProtocol):
    """ Extract 2D series of particles from the tilt series. """
    _label = 'extract tilt-series particles'
    _possibleOutputs = outputs

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam("inputTS", params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label='Tilt series with alignment, non-interpolated',
                      important=True)

        form.addParam("inputCoords", params.PointerParam,
                      label="Coordinates",
                      important=True,
                      pointerClass='SetOfCoordinates3D',
                      help='Input 3D coordinates.')

        form.addParam("inputCTF", params.PointerParam,
                      label="CTF tomo series",
                      pointerClass='SetOfCTFTomoSeries',
                      allowsNull=True,
                      help='Estimated CTF for the tilt series.')

        form.addParam("boxSize", params.IntParam,
                      default=200,
                      label='Box size, unbinned (pix)')

        form.addParam('binFactor', params.IntParam,
                      default=4,
                      label='Binning factor')

    # -------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.runStep)
        self._insertFunctionStep(self.closeStep)

    # --------------------------- STEPS functions ------------------------------
    def runStep(self):
        # Load input TS
        t = Timer()
        inputTS = self._getInputTS()
        tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in inputTS}
        acq = self._getInputTS().getAcquisition()
        pix = self._getOutputSampling()
        inputTS.close()
        t.toc(f'Loading TS completed')

        # Load input CTF
        t = Timer()
        inputCTF = self.inputCTF.get()
        ctfDict = {ctf.getTsId(): ctf.clone(ignoreAttrs=[]) for ctf in inputCTF}
        inputCTF.close()
        t.toc(f'Loading CTF completed')

        # Create output set
        out = self.getOutput()
        out.setAcquisition(acq)
        out.setSamplingRate(pix)

        for volId, coord in enumerate(self.inputCoords.get().iterCoordinates()):
            t = Timer()
            tsId = coord.getTomoId()
            #matrix = coord.getMatrix()
            ts = tsDict[tsId]
            ctfseries = ctfDict.get(tsId)

            i = 1
            for tiltImage, ctf in zip(ts.iterItems(orderBy="_index"),
                                      ctfseries.iterItems(orderBy="_index")):
                tp = TiltParticle1()
                tp.setLocation((i, f"{tsId}_{volId}_stack.mrcs"))
                tp.setSamplingRate(pix)

                acq.setAccumDose(tiltImage.getAcquisition().getAccumDose())
                tp.setAcquisition(acq)
                tp.setAcquisitionOrder(tiltImage.getAcquisitionOrder())
                tp.setTiltAngle(tiltImage.getTiltAngle())

                tp.setTsId(tsId)
                tp.setVolId(volId)

                tp.setCoordinate(coord)
                tp.setTransform(tiltImage.getTransform())
                tp.setCTF(ctf)

                out.append(tp)

            t.toc(f'Iter over coord (loop {volId+1}) completed')

        out.write()
        self._store(out)

    def closeStep(self):
        self.getOutput().setStreamState(Set.STREAM_CLOSED)
        self._store()

    def getOutput(self):
        output = outputs.outputTSP.name
        if hasattr(self, output):
            getattr(self, output).enableAppend()
        else:
            outSet = SetOfTiltSeriesParticles.create(self._getPath())
            outSet.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outputs.outputTSP.name: outSet})
            self._defineSourceRelation(self._getInputTS(pointer=True), outSet)

        return getattr(self, output)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        validate = []

        ts = self._getInputTS()
        if not ts.getFirstItem().getFirstItem().hasTransform():
            validate.append("Input tilt-series are missing "
                            "a transformation matrix.")

        ctf = self.inputCTF.get()
        if len(ts) != len(ctf):
            validate.append("Number of input TS and CTFs does not match.")

        return validate

    def _summary(self):
        return []

    # --------------------------- UTILS functions -----------------------------
    def _getInputTS(self, pointer=False):
        return self.inputTS if pointer else self.inputTS.get()

    def _getOutputSampling(self):
        return self._getInputTS().getSamplingRate() * self.binFactor.get()
