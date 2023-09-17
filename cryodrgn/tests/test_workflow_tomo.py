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
from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr
from pwem.tests.workflows.test_workflow import TestWorkflow

from tomo.protocols import ProtImportTs
from reliontomo.protocols import ProtImportCoordinates3DFromStar
from imod.protocols import (ProtImodAutomaticCtfEstimation,
                            ProtImodImportTransformationMatrix,
                            ProtImodTSNormalization,
                            ProtImodTomoReconstruction)

from . import DataSet
from ..protocols import CryoDrgnExtractTSParticles


class TestTomoDrgn(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.tsPath = cls.dataset.getFile("empiar-ts")
        cls.coordsFn = cls.dataset.getFile("empiar-coords")
        cls.pixSize = 1.35
        cls.tiltAxis = 85.3
        cls.dose = 3.05
        cls.boxSize = 192
        cls.defocus = 1900

    def _runImportTiltSeries(self, pixSize, tiltAxis, dose):
        print(magentaStr("\n==> Importing data - tilt series:"))
        protImport = self.newProtocol(
            ProtImportTs,
            filesPath=self.tsPath,
            filesPattern="*.mdoc",
            stepAngle=3.0,
            voltage=300.0,
            sphericalAberration=2.7,
            amplitudeContrast=0.07,
            magnification=105000,
            samplingRate=pixSize,
            tiltAxisAngle=tiltAxis,
            dosePerFrame=dose
            )
        self.launchProtocol(protImport)
        outputTS = getattr(protImport, 'outputTiltSeries')
        self.assertIsNotNone(outputTS, 'No tilt series were imported')
        self.assertSetSize(outputTS, 2)

        return outputTS

    def _runImportCoords(self, starFile, pixSize, boxSize, inTomos=None):
        print(magentaStr("\n==> Importing data - coords 3D:"))
        protImport = self.newProtocol(ProtImportCoordinates3DFromStar,
                                      starFile=starFile,
                                      inTomos=inTomos,
                                      samplingRate=pixSize,
                                      boxSize=boxSize)

        self.launchProtocol(protImport)
        outCoords = getattr(protImport, protImport._possibleOutputs.coordinates.name)
        self.assertSetSize(outCoords, 5693)
        self.assertIsNotNone(outCoords, 'No coordinates were imported')

        return outCoords

    def _runImportTransformationMatrix(self, filesPath, pattern, inputTS=None):
        print(magentaStr("\n==> Importing data - TS alignment:"))
        protImport = self.newProtocol(ProtImodImportTransformationMatrix,
                                      filesPath=filesPath,
                                      filesPattern=pattern,
                                      inputSetOfTiltSeries=inputTS)
        self.launchProtocol(protImport)
        outTS = protImport.TiltSeries
        self.assertIsNotNone(outTS, "No alignment was imported")
        self.assertSetSize(outTS, size=2)

        for ts in outTS:
            self.assertTrue(ts.getFirstItem().hasTransform())

        return outTS

    def _runImodPreprocess(self, inputTS=None, binning=1):
        print(magentaStr("\n==> Running imod - preprocess TS:"))
        prot = self.newProtocol(ProtImodTSNormalization,
                                inputSetOfTiltSeries=inputTS,
                                binning=binning)
        self.launchProtocol(prot)
        outTS = getattr(prot, "TiltSeries")
        self.assertIsNotNone(outTS, "TS binning has failed")

        inSamplingRate = prot.inputSetOfTiltSeries.get().getSamplingRate()
        outSamplingRate = outTS.getSamplingRate()
        self.assertEqual(inSamplingRate * binning, outSamplingRate)

        return outTS

    def _runImodReconstruct(self, inputTS, thickness):
        print(magentaStr("\n==> Running imod - reconstruct tomo:"))
        prot = self.newProtocol(ProtImodTomoReconstruction,
                                inputSetOfTiltSeries=inputTS,
                                tomoThickness=thickness,
                                useGpu=False)
        self.launchProtocol(prot)
        outTomos = getattr(prot, "Tomograms")
        self.assertIsNotNone(outTomos, "SetOfTomograms has not been produced.")

        return outTomos

    def _runCTFEstimation(self, inputTS, defocus):
        print(magentaStr("\n==> Running imod - CTF estimation:"))
        prot = self.newProtocol(ProtImodAutomaticCtfEstimation,
                                inputSet=inputTS,
                                expectedDefocusValue=defocus,
                                angleStep=0,
                                searchAstigmatism=0)
        self.launchProtocol(prot)
        outCTF = getattr(prot, "CTFTomoSeries")
        self.assertIsNotNone(outCTF, 'CTFTomoSeries has not been produced')

        return outCTF

    def test_tomodrgn(self):
        importedTS = self._runImportTiltSeries(self.pixSize,
                                               self.tiltAxis,
                                               self.dose)
        ctf = self._runCTFEstimation(importedTS, self.defocus)
        alignedTS = self._runImportTransformationMatrix(self.tsPath,
                                                        "*.xf",
                                                        inputTS=importedTS)
        binnedTS = self._runImodPreprocess(alignedTS, binning=8)
        tomograms = self._runImodReconstruct(inputTS=binnedTS,
                                             thickness=162)
        importedCoords = self._runImportCoords(self.coordsFn,
                                               self.pixSize,
                                               self.boxSize//8,
                                               tomograms)

        print(magentaStr("\n==> Running tomodrgn - extract tilt series particles:"))
        protExtract = self.newProtocol(CryoDrgnExtractTSParticles,
                                       inputTS=alignedTS,
                                       inputCoords=importedCoords,
                                       inputCTF=ctf,
                                       boxSize=self.boxSize,
                                       binFactor=4)

        self.launchProtocol(protExtract)
        outTomos = getattr(protExtract, "outputTSP")
        self.assertIsNotNone(outTomos, "SetOfTiltSeriesParticles has not been produced.")
