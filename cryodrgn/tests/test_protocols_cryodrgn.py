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

from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportParticles
from pwem.tests.workflows import TestWorkflow

from cryodrgn.protocols import (CryoDrgnProtPreprocess, CryoDrgnProtTrain,
                                CryoDrgnProtAbinitio, CryoDrgnProtAnalyze)


class TestWorkflowCryoDrgn(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/refine3d_case2/relion_data.star')

    def _importParticles(self, partStar, mag, samplingRate):
        """ Import particles from Relion star file. """
        print(magentaStr("\n==> Importing data - particles from star:"))
        protImport = self.newProtocol(ProtImportParticles,
                                      importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                      starFile=partStar,
                                      magnification=mag,
                                      samplingRate=samplingRate,
                                      haveDataBeenPhaseFlipped=False)

        return self.launchProtocol(protImport)

    def _runPreprocess(self, protImport, protLabel, **kwargs):
        print(magentaStr(f"\n==> Testing cryoDRGN - {protLabel}:"))
        protPreprocess = self.newProtocol(CryoDrgnProtPreprocess,
                                          objLabel=protLabel, **kwargs)
        protPreprocess.inputParticles.set(protImport.outputParticles)

        return self.launchProtocol(protPreprocess)

    def _runTraining(self, protPreprocess, **kwargs):
        print(magentaStr("\n==> Testing cryoDRGN - training vae:"))
        protTrain = self.newProtocol(CryoDrgnProtTrain, **kwargs)
        protTrain.inputParticles.set(protPreprocess.Particles)

        return self.launchProtocol(protTrain)

    def _runAbinitio(self, protPreprocess, protLabel, **kwargs):
        print(magentaStr(f"\n==> Testing cryoDRGN - ab initio ({protLabel}):"))
        protAbinitio = self.newProtocol(CryoDrgnProtAbinitio,
                                        objLabel=protLabel, **kwargs)
        protAbinitio.inputParticles.set(protPreprocess.Particles)

        return self.launchProtocol(protAbinitio)

    def _runAnalyze(self, protTrain, **kwargs):
        print(magentaStr("\n==> Testing cryoDRGN - analyze results:"))
        protAnalyze = self.newProtocol(CryoDrgnProtAnalyze, **kwargs)
        protAnalyze.inputProt.set(protTrain)

        return self.launchProtocol(protAnalyze)

    def testWorkflow(self):
        protImport = self._importParticles(self.partFn, 50000, 3.54)

        protPreprocess1 = self._runPreprocess(protImport,
                                              "downsample scale=64",
                                              scaleSize=64)
        self.assertIsNotNone(protPreprocess1._possibleOutputs.Particles.name)

        protPreprocess2 = self._runPreprocess(protImport,
                                              "downsample scale=48 with chunks",
                                              scaleSize=48, chunk=200)
        self.assertIsNotNone(protPreprocess2._possibleOutputs.Particles.name)

        protTraining = self._runTraining(protPreprocess2, numEpochs=3, zDim=2)
        self.assertIsNotNone(protTraining._possibleOutputs.Particles.name)

        protAbinitio = self._runAbinitio(protPreprocess2, "hetero",
                                         numEpochs=2, zDim=2)
        self.assertIsNotNone(protAbinitio._possibleOutputs.Particles.name)

        protAbinitio2 = self._runAbinitio(protPreprocess2, "homo",
                                          numEpochs=2, zDim=1, protType=0)
        self.assertIsNotNone(protAbinitio2._possibleOutputs.Volumes.name)

        protAnalyze = self._runAnalyze(protTraining)
        self.assertIsNotNone(protAnalyze._possibleOutputs.Volumes.name)
