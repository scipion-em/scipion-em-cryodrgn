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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportParticles

from cryodrgn.protocols import (CryoDrgnProtPreprocess, CryoDrgnProtTrain,
                                CryoDrgnProtAbinitio)


class TestCryoDrgn(BaseTest):
    @classmethod
    def runImportParticlesStar(cls, partStar, mag, samplingRate):
        """ Import particles from Relion star file. """
        print(magentaStr("\n==> Importing data - particles from star:"))
        protImport = cls.newProtocol(ProtImportParticles,
                                     importFrom=ProtImportParticles.IMPORT_FROM_RELION,
                                     starFile=partStar,
                                     magnification=mag,
                                     samplingRate=samplingRate,
                                     haveDataBeenPhaseFlipped=False)
        cls.launchProtocol(protImport)
        cls.assertIsNotNone(protImport.outputParticles,
                            "SetOfParticles has not been produced.")

        return protImport

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('relion_tutorial')
        cls.partFn = cls.dataset.getFile('import/refine3d_case2/relion_data.star')
        cls.protImport = cls.runImportParticlesStar(cls.partFn, 50000, 3.54)

    def runPreprocess(self, protLabel, particles, **kwargs):
        print(magentaStr(f"\n==> Testing cryoDRGN - {protLabel}:"))
        protPreprocess = self.newProtocol(CryoDrgnProtPreprocess,
                                          objLabel=protLabel, **kwargs)
        protPreprocess.inputParticles.set(particles)
        return self.launchProtocol(protPreprocess)

    def checkPreprocessOutput(self, preprocessProt):
        output = getattr(preprocessProt,
                         preprocessProt._possibleOutputs.Particles.name,
                         None)
        self.assertIsNotNone(output)

    def checkTrainOutput(self, trainProt):
        output = getattr(trainProt, trainProt._possibleOutputs.Particles.name, None)
        self.assertIsNotNone(output)

    def testPreprocess(self):
        parts = self.protImport.outputParticles

        preprocess1 = self.runPreprocess("downsample scale=64", parts, scaleSize=64)
        self.checkPreprocessOutput(preprocess1)

        preprocess2 = self.runPreprocess("downsample scale=50 with chunks", parts,
                                         scaleSize=50, chunk=200)
        self.checkPreprocessOutput(preprocess2)

    def testTraining(self):
        parts = self.protImport.outputParticles
        preprocess = self.runPreprocess("downsample scale=64", parts, scaleSize=64)

        print(magentaStr("\n==> Testing cryoDRGN - training vae:"))
        protTrain = self.newProtocol(CryoDrgnProtTrain, numEpochs=3, zDim=2)
        protTrain.inputParticles.set(preprocess.Particles)
        self.launchProtocol(protTrain)
        self.checkTrainOutput(protTrain)

        print(magentaStr("\n==> Testing cryoDRGN - ab initio (het):"))
        protAbinitio = self.newProtocol(CryoDrgnProtAbinitio, numEpochs=2, zDim=2)
        protAbinitio.inputParticles.set(preprocess.Particles)
        self.launchProtocol(protAbinitio)
        self.checkTrainOutput(protAbinitio)
