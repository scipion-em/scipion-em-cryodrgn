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
from pwem.protocols import ProtImportParticles
from pyworkflow.utils import magentaStr

from ..protocols import *


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
        cls.protPreprocess = cls.runPreprocess(cls.protImport.outputParticles)

    @classmethod
    def runPreprocess(cls, particles):
        print(magentaStr("\n==> Testing cryoDRGN - preprocess:"))
        protPreprocess = cls.newProtocol(CryoDrgnProtPreprocess,
                                         scaleSize=64)
        protPreprocess._createFilenameTemplates()
        protPreprocess.inputParticles.set(particles)
        cls.launchProtocol(protPreprocess)

        return protPreprocess

    def testTraining(self):
        print(magentaStr("\n==> Testing cryoDRGN - training:"))
        protTrain = self.newProtocol(CryoDrgnProtTrain, numEpochs=3)
        protTrain._createFilenameTemplates()
        protTrain.protPreprocess.set(self.protPreprocess)
        self.launchProtocol(protTrain)
