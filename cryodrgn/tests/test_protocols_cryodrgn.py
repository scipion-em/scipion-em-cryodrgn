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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from pwem.protocols import ProtImportParticles

from ..protocols import CryoDrgnProtPreprocess, CryoDrgnProtTrain


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
        # cls.protPreprocess = cls.runPreprocess(cls.protImport.outputParticles)

    def runPreprocess(self, particles, **kwargs):
        print(magentaStr("\n==> Testing cryoDRGN - preprocess:"))
        protPreprocess = self.newProtocol(CryoDrgnProtPreprocess, **kwargs)
        # protPreprocess._createFilenameTemplates()
        protPreprocess.inputParticles.set(particles)
        return self.launchProtocol(protPreprocess)

    def checkPreprocessOutput(self, preprocessProt):
        """ Do some check on the output of the preprocess. """
        output = getattr(preprocessProt, 'outputCryoDrgnParticles', None)
        self.assertIsNotNone(output)

        filename = output.filename.get()
        poses = output.poses.get()
        ctfs = output.ctfs.get()

        for f in [filename, poses, ctfs]:
            fn = os.path.join(self.proj.path, f)
            self.assertTrue(os.path.exists(fn))

        self.assertTrue(filename.endswith('.txt'))

    def testPreprocess(self):
        parts = self.protImport.outputParticles

        preprocess1 = self.runPreprocess(parts, scaleSize=64)
        self.checkPreprocessOutput(preprocess1)

        preprocess2 = self.runPreprocess(parts, scaleSize=50, chunk=200)
        self.checkPreprocessOutput(preprocess2)

    def testTraining(self):
        parts = self.protImport.outputParticles

        preprocess = self.runPreprocess(parts, scaleSize=64)

        print(magentaStr("\n==> Testing cryoDRGN - training:"))

        protTrain = self.newProtocol(CryoDrgnProtTrain, numEpochs=3)
        protTrain.inputParticles.set(preprocess.outputCryoDrgnParticles)
        self.launchProtocol(protTrain)
