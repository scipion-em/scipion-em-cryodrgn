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
from emtable import Table

import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtProcessParticles

from cryodrgn import Plugin
from cryodrgn.objects import CryoDrgnParticles

convert = Domain.importFromPlugin('relion.convert', doRaise=True)


class CryoDrgnProtPreprocess(ProtProcessParticles):
    """ Protocol to downsample a particle stack and prepare alignment/CTF parameters.

    Find more information at https://github.com/zhonge/cryodrgn
    """
    _label = 'preprocess'
    _devStatus = BETA

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        def out(*p):
            return os.path.join(self._getPath('output_particles'), *p)

        myDict = {
            'input_parts': self._getExtraPath('input_particles.star'),
            'output_folder': out(),
            'output_parts': out('particles.%d.mrcs' % self._getBoxSize()),
            'output_poses': out('poses.pkl'),
            'output_ctfs': out('ctfs.pkl'),
        }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select a set of particles from a consensus C1 '
                           '3D refinement.')
        form.addParam('doScale', params.BooleanParam, default=True,
                      label='Downsample particles?')
        form.addParam('scaleSize', params.IntParam, default=128,
                      condition='doScale',
                      validators=[params.Positive],
                      label='New box size (px)',
                      help='New box size in pixels, must be even.')

        form.addParam('chunk', params.IntParam, default=0,
                      label='Split in chunks',
                      help='If this value is greater than 0, the output stack'
                           'will be saved into parts of this size. This will '
                           'avoid out-of-memory errors when saving out a large '
                           'particle stack. (param **--chunk**)\n'
                           'For example, use --chunk 50000 to chunk the output '
                           'into separate .mrcs with 50k images each. ')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')

        if self.doScale or self.chunk > 0:
            self._insertFunctionStep('runDownSampleStep')

        self._insertFunctionStep('runParsePosesStep')
        self._insertFunctionStep('runParseCtfStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a star file as expected by cryoDRGN."""
        outputFolder = self._getFileName('output_folder')
        pwutils.cleanPath(outputFolder)
        pwutils.makePath(outputFolder)

        imgSet = self.inputParticles.get()
        # Create links to binary files and write the relion .star file
        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_parts'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)

    def runDownSampleStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        self._runProgram('downsample', self._getDownsampleArgs())

    def runParsePosesStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        self._runProgram('parse_pose_star', self._getParsePosesArgs())

    def runParseCtfStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        self._runProgram('parse_ctf_star', self._getParseCtfArgs())

    def createOutputStep(self):
        outputParts = self._getFileName('output_parts')
        if self.chunk > 0:
            outputParts = outputParts.replace('.mrcs', '.txt')

        if not os.path.exists(outputParts):
            outputParts = self._getFileName('input_parts')

        output = CryoDrgnParticles(filename=outputParts,
                                   poses=self._getFileName('output_poses'),
                                   ctfs=self._getFileName('output_ctfs'),
                                   dim=self._getBoxSize(),
                                   samplingRate=self._getSamplingRate())

        self._defineOutputs(outputCryoDrgnParticles=output)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        self._createFilenameTemplates()
        if not os.path.exists(self._getFileName("output_ctfs")):
            summary.append("Output not ready")
        else:
            summary.append("Created poses and ctf files for cryoDRGN.")

        return summary

    def _validate(self):
        errors = []

        particles = self._getInputParticles()
        if not particles.hasCTF():
            errors.append("The input has no CTF values!")

        if self.doScale and self.scaleSize > particles.getXDim():
            errors.append("You cannot upscale particles!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getDownsampleArgs(self):
        args = ['%s ' % self._getFileName('input_parts'),
                '-o %s ' % self._getFileName('output_parts'),
                '--datadir %s' % self._getDataDir(),
                '--relion31',
                '-D %d' % self._getBoxSize()]

        if self.chunk > 0:
            args.append('--chunk %d ' % self.chunk)

        return args

    def _getParsePosesArgs(self):
        args = ['%s ' % self._getFileName('input_parts'),
                '-o %s ' % self._getFileName('output_poses'),
                '-D %d' % self._getBoxSize(),
                '--Apix %0.3f' % self._getSamplingRate(),
                '--relion31']

        return args

    def _getParseCtfArgs(self):
        args = ['%s ' % self._getFileName('input_parts'),
                '-o %s ' % self._getFileName('output_ctfs'),
                '-D %d' % self._getBoxSize(),
                '--Apix %0.3f' % self._getSamplingRate(),
                '--relion31']

        cs, amp, kv, ps = self._getExtraCtfParams()
        args.extend(['--kv %f ' % kv,
                     '--cs %f' % cs,
                     '-w %f' % amp])

        if ps == 0:  # no phase shift found
            args.append('--ps 0')

        return args

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _getBoxSize(self):
        if self.doScale:
            return self.scaleSize.get()
        else:
            return self._getInputParticles().getXDim()

    def _getSamplingRate(self):
        inputSet = self._getInputParticles()
        oldSampling = inputSet.getSamplingRate()
        scaleFactor = self._getScaleFactor(inputSet)
        newSampling = oldSampling * scaleFactor

        return newSampling

    def _getScaleFactor(self, inputSet):
        return inputSet.getXDim() / float(self._getBoxSize())

    def _getDataDir(self):
        """ We assume all mrcs stacks are in the same folder. """
        part = self._getInputParticles().getFirstItem()
        _, fn = part.getLocation()

        return os.path.dirname(fn)

    def _getExtraCtfParams(self):
        """Remove once optics parsing is implemented in parse_ctf_star"""
        mdOptics = Table(fileName=self._getFileName('input_parts'),
                         tableName='optics')[0]
        cs = mdOptics.rlnSphericalAberration
        amp = mdOptics.rlnAmplitudeContrast
        kv = mdOptics.rlnVoltage

        mdParts = Table(fileName=self._getFileName('input_parts'),
                        tableName='particles')[0]
        ps = getattr(mdParts, 'rlnCtfPhaseShift', 0.)

        return cs, amp, kv, ps

    def _runProgram(self, program, args):
        self.runJob(Plugin.getProgram(program), ' '.join(args))
