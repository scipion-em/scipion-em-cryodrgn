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

import pyworkflow.utils as pwutils
from pyworkflow.plugin import Domain
import pyworkflow.protocol.params as params
from pwem.constants import ALIGN_PROJ
from pwem.protocols import ProtProcessParticles

from cryodrgn import Plugin

convert = Domain.importFromPlugin('relion.convert', doRaise=True)
md = Domain.importFromPlugin('relion.convert.metadata', doRaise=True)


class CryoDrgnProtPreprocess(ProtProcessParticles):
    """
    Protocol to downsample a particle stack and prepare alignment/CTF parameters.

    Find more information at https://github.com/zhonge/cryodrgn
    """
    _label = 'preprocess'

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_parts': self._getExtraPath('input_particles.star'),
            'output_parts': self._getExtraPath('downsampled_parts.mrcs'),
            'output_poses': self._getExtraPath('pose.pkl'),
            'output_ctf': self._getExtraPath('ctf.pkl'),
        }

        self._updateFilenamesDict(myDict)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select a set of particles from a consensus '
                           '3D refinement.')
        form.addParam('doScale', params.BooleanParam, default=True,
                      label='Downsample particles?')
        form.addParam('scaleSize', params.IntParam, default=128,
                      condition='doScale',
                      validators=[params.Positive],
                      label='New box size (px)',
                      help='New box size in pixels, must be even.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep')

        if self.doScale:
            self._insertFunctionStep('runDownSampleStep')

        self._insertFunctionStep('runParsePosesStep')
        self._insertFunctionStep('runParseCtfStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Create a star file as expected by cryoDRGN."""
        imgSet = self.inputParticles.get()
        # Create links to binary files and write the relion .star file
        convert.writeSetOfParticles(
            imgSet, self._getFileName('input_parts'),
            outputDir=self._getExtraPath(), alignType=ALIGN_PROJ)

    def runDownSampleStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        params = ' '.join(self._getDownsampleArgs())
        program = Plugin.getProgram('downsample')
        self.runJob(program, params, env=Plugin.getEnviron(), cwd=None)

    def runParsePosesStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        params = ' '.join(self._getParsePosesArgs())
        program = Plugin.getProgram('parse_pose_star')
        self.runJob(program, params, env=Plugin.getEnviron(), cwd=None)

    def runParseCtfStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        params = ' '.join(self._getParseCtfArgs())
        program = Plugin.getProgram('parse_ctf_star')
        self.runJob(program, params, env=Plugin.getEnviron(), cwd=None)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        self._createFilenameTemplates()
        if not pwutils.exists(self._getFileName("output_ctf")):
            summary.append("Output not ready")
        else:
            summary.append("Created poses and ctf files for cryoDRGN.")

        return summary

    def _validate(self):
        errors = []

        particles = self._getInputParticles()
        if not particles.hasCTF():
            errors.append("The input has no CTF values!")

        if self.doScale and self.scaleSize.get() >= particles.getDim()[0]:
            errors.append("You cannot upscale particles!")

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getDownsampleArgs(self):
        args = ['%s ' % self._getFileName('input_parts'),
                '-o %s ' % self._getFileName('output_parts'),
                '--datadir %s' % self._getDataDir(),
                '--relion31']

        if self.doScale:
            args.append('-D %d ' % self.scaleSize.get())

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
                '-o %s ' % self._getFileName('output_ctf'),
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
            return self._getInputParticles().getDim()[0]

    def _getSamplingRate(self):
        inputSet = self._getInputParticles()
        oldSampling = inputSet.getSamplingRate()
        scaleFactor = self._getScaleFactor(inputSet)
        newSampling = oldSampling * scaleFactor

        return newSampling

    def _getScaleFactor(self, inputSet):
        xdim = inputSet.getDim()[0]
        scaleFactor = xdim / float(
            self.scaleSize.get() if self.doScale else xdim)

        return scaleFactor

    def _getDataDir(self):
        """ We assume all mrcs stacks are in the same folder. """
        part = self._getInputParticles().getFirstItem()
        _, fn = part.getLocation()

        return os.path.dirname(fn)

    def _getExtraCtfParams(self):
        """Remove once optics parsing is implemented in parse_ctf_star"""
        mdOptics = md.Table(fileName=self._getFileName('input_parts'),
                         tableName='optics')[0]
        cs = mdOptics.rlnSphericalAberration
        amp = mdOptics.rlnAmplitudeContrast
        kv = mdOptics.rlnVoltage

        mdParts = md.Table(fileName=self._getFileName('input_parts'),
                           tableName='particles')[0]
        ps = getattr(mdParts, 'rlnCtfPhaseShift', 0.)

        return cs, amp, kv, ps
