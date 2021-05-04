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
from glob import glob
import re
from emtable import Table

import pyworkflow.utils as pwutils
from pyworkflow.constants import VERSION_3_0
import pyworkflow.protocol.params as params
from pwem.protocols import ProtProcessParticles

from cryodrgn import Plugin


class CryoDrgnProtTrain(ProtProcessParticles):
    """
    Protocol to train cryoDRGN neural network.

    Find more information at https://github.com/zhonge/cryodrgn
    """
    _lastUpdateVersion = VERSION_3_0
    _label = 'training'

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    def _initialize(self):
        self._createFilenameTemplates()
        self._createEpochTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'output_dir': self._getExtraPath('output'),
            'output_notebook': self._getExtraPath('output/analyze.%(epoch)d/cryoDRGN_viz.ipynb'),
            'output_hist': self._getExtraPath('output/analyze.%(epoch)d/z_hist.png'),
            'output_dist': self._getExtraPath('output/analyze.%(epoch)d/z.png'),
            'output_umap': self._getExtraPath('output/analyze.%(epoch)d/umap.png'),
            'output_umaphex': self._getExtraPath('output/analyze.%(epoch)d/umap_hexbin.png'),
            'output_pca': self._getExtraPath('output/analyze.%(epoch)d/z_pca.png'),
            'output_pcahex': self._getExtraPath('output/analyze.%(epoch)d/z_pca_hexbin.png'),
            'output_vol': self._getExtraPath('output/analyze.%(epoch)d/vol_%(id)03d.mrc'),
            'output_volN': self._getExtraPath('output/analyze.%(epoch)d/kmeans20/vol_%(id)03d.mrc'),
            'output_z': self._getExtraPath('output/z.%(z)d.pkl')
        }

        self._updateFilenamesDict(myDict)

    def _createEpochTemplates(self):
        """ Setup the regex on how to find epochs. """
        self._epochTemplate = self._getFileName('output_z', z=0).replace('z.0', 'z.*')
        self._epochRegex = re.compile('z\.(\d+)\.pkl')

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
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

        form.addParam('zDim', params.IntParam, default=1,
                      validators=[params.Positive],
                      label='Dimension of latent variable',
                      help='It is recommended to first train on lower '
                           'resolution images (e.g. D=128) with '
                           '--zdim 1 and with --zdim 10 using the '
                           'default architecture (fast).')
        form.addParam('numEpochs', params.IntParam, default=20,
                      label='Number of epochs',
                      help='The number of epochs refers to the number '
                           'of full passes through the dataset for '
                           'training, and should be modified depending '
                           'on the number of particles in the dataset. '
                           'For a 100k particle dataset, the above '
                           'settings required ~6 min per epoch for D=128 '
                           'images + default architecture, ~12 min/epoch '
                           'for D=128 images + large architecture, and ~47 '
                           'min per epoch for D=256 images + large architecture.')
        form.addParam('doInvert', params.BooleanParam, default=True,
                      label="Are particles white?")

        form.addSection(label='Advanced')
        group = form.addGroup('Encoder')
        group.addParam('qLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('qDim', params.IntParam, default=256,
                       label='Number of nodes in hidden layers')

        group = form.addGroup('Decoder')
        group.addParam('pLayers', params.IntParam, default=3,
                       label='Number of hidden layers')
        group.addParam('pDim', params.IntParam, default=256,
                       label='Number of nodes in hidden layers')

        form.addParam('extraParams', params.StringParam, default="",
                      label="Extra params",
                      help="Here you can provide all extra command-line "
                           "parameters. See *cryodrgn train_vae -h* for help.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('runTrainingStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """ Copy files as expected by cryoDRGN."""
        pwutils.cleanPath(self._getFileName('output_dir'))
        pwutils.makePath(self._getFileName('output_dir'))

    def runTrainingStep(self):
        """ Call cryoDRGN with the appropriate parameters. """
        params = ' '.join(self._getTrainingArgs())
        gpus = ','.join(str(i) for i in self.getGpuList())
        program = Plugin.getProgram('train_vae', gpus=gpus)
        self.runJob(program, params)

    def createOutputStep(self):
        pass

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = ["Training VAE for %d epochs." % self.numEpochs.get()]

        return summary

    def _validate(self):
        errors = []

        return errors

    # --------------------------- UTILS functions -----------------------------
    def _getTrainingArgs(self):
        parts = self.inputParticles.get()

        args = [
            parts.filename.get(),
            '--poses %s' % parts.poses,
            '--ctf %s' % parts.ctfs,
            '-o %s ' % self._getFileName('output_dir'),
            '--zdim %d' % self.zDim,
            '-n %d' % self.numEpochs,
            '--relion31 '
        ]

        if len(self.getGpuList()) > 1:
            args.append('--multigpu')
        if not self.doInvert:
            args.append('--uninvert-data')

        args.extend([
            '--enc-layers %d' % self.qLayers,
            '--enc-dim %d' % self.qDim,
            '--dec-layers %d' % self.pLayers,
            '--dec-dim %d' % self.pDim
        ])

        return args

    def _getDataDir(self):
        """ We assume all mrcs stacks are in the same folder. """
        mdOptics = Table(fileName=self._getFileName('input_parts'),
                         tableName='particles')
        row = mdOptics[0]
        location = str(row.rlnImageName)

        return os.path.dirname(location.split('@')[1])

    def _getEpochNumber(self, index):
        """ Return the list of epoch files, given the epochTemplate. """
        result = None
        files = sorted(glob(self._epochTemplate),
                       key=lambda x: int(re.findall("z\.(\d+)\.", x)[0]))
        if files:
            f = files[index]
            s = self._epochRegex.search(f)
            if s:
                result = int(s.group(1))  # group 1 is a digit epoch number
        return result

    def _lastIter(self):
        return self._getEpochNumber(-1)
