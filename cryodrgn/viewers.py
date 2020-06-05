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
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from pyworkflow.protocol.params import LabelParam, EnumParam, IntParam
from pyworkflow.protocol.executor import StepExecutor
import pyworkflow.utils as pwutils
from pyworkflow.viewer import DESKTOP_TKINTER
from pwem.viewers import (ChimeraClientView, ObjectView,
                          ChimeraView, EmProtocolViewer)

from cryodrgn import Plugin
from .protocols import CryoDrgnProtTrain
from .constants import *


class CryoDrgnViewer(EmProtocolViewer):
    """ Visualization of cryoDRGN results. """
           
    _environments = [DESKTOP_TKINTER]
    _targets = [CryoDrgnProtTrain]
    _label = 'analyze results'

    def __init__(self, **kwargs):
        EmProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewEpoch', EnumParam,
                      choices=['last', 'selection'], default=EPOCH_LAST,
                      display=EnumParam.DISPLAY_LIST,
                      label="Epoch to analyze")
        form.addParam('epochNum', IntParam,
                      condition='viewEpoch==%d' % EPOCH_SELECTION,
                      label="Epoch number")
        form.addParam('runAnalyze', LabelParam, important=True,
                      label="RUN ANALYSIS!",
                      help="This script runs a series of standard analyses:\n\n"
                           "- PCA of the latent space\n"
                           "- UMAP embedding of the latent space\n"
                           "- Generation of volumes from the latent space\n"
                           "- Generation of trajectories along the first "
                           "and second principal components\n"
                           "- Generation of a template jupyter notebook "
                           "that may be used for further interactive "
                           "analyses and visualization")

        group = form.addGroup('Visualize')
        group.addParam('displayVol', EnumParam,
                       choices=['slices', 'chimera'],
                       default=VOLUME_SLICES,
                       display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('doShowDistr', LabelParam,
                       label='Show latent coordinates distribution')
        group.addParam('doShowHistogram', LabelParam,
                       label="Show latent coordinates histogram")
        group.addParam('doShowNotebook', LabelParam,
                       label="Show Jupyter notebook")

    def _getVisualizeDict(self):
        self.protocol._initialize()  # Load filename templates
        self._loadEpochs()
        return {'runAnalyze': self._runAnalysis,
                'displayVol': self._showVolumes,
                'doShowDistr': self._showDistribution,
                'doShowHistogram': self._showHistogram,
                'doShowNotebook': self._showNotebook
                }

    def _runAnalysis(self, paramName=None):
        program = Plugin.getProgram('analyze')
        args = ' %s %d --Apix %0.3f' % (
            self.protocol._getFileName("output_dir"),
            self._epoch,
            self.protocol._getSampling())

        hostConfig = self.protocol.getHostConfig()
        # Create the steps executor
        executor = StepExecutor(hostConfig)
        self.protocol.setStepsExecutor(executor)
        # Finally run the protocol
        self.protocol.runJob(program, args, env=Plugin.getEnviron(),
                             cwd=None)

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()

        if len(volumes) > 1:
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cmd')
            with open(cmdFile, 'w+') as f:
                for vol in volumes:
                    localVol = os.path.basename(vol)
                    if pwutils.exists(vol):
                        f.write("open %s\n" % localVol)
                f.write('tile\n')
            view = ChimeraView(cmdFile)
        else:
            view = ChimeraClientView(volumes[0])

        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('viewer_volumes.sqlite')
        samplingRate = self.protocol._getSampling()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            print("Adding vol: %s" % volFn)
            files.append(volFn)

        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showPlot(self, fn, epoch):
        fn = self.protocol._getFileName(fn, epoch=epoch)
        if pwutils.exists(fn):
            img = mpimg.imread(fn)
            imgplot = plt.imshow(img)
            plt.axis('off')
            plt.show()
            return [imgplot]
        else:
            self.showError('File %s not found! Have you run analysis?' % fn)

    def _showHistogram(self, paramName=None):
        self._showPlot('output_hist', epoch=self._epoch)

    def _showDistribution(self, paramName=None):
        self._showPlot('output_dist', epoch=self._epoch)

    def _showNotebook(self, paramName=None):
        program = Plugin.getProgram('').split()[:-1]  # remove cryodrgn command
        fn = self.protocol._getFileName('output_notebook',
                                        epoch=self._epoch)
        program.append('jupyter notebook %s' % os.path.basename(fn))

        if pwutils.exists(fn):
            fnDir = os.path.dirname(fn)
            hostConfig = self.protocol.getHostConfig()
            executor = StepExecutor(hostConfig)
            self.protocol.setStepsExecutor(executor)
            self.protocol.runJob(" ".join(program), '',
                                 env=Plugin.getEnviron(),
                                 cwd=fnDir)
        else:
            self.showError('Jupyter notebook not found! Have you run analysis?')


    def _getVolumeNames(self):
        vols = []
        for volId in range(10):  # FIXME: is it always 10 volumes?
            volFn = self.protocol._getFileName('output_vol', epoch=self._epoch,
                                               id=volId)
            if pwutils.exists(volFn):
                vols.append(volFn)
            else:
                raise Exception("Volume %s does not exists. \n"
                                "Please select a valid epoch "
                                "number OR *Run analysis!* first." % volFn)
        return vols

    def _loadEpochs(self):
        if self.viewEpoch.get() == EPOCH_LAST:
            self._epoch = self.protocol._lastIter()
            self.protocol._getEpochNumber(-1)
        else:
            self._epoch = self.epochNum.get()
