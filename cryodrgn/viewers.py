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

from pyworkflow.protocol.params import (LabelParam, EnumParam,
                                        IntParam, BooleanParam)
from pyworkflow.protocol.executor import StepExecutor
from pyworkflow.viewer import DESKTOP_TKINTER
from pwem.viewers import ObjectView, ChimeraView, EmProtocolViewer

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

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        prot = self.protocol
        self._epoch = self.getEpoch()
        epochDir = os.path.join(prot.getOutputDir(), 'analyze.%d' % self._epoch)

        def _out(f):
            return os.path.join(epochDir, f)

        self._updateFilenamesDict({
            'output_notebook': _out('cryoDRGN_viz.ipynb'),
            'output_hist': _out('z_hist.png'),
            'output_dist': _out('z.png'),
            'output_umap': _out('umap.png'),
            'output_umaphex': _out('umap_hexbin.png'),
            'output_pca': _out('z_pca.png'),
            'output_pcahex': _out('z_pca_hexbin.png'),
            'output_vol': _out('vol_%(id)03d.mrc'),
            'output_volN': _out('kmeans20/vol_%(id)03d.mrc'),
        })

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewEpoch', EnumParam,
                      choices=['last', 'selection'], default=EPOCH_LAST,
                      display=EnumParam.DISPLAY_LIST,
                      label="Epoch to analyze")
        form.addParam('epochNum', IntParam,
                      condition='viewEpoch==%d' % EPOCH_SELECTION,
                      label="Epoch number")

        form.addParam('displayVol', EnumParam,
                      choices=['slices', 'chimera'],
                      default=VOLUME_SLICES,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Display volume with',
                      help='*slices*: display volumes as 2D slices along z axis.\n'
                           '*chimera*: display volumes as surface with Chimera.')
        if self.hasMultLatentVars():
            form.addParam('useHexBin', BooleanParam, default=False,
                          label="Use hexagonal bins for the plots?")
            form.addParam('doShowPCA', LabelParam,
                          label='Show PCA projection of latent space encodings')
            form.addParam('doShowUMAP', LabelParam,
                          label="Show UMAP visualization of latent space encodings")
        else:
            form.addParam('doShowDistr', LabelParam,
                          label='Show latent coordinates distribution')
            form.addParam('doShowHistogram', LabelParam,
                          label="Show latent coordinates histogram")

        group = form.addGroup('Jupyter')
        group.addParam('serverMode', BooleanParam, default=False,
                       label='Launch Jupyter in server mode?',
                       help="If yes, the Jupyter notebook will be initialized "
                            "with option *--no-browser*.\n"
                            "Then you might be able to connect to the local server. "
                            "One can also access the server remotely by setting "
                            "a SSH tunnel:\n"
                            "$ ssh -N -f -L localhost:8888:localhost:8888 remote_username@remote_host_name "
                            "# replace remote_username and remote_host_name with your login information")
        group.addParam('doShowNotebook', LabelParam,
                       label="Show Jupyter notebook")

    def _getVisualizeDict(self):
        self._createFilenameTemplates()

        visDict = {
            'displayVol': self._showVolumes,
            'doShowNotebook': self._showNotebook
        }

        if self.hasMultLatentVars():
            visDict.update({
                'doShowPCA': self._showPCA,
                'doShowUMAP': self._showUMAP
            })
        else:
            visDict.update({
                'doShowDistr': self._showDistribution,
                'doShowHistogram': self._showHistogram
            })

        return visDict

    def _showVolumes(self, paramName=None):
        try:
            if self.displayVol == VOLUME_CHIMERA:
                return self._showVolumesChimera()
            elif self.displayVol == VOLUME_SLICES:
                return self._createVolumesSqlite()
        except Exception as e:
            self.showError(str(e))

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        prot = self.protocol
        volumes = self._getVolumeNames()
        extra = prot._getExtraPath()
        cmdFile = prot._getExtraPath('chimera_volumes.cxc')

        with open(cmdFile, 'w+') as f:
            for vol in volumes:
                localVol = os.path.relpath(vol, extra)
                if os.path.exists(vol):
                    f.write("open %s\n" % localVol)
            f.write('tile\n')

        view = ChimeraView(cmdFile)

        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """
        prot = self.protocol  # shortcut
        path = prot._getExtraPath('viewer_volumes.sqlite')
        samplingRate = prot.inputParticles.get().getSamplingRate()
        files = self._getVolumeNames()
        self.createVolumesSqlite(files, path, samplingRate)

        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showPlot(self, fn, epoch):
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt
        fn = self._getFileName(fn, epoch=epoch)
        if os.path.exists(fn):
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

    def _showPCA(self, paramName=None):
        fn = 'output_pcahex' if self.useHexBin else 'output_pca'
        self._showPlot(fn, epoch=self._epoch)

    def _showUMAP(self, paramName=None):
        fn = 'output_umaphex' if self.useHexBin else 'output_umap'
        self._showPlot(fn, epoch=self._epoch)

    def _showNotebook(self, paramName=None):
        """ Open jupyter notebook with results in a browser. """

        def _extraWork():
            program = Plugin.getProgram('').split()[:-1]  # remove cryodrgn command
            fn = self._getFileName('output_notebook', epoch=self._epoch)

            if self.serverMode:
                args = '--no-browser --port 8888 '
            else:
                args = '%s ' % os.path.basename(fn)

            program.append('jupyter notebook %s' % args)

            if os.path.exists(fn):
                fnDir = os.path.dirname(fn)
                hostConfig = self.protocol.getHostConfig()
                executor = StepExecutor(hostConfig)
                self.protocol.setStepsExecutor(executor)
                self.protocol.runJob(" ".join(program), '',
                                     env=Plugin.getEnviron(),
                                     cwd=fnDir)
            else:
                self.showError('Jupyter notebook not found! Have you run analysis?')

        from threading import Thread
        thread = Thread(target=_extraWork)
        thread.start()

    def _getVolumeNames(self):
        vols = []
        if self.hasMultLatentVars():
            fn = 'output_volN'
            num = 20
        else:
            fn = 'output_vol'
            num = 10

        for volId in range(num):
            volFn = self._getFileName(fn, epoch=self._epoch, id=volId)
            if os.path.exists(volFn):
                vols.append(volFn)
            else:
                raise FileNotFoundError("Volume %s does not exists. \n"
                                        "Please select a valid epoch "
                                        "number." % volFn)
        return vols

    def getEpoch(self):
        """ Get the epoch number for visualisation. """
        if self.viewEpoch == EPOCH_LAST:
            return self.protocol.getLastEpoch()
        else:
            return self.epochNum.get()

    def hasMultLatentVars(self):
        return self.protocol.zDim > 1
