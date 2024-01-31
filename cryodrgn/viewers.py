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

from pyworkflow.protocol.params import LabelParam, EnumParam, BooleanParam
from pyworkflow.protocol.executor import StepExecutor
from pyworkflow.viewer import DESKTOP_TKINTER
from pwem.viewers import ObjectView, ChimeraView, EmProtocolViewer

from cryodrgn import Plugin
from cryodrgn.protocols import CryoDrgnProtAnalyze
from cryodrgn.constants import VOLUME_SLICES, VOLUME_CHIMERA


class CryoDrgnViewer(EmProtocolViewer):
    """ Visualization of cryoDRGN results. """

    _environments = [DESKTOP_TKINTER]
    _targets = [CryoDrgnProtAnalyze]
    _label = 'analyze results'

    def __init__(self, **kwargs):
        EmProtocolViewer.__init__(self, **kwargs)

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        if self.protocol.hasMultLatentVars():
            path = glob(self.protocol.getOutputDir("analyze.*/kmeans*"))[0]
            self._updateFilenamesDict({
                'output_umap': os.path.join(path, 'umap.png'),
                'output_umaphex': os.path.join(path, 'umap_hex.png'),
                'output_pca': os.path.join(path, 'z_pca.png'),
                'output_pcahex': os.path.join(path, 'z_pca_hex.png'),
                'output_notebook': os.path.join(path, '../cryoDRGN_viz.ipynb')
            })
        else:
            path = glob(self.protocol.getOutputDir("analyze.*"))[0]
            self._updateFilenamesDict({
                'output_hist': os.path.join(path, 'z_hist.png'),
                'output_dist': os.path.join(path, 'z.png'),
                'output_notebook': os.path.join(path, 'cryoDRGN_viz.ipynb')
            })

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayVol', EnumParam,
                      choices=['slices', 'chimera'],
                      default=VOLUME_SLICES,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Display K-means volumes with',
                      help="Cryodrgn analyze uses the k-means clustering algorithm to "
                           "partition the latent space into k regions (by default k=20). "
                           "A density map is generated from the center of each of these "
                           "clusters. The kmeans20 volumes provide an initial set of "
                           "representative density maps to visually inspect (and not "
                           "necessarily to assign classes).")
        if self.protocol.hasMultLatentVars():
            form.addParam('useHexBin', BooleanParam, default=False,
                          label="Use hexagonal bins for the plots?")
            form.addParam('doShowPCA', LabelParam,
                          label='Show PCA projection of latent space encodings')
            if not self.protocol.skipUmap:
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

        if self.protocol.hasMultLatentVars():
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
                return self._showVolumeSlices()
        except Exception as e:
            self.showError(str(e))

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        prot = self.protocol
        extra = prot._getExtraPath()
        cmdFile = prot._getExtraPath('chimera_volumes.cxc')

        with open(cmdFile, 'w+') as f:
            for vol in prot.getFiles():
                localVol = os.path.relpath(vol, extra)
                if os.path.exists(vol):
                    f.write(f"open {localVol}\n")
            f.write('tile\n')

        view = ChimeraView(cmdFile)

        return [view]

    def _showVolumeSlices(self):
        """ Open a sqlite with all volumes selected for visualization. """
        path = self.protocol._getExtraPath('volumes.sqlite')
        return [ObjectView(self._project, self.protocol.strId(), path)]

    def _showPlot(self, fn):
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt
        fn = self._getFileName(fn)
        if os.path.exists(fn):
            img = mpimg.imread(fn)
            imgplot = plt.imshow(img)
            plt.axis('off')
            plt.show()
            return [imgplot]
        else:
            self.showError(f"File {fn} not found! Have you run analysis?")

    def _showHistogram(self, paramName=None):
        self._showPlot('output_hist')

    def _showDistribution(self, paramName=None):
        self._showPlot('output_dist')

    def _showPCA(self, paramName=None):
        fn = 'output_pcahex' if self.useHexBin else 'output_pca'
        self._showPlot(fn)

    def _showUMAP(self, paramName=None):
        fn = 'output_umaphex' if self.useHexBin else 'output_umap'
        self._showPlot(fn)

    def _showNotebook(self, paramName=None):
        """ Open jupyter notebook with results in a browser. """

        def _extraWork():
            program = Plugin.getProgram('').split()[:-1]  # remove cryodrgn command
            fn = self._getFileName('output_notebook')

            if self.serverMode:
                args = '--no-browser --port 8888 '
            else:
                args = f"{os.path.basename(fn)}"

            program.append(f"jupyter notebook {args}")

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
