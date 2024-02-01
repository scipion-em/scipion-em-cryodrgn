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

from pyworkflow.protocol.params import (LabelParam, EnumParam, BooleanParam,
                                        IntParam)
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
            out = lambda p: os.path.join(path, p)
            self._updateFilenamesDict({
                'output_umap': out('umap.png'),
                'output_umaphex': out('umap_hex.png'),
                'output_pca': out('z_pca.png'),
                'output_pcahex': out('z_pca_hex.png'),
                'output_pca_volN': out('../pc%(pc)d/vol_%(id)03d.mrc'),
                'output_umap_pcN': out('../pc%(pc)d/umap.png'),
                'output_umap_pcN_traversal': out('../pc%(pc)d/umap_traversal_connected.png'),
                'output_graph_vol': out('../graph_traversal/vol_000.mrc'),
                'output_notebook': out('../cryoDRGN_viz.ipynb')
            })
        else:
            path = glob(self.protocol.getOutputDir("analyze.*"))[0]
            out = lambda p: os.path.join(path, p)
            self._updateFilenamesDict({
                'output_hist': out('z_hist.png'),
                'output_dist': out('z.png'),
                'output_notebook': out('cryoDRGN_viz.ipynb')
            })

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayType', EnumParam,
                      choices=['slices', 'chimera'],
                      default=VOLUME_SLICES,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Display volumes with')
        form.addParam('displayVolKmeans', LabelParam,
                      label='Show K-means volumes',
                      help="Cryodrgn analyze uses the k-means clustering algorithm to "
                           "partition the latent space into k regions (by default k=20). "
                           "A density map is generated from the center of each of these "
                           "clusters. The kmeans20 volumes provide an initial set of "
                           "representative density maps to visually inspect (and not "
                           "necessarily to assign classes).")
        if self.protocol.hasMultLatentVars():
            if self.protocol.doGraphTraversal:
                form.addParam('displayVolGraph', LabelParam,
                              label='Show graph traversal movie in ChimeraX')

            form.addParam('useHexBin', BooleanParam, default=False,
                          label="Use hexagonal bins for the plots?")

            group = form.addGroup('PCA & UMAP')
            group.addParam('doShowPCA', LabelParam,
                           label='Show PCA projection of latent space encodings')
            group.addParam('pcNum', IntParam, default=1,
                           label="Which PC to plot",
                           help="The plots within the pcX subdirectories show the UMAP "
                                "embedding colored by each particle's projected value "
                                "along PCX. This helps give a sense of the layout of the "
                                "latent space and how the UMAP embedding (e.g. a nonlinear "
                                "8D → 2D embedding) is related to the PCA projection (a "
                                "linear projection from 8D → 2D).")
            group.addParam('displayVolPCA', LabelParam,
                           label='Show PCX volumes with',
                           help="By default, the 10 volumes are generated at equally spaced "
                                "points between the first and 99th percentile of the data "
                                "distribution projected onto each principal component.")

            if not self.protocol.skipUmap:
                group.addParam('doShowUMAP', LabelParam,
                               label="Show UMAP visualization of latent space encodings")

            group.addParam('doShowPcaToUmap', LabelParam,
                           label="Show UMAP embeddings colored by PCX value")
            group.addParam('doShowPcaToUmapTrav', LabelParam,
                           label="Show UMAP connected traversal along PCX")

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
            'displayVolKmeans': lambda paramName: self.showVolumes(key='kmeans'),
            'doShowNotebook': self.showNotebook
        }

        if self.protocol.hasMultLatentVars():
            visDict.update({
                'displayVolPCA': lambda paramName: self.showVolumes(key='pca'),
                'doShowPCA': self.showPCA,
                'doShowUMAP': self.showUMAP,
                'doShowPcaToUmap': self.showPCAToUMAP,
                'doShowPcaToUmapTrav': self.showPCAToUMAPTrav
            })
            if self.protocol.doGraphTraversal:
                visDict['displayVolGraph'] = lambda paramName: self.showVolumes(key='graph')
        else:
            visDict.update({
                'doShowDistr': self.showDistribution,
                'doShowHistogram': self.showHistogram
            })

        return visDict

    def showVolumes(self, key):
        try:
            if self.displayType == VOLUME_CHIMERA:
                return self.showVolumesChimera(key)
            elif self.displayType == VOLUME_SLICES:
                return self.showVolumeSlices(key)
        except Exception as e:
            self.showError(str(e))

    def showVolumesChimera(self, key='kmeans'):
        """ Create a chimera script to visualize selected volumes. """
        prot = self.protocol
        if key == 'kmeans':
            vols = list(prot.Volumes.getFiles())
        elif key == 'pca':
            vols = self._getVolumesNamesPCA()
        elif key == 'graph':
            vols = [self._getFileName('output_graph_vol')]
        else:
            raise KeyError("Unknown volume type")

        cmdFile = prot._getExtraPath('chimera_volumes.cxc')
        localVol = os.path.relpath(vols[0], prot._getExtraPath())

        if os.path.exists(vols[0]):
            with open(cmdFile, 'w+') as f:
                f.write(f"open {os.path.dirname(localVol)}/vol_*.mrc vseries true\n")
                if key == 'graph':
                    f.write("vol all color cornflowerblue\n"
                            "mov record\n"
                            "mseries all\n"
                            "mov encode graph_traversal.mp4 framerate 15\n")
        else:
            raise FileNotFoundError(f"File {vols[0]} not found!")

        view = ChimeraView(cmdFile)

        return [view]

    def showVolumeSlices(self, key='kmeans'):
        """ Open a sqlite with all volumes selected for visualization. """
        if key == 'pca':
            path = self.protocol._getExtraPath(f"volumes_pc{self.pcNum}.sqlite")
            files = self._getVolumesNamesPCA()
            samplingRate = self.protocol._getOutputSampling()
            self.createVolumesSqlite(files, path, samplingRate)
        elif key == 'kmeans':
            path = self.protocol._getExtraPath('volumes.sqlite')
        else:  # only chimerax allowed for graph volumes
            return self.showVolumesChimera(key='graph')

        return [ObjectView(self._project, self.protocol.strId(), path)]

    def showHistogram(self, paramName=None):
        self._showPlot('output_hist')

    def showDistribution(self, paramName=None):
        self._showPlot('output_dist')

    def showPCA(self, paramName=None):
        fn = 'output_pcahex' if self.useHexBin else 'output_pca'
        self._showPlot(fn)

    def showUMAP(self, paramName=None):
        fn = 'output_umaphex' if self.useHexBin else 'output_umap'
        self._showPlot(fn)

    def showPCAToUMAP(self, paramName=None):
        self._showPlot('output_umap_pcN', pc=self.pcNum)

    def showPCAToUMAPTrav(self, paramName=None):
        self._showPlot('output_umap_pcN_traversal', pc=self.pcNum)

    def showNotebook(self, paramName=None):
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
                self.showError(f"Jupyter notebook {fn} not found!")

        from threading import Thread
        thread = Thread(target=_extraWork)
        thread.start()

    # --------------------------- UTILS functions -----------------------------
    def _getVolumesNamesPCA(self):
        """ Get filenames for 10 output volumes along PCX. """
        names = []
        vols = [self._getFileName('output_pca_volN',
                                  pc=self.pcNum,
                                  id=i) for i in range(10)]
        for fn in vols:
            if os.path.exists(fn):
                names.append(fn)
            else:
                raise FileNotFoundError(f"File {fn} not found!")

        return names

    def _showPlot(self, fn, **kwargs):
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt
        fn = self._getFileName(fn, **kwargs)
        if os.path.exists(fn):
            img = mpimg.imread(fn)
            imgplot = plt.imshow(img)
            plt.axis('off')
            plt.show()
            return [imgplot]
        else:
            self.showError(f"File {fn} not found!")
