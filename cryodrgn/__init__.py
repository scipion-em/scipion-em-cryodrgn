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
import pwem
import pyworkflow.utils as pwutils
from pyworkflow import Config

from .constants import *


__version__ = '3.4'
_references = ['Zhong2020a', 'Zhong2020b']
_logo = "cryodrgn_logo.png"


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-cryodrgn"
    _supportedVersions = VERSIONS

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CRYODRGN_ENV_ACTIVATION, DEFAULT_ACTIVATION_CMD)

    @classmethod
    def getCryoDrgnEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(CRYODRGN_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch cryoDRGN. """
        environ = pwutils.Environ(os.environ)
        if 'PYTHONPATH' in environ:
            # this is required for python virtual env to work
            del environ['PYTHONPATH']
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
        for ver in VERSIONS:
            cls.addCryoDrgnPackage(env, ver,
                                   default=ver == CRYODRGN_DEFAULT_VER_NUM)

    @classmethod
    def addCryoDrgnPackage(cls, env, version, default=False):
        CRYODRGN_INSTALLED = 'cryodrgn_%s_installed' % version
        ENV_NAME = getCryoDrgnEnvName(version)
        # try to get CONDA activation command
        installCmds = [
            cls.getCondaActivationCmd(),
            'conda create -y -n %s python=3.7;' % ENV_NAME,  # Create the environment
            'conda activate %s;' % ENV_NAME,  # Activate the new environment and install downloaded code
            'conda install -y "pytorch<1.9.0" cudatoolkit -c pytorch &&',  # PyTorch 1.9+ is not working for now
            'conda install -y pandas seaborn scikit-learn &&',
            'conda install -y umap-learn jupyterlab ipywidgets cufflinks-py "nodejs>=15.12.0" -c conda-forge &&',
            'jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build &&',
            'jupyter labextension install jupyterlab-plotly --no-build &&',
            'jupyter labextension install plotlywidget --no-build &&',
            'jupyter lab build &&',
            'pip install -e . &&',
            'touch %s' % CRYODRGN_INSTALLED  # Flag installation finished
        ]

        cryodrgnCmds = [(" ".join(installCmds), CRYODRGN_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        tarPrefix = VERSION_PREFIX.get(version, '')
        env.addPackage('cryodrgn', version=version,
                       url='https://github.com/zhonge/cryodrgn/archive/refs/tags/%s%s.tar.gz' % (tarPrefix, version),
                       commands=cryodrgnCmds,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getCryoDrgnEnvActivation())

    @classmethod
    def getProgram(cls, program, gpus='0'):
        """ Create cryoDRGN command line. """
        fullProgram = '%s && CUDA_VISIBLE_DEVICES=%s cryodrgn %s' % (
            cls.getActivationCmd(), gpus, program)

        return fullProgram

    @classmethod
    def getActiveVersion(cls, *args):
        """ Return the env name that is currently active. """
        envVar = cls.getVar(CRYODRGN_ENV_ACTIVATION)
        return envVar.split()[-1].split("-")[-1]

    @classmethod
    def versionGE(cls, version):
        """ Return True if current version of cryodrgn is newer
         or equal than the input argument.
         Params:
            version: string version (semantic version, e.g 0.3.3b)
        """
        v1 = cls.getActiveVersion()
        if v1 not in VERSIONS:
            raise Exception("This version of cryoDRGN is not supported: ", v1)

        if VERSIONS.index(v1) < VERSIONS.index(version):
            return False
        return True
