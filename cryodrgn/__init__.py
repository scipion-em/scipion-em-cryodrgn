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


__version__ = '3.0.3'
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
        installCmd = [cls.getCondaActivationCmd()]

        # Create the environment
        installCmd.append('conda create -y -n %s python=3.7;' % ENV_NAME)

        # Activate the new environment
        installCmd.append('conda activate %s;' % ENV_NAME)

        # Install downloaded code
        installCmd.extend(['conda install -y pytorch cudatoolkit=10.1 -c pytorch &&',
                           'conda install -y pandas seaborn scikit-learn &&',
                           'conda install -y -c conda-forge umap-learn jupyterlab &&',
                           'pip install ipywidgets cufflinks &&',
                           'pip install -e . &&'])

        # Flag installation finished
        installCmd.append('touch %s' % CRYODRGN_INSTALLED)

        cryodrgn_commands = [(" ".join(installCmd), CRYODRGN_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('cryodrgn', version=version,
                       url='https://github.com/zhonge/cryodrgn/archive/%s.tar.gz' % version,
                       commands=cryodrgn_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getProgram(cls, program, gpus='0'):
        """ Create cryoDRGN command line. """
        fullProgram = '%s %s && CUDA_VISIBLE_DEVICES=%s cryodrgn %s' % (
            cls.getCondaActivationCmd(), cls.getCryoDrgnEnvActivation(),
            gpus, program)

        return fullProgram

    @classmethod
    def getActiveVersion(cls, *args):
        """ Return the env name that is currently active. """
        envVar = cls.getVar(CRYODRGN_ENV_ACTIVATION)
        return envVar.split()[-1]

    @classmethod
    def IS_V03(cls):
        return cls.getActiveVersion().startswith(getCryoDrgnEnvName('0.3'))
