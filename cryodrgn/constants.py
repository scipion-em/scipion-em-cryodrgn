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


def getCryoDrgnEnvName(version):
    return "cryodrgn-%s" % version

def getCryoDrgnVersion():
    """ Return the currently used version of cryodrgn. """
    # FIXME: Read the real used version from somewhere (config, environ, etc)
    return CRYODRGN_DEFAULT_VER_NUM

def isCryoDrgnGT033():
    """ Return True if current version is greater that 0.3.3b. """
    # FIXME: Change it for more general solution for future versions
    return getCryoDrgnVersion() == V0_3_3b


V0_3_0b = "0.3.0b"
V0_3_1 = "0.3.1"
V0_3_2b = "0.3.2b"
V0_3_2 = "0.3.2"
V0_3_3b = "0.3.3b"

VERSIONS = [V0_3_1, V0_3_2, V0_3_3b]
# They set release tag to 'v0.3.2' instead of '0.3.2' since that version
VERSION_PREFIX = {
    V0_3_2: 'v'
}
CRYODRGN_DEFAULT_VER_NUM = VERSIONS[-1]

DEFAULT_ENV_NAME = getCryoDrgnEnvName(CRYODRGN_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
CRYODRGN_ENV_ACTIVATION = 'CRYODRGN_ENV_ACTIVATION'

# Viewer constants
EPOCH_LAST = 0
EPOCH_SELECTION = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1
