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


V1_0_0 = "1.0.0"
V1_1_0 = "1.1.0"

VERSIONS = [V1_0_0, V1_1_0]
CRYODRGN_DEFAULT_VER_NUM = VERSIONS[-1]

DEFAULT_ENV_NAME = getCryoDrgnEnvName(CRYODRGN_DEFAULT_VER_NUM)
DEFAULT_ACTIVATION_CMD = 'conda activate ' + DEFAULT_ENV_NAME
CRYODRGN_ENV_ACTIVATION = 'CRYODRGN_ENV_ACTIVATION'

# Viewer constants
EPOCH_LAST = 0
EPOCH_SELECTION = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

# extra metadata attrs
Z_VALUES = "_cryodrgnZValues"
WEIGHTS = "_cryodrgnWeights"
CONFIG = "_cryodrgnConfig"
