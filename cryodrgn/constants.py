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


V3_1_0 = "3.1.0"
V3_3_2 = "3.3.2"
V3_4_0 = "3.4.0"

VERSIONS = [V3_1_0, V3_3_2, V3_4_0]
CRYODRGN_DEFAULT_VER_NUM = V3_4_0

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

# ab initio type
AB_INITIO_HOMO = 0
AB_INITIO_HETERO = 1

# Linkage for agglomerative clustering
CLUSTER_AVERAGE = 0
CLUSTER_WARD = 1

# FlexHub program
CRYODRGN = "CryoDRGN"
