# **************************************************************************
# *
# * Authors:     Eduardo García (eduardo.garcia@cnb.csic.es)     [2]
# *
# * [1] MRC Laboratory of Molecular Biology, MRC-LMB
# * [2] Unidad de  Biocomputacion, Centro Nacional de Biotecnologia, CSIC (CNB-CSIC)
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

import numpy as np
import torch
import torch.nn as nn
from cryodrgn import config
from cryodrgn.commands.eval_vol import reset_origin, postprocess_vol
from cryodrgn.models import HetOnlyVAE, load_decoder
from cryodrgn.source import write_mrc

class HeterogeneityProgramInterface:
    def __init__(self, _path_template: str, _program_loading_params: dict):
        self.model = self.prepare_heterogeneity_program(**_program_loading_params)
        self.path_template = _path_template

    def prepare_heterogeneity_program(self, **kwargs) -> object:
        gpu_id = kwargs.pop("gpu_id", None)
        config = kwargs.pop("config", None)
        load = kwargs.pop("load", None)
        self.device = "cpu" if gpu_id is None else 'cuda:' + str(int(gpu_id))

        cfg = config.load(config)

        D = cfg["lattice_args"]["D"]
        zdim = cfg["model_args"]["zdim"]
        self.norm = [float(x) for x in cfg["dataset_args"]["norm"]]
        self.Apix = cfg["model_args"]["Apix"]

        if "players" in cfg["model_args"]:
            model, self.lattice = HetOnlyVAE.load(cfg, load, device=self.device)
            decoder = model.decoder
        else:
            decoder, self.lattice = load_decoder(cfg, load, device=self.device)
        decoder.eval()
        return decoder

    def decode_state_from_latent(self, latent: np.array) -> None:
        latent = torch.from_numpy(latent.astype(np.float32)).to(self.device)
        for idx, zz in enumerate(latent):
            self.model.eval_volume(self.lattice.coords, self.lattice.D, self.lattice.extent, self.norm, zz)
        out_mrc = "{}/{}{:03d}.mrc".format(args.o, args.prefix, i)
        write_mrc(out_mrc, np.array(vol.cpu()).astype(np.float32), Apix=self.Apix)