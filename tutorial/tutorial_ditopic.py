import os
import shutil
import subprocess
from rdkit.Chem import AllChem
import py3Dmol

from importlib_resources import files

import silc.util as util
from silc.build_molecule import binding_molecule
from silc.force_field import gaff2

core1_smiles = "*OC(=O)c1ccc2cc(ccc2c1)C(=O)O*"
core2_smiles = "*Oc6ccc(c2ccc3ccc5c(c1ccc(O*)cc1)ccc4ccc2c3c45)cc6"
core3_smiles = "*C(c2ccc3ccc5c(C*)ccc4ccc2c3c45)"

tail1_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "*"
tail2_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "N" + "*"
tail3_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "[N+]" + "*"

bridge1_smiles = lambda nEO: "O=C(O)c1cc(c2ccc(C" + "OCC" * nEO + "*)cc2)cc(c2ccc(C" + "OCC" * nEO + "*)cc2)c1"
bridge2_smiles = lambda nEO: "CC1(C)c2ccc(C" + "OCC" * nEO + "*)cc2Oc2cc(C" + "OCC" * nEO + "*)ccc21"


if os.path.exists("ditopic_solv"):
    shutil.rmtree("ditopic_solv")
os.mkdir("ditopic_solv")
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('COR'), os.path.join("ditopic_solv", "COR"))
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('BRD'), os.path.join("ditopic_solv", "BRD"))
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('TLA_TLB'), os.path.join("ditopic_solv", "TLA_TLB"))

dm = binding_molecule()
dm.set_core_smiles(core2_smiles, "C")
dm.set_tail_smiles(tail1_smiles(3), "O")
dm.set_bridge_smiles(bridge1_smiles(1), "O")
dm.set_core_num_confs_for_charge(1)
dm.set_tail_num_confs_for_charge(5)
dm.set_bridge_num_confs_for_charge(5)
dm.set_work_path("ditopic_solv")
#dm.clear_work_path()    # this will remove the work_path, it is only for building new molecules from new segments
dm.create_ditopic_molecule(solvate=True)
