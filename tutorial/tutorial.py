import os
import shutil
from rdkit.Chem import AllChem

from importlib_resources import files

import silc.util as util
from silc.dock import dock
from silc.force_field import gaff2
from silc.build_molecule import binding_molecule, complex

core_smiles = "*Oc6ccc(c2ccc3ccc5c(c1ccc(O*)cc1)ccc4ccc2c3c45)cc6"
tail_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "*"
bridge_smiles = lambda nEO: "O=C(O)c1cc(c2ccc(C" + "OCC" * nEO + "*)cc2)cc(c2ccc(C" + "OCC" * nEO + "*)cc2)c1"

'''
if os.path.exists("molecules"):
    shutil.rmtree("molecules")
os.mkdir("molecules")
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('COR'), os.path.join("molecules", "COR"))
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('BRD'), os.path.join("molecules", "BRD"))
shutil.copytree(files('silc.data.tutorial.build_molecules').joinpath('TLA_TLB'), os.path.join("molecules", "TLA_TLB"))
'''

bm = binding_molecule()
bm.set_core_smiles(core_smiles, "C")
bm.set_tail_smiles(tail_smiles(3), "O")
bm.set_bridge_smiles(bridge_smiles(1), "O")
bm.set_core_num_confs_for_charge(1)
bm.set_tail_num_confs_for_charge(1) # 5
bm.set_bridge_num_confs_for_charge(1) # 5
bm.set_work_path("molecules")
bm.create_binding_motif(solvate=True)
bm.create_binding_motif(solvate=True, nmol=2, translate=[0,10,10])
bm.create_ditopic_molecule(solvate=True)
bm.create_ditopic_molecule(solvate=True, nmol=2, translate=[0,10,10])


d = dock(sf_name='ad4', receptor_name='GCDOH')
#d.prepare_ligand_from_smiles(motif_smiles(0))
#d.run(n_ligand=1)
#d.run(n_ligand=2)
d.load_results(files('silc.data.tutorial.docking').joinpath('GCDOH_2_ligand_ad4.pdbqt'), n_ligand=2)

# assign gaff2 force field
c = complex()
c.set_work_path("complexes")
c.set_binding_molecule(bm)
c.set_dock(d)
#c.create_receptor_motif_complex(n_motif=2, dock_pose_id=0)    # parrallel, "gauche"
c.create_receptor_motif_complex(n_motif=2, dock_pose_id=19, solvate=True)    # X-shape, "trans"
c.create_receptor_ditopic_complex(solvate=True)
