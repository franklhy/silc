import os
from rdkit.Chem import AllChem

from importlib_resources import files

import silc.util as util
from silc.dock import dock
from silc.force_field import gaff2
from silc.build_molecules import complex

core_smiles = "*Oc6ccc(c2ccc3ccc5c(c1ccc(O*)cc1)ccc4ccc2c3c45)cc6"
tail_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "*"
bridge_smiles = lambda nEO: "O=C(O)c1cc(c2ccc(C" + "OCC" * nEO + "*)cc2)cc(c2ccc(C" + "OCC" * nEO + "*)cc2)c1"

dm = binding_molecule()
dm.set_core_smiles(core_smiles, "C")
dm.set_tail_smiles(tail_smiles(3), "O")
dm.set_bridge_smiles(bridge_smiles(2), "O")
dm.set_core_num_confs_for_charge(3)
dm.set_tail_num_confs_for_charge(3)
dm.set_bridge_num_confs_for_charge(3)
dm.set_work_path("building_moleclue")
dm.create_binding_motif(solvate=True)
dm.create_ditopic_molecule(solvate=True)

d = dock(sf_name='ad4', receptor_name='GCDOH')
#d.prepare_ligand_from_smiles(motif_smiles(0))
#d.run(n_ligand=1)
#d.run(n_ligand=2)
d.load_results(files('silc.data.tutorial.docking').joinpath('GCDOH_2_ligand_ad4.pdbqt'), n_ligand=2)


# Calculate partial charge
numConfs = 5
chg = gaff2.charge(read_only=True)    # set read_only = False (default) for new molecules
chg.set_work_path(files('silc.data.tutorial').joinpath("amber_charge"))    # set the work_path
#chg = gaff2.charge()
#chg.set_work_path("amber_charge")
chg.set_molecule_from_smile(motif_smiles(nEO))
chg.set_num_confs(numConfs)
chg.run()
chf = chg.write_charge(file_name=os.path.join("docking_motif","partial_charge.txt"))


# assign gaff2 force field
c = complex()
c.set_work_path("complex")
c.set_binding_molecule()
c.set_dock(d)
#c.create_receptor_motif_complex(n_motif=2, dock_pose_id=0)    # parrallel, "gauche"
c.create_receptor_motif_complex(n_motif=2, dock_pose_id=19)    # X-shape, "trans"
c.create_receptor_ditopic_complex()
