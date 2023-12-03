from importlib_resources import files

from silc.dock import dock
from silc.build_molecule import binding_molecule, complex

core_smiles = "[1*]Oc6ccc(c2ccc3ccc5c(c1ccc(O[2*])cc1)ccc4ccc2c3c45)cc6"
tail_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "[1*]"
bridge_smiles = lambda nEO: "O=C(O)c1cc(c2ccc(C" + "OCC" * nEO + "[1*])cc2)cc(c2ccc(C" + "OCC" * nEO + "[2*])cc2)c1"

bm = binding_molecule()
bm.set_core_smiles(core_smiles, ["C", "C"])
bm.set_tail_smiles(tail_smiles(3), ["O"])
bm.set_bridge_smiles(bridge_smiles(1), ["O", "O"])
bm.set_core_num_confs_for_charge(5)
bm.set_tail_num_confs_for_charge(5)
bm.set_bridge_num_confs_for_charge(5)
bm.set_work_path("molecules")
bm.clear_work_path()
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
c.clear_work_path()
c.set_binding_molecule(bm)
c.set_dock(d)
#c.create_receptor_motif_complex(n_motif=2, dock_pose_id=0)    # parrallel, "gauche"
c.create_receptor_motif_complex(n_motif=2, dock_pose_id=19, solvate=True)    # X-shape, "trans"
c.create_receptor_ditopic_complex(solvate=True)
