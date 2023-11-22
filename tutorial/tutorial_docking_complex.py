import os
from rdkit.Chem import AllChem

from importlib_resources import files

import silc.util as util
from silc.dock import dock
from silc.force_field import gaff2

motif_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "Oc6ccc(c2ccc3ccc5c(c1ccc(O" + "CCO"*nEO + "CC=C)cc1)ccc4ccc2c3c45)cc6"
diphenylpyrene = "c6ccc(c2ccc3ccc5c(c1ccccc1)ccc4ccc2c3c45)cc6"
core = util.gen_mol_from_smiles(diphenylpyrene)

d = dock(sf_name='ad4', receptor_name='GCDOH')
#d_ad4.prepare_ligand_from_smiles(motif_smiles(0))
##d_ad4.run(n_ligand=1)
#d_ad4.run(n_ligand=2)
d.load_results(files('silc.data.tutorial.docking').joinpath('GCDOH_2_ligand_ad4.pdbqt'), n_ligand=2)
#ligand0, ligand1 = d.ligand_mol(n_ligand=2, pose_id=0)    # parrallel, "gauche"
ligand0, ligand1 = d.ligand_mol(n_ligand=2, pose_id=19)    # X-shape, "trans"

expanded_core0 = util.expand_substructure(ligand0, core)
expanded_core1 = util.expand_substructure(ligand1, core)

nEO = 3
stretch = True
motif = util.gen_mol_from_smiles(motif_smiles(nEO))
motif0 = util.align_to_substructure(motif, expanded_core0, stretch=stretch)
motif = util.gen_mol_from_smiles(motif_smiles(nEO))
motif1 = util.align_to_substructure(motif, expanded_core1, stretch=stretch)

receptor = AllChem.MolFromPDBFile(str(files('silc.data.receptor').joinpath('GCDOH.pdbqt')))
receptor, motif_list = util.optimize_complex(receptor, [motif0, motif1], [expanded_core0, expanded_core1])
motif0, motif1 = motif_list
motif0_h = util.add_Hs(motif0)
motif1_h = util.add_Hs(motif1)
if not os.path.exists("docking_motif"):
    os.mkdir("docking_motif")
util.save_pdb(motif0_h, os.path.join("docking_motif", "motif0.pdb"))
util.save_pdb(motif1_h, os.path.join("docking_motif", "motif1.pdb"))
util.save_pdb(motif0_h, os.path.join("docking_motif", "motif0_nobond.pdb"), bond=False)    # save for antechamber
util.save_pdb(motif1_h, os.path.join("docking_motif", "motif1_nobond.pdb"), bond=False)    # save for antechamber


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
ff = gaff2.setup()
ff.set_work_path("docking_complex")
ff.add_receptor(receptor_name="GCDOH")
ff.add_motif(pdb_file=os.path.join("docking_motif", "motif0_nobond.pdb"),charge_file=chf)
ff.add_motif(pdb_file=os.path.join("docking_motif", "motif1_nobond.pdb"),charge_file=chf)
ff.run()
