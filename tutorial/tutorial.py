import os
from rdkit.Chem import AllChem
import meeko
from vina import Vina
import py3Dmol

from importlib_resources import files

import silc.util as util
from silc.dock import dock

### motif1: diphenylpyreneEO with OH chain end
### motif2: diphenylpyreneEO with CC=C chain end
### motif3: diphenylpyreneCNEO with CC=C chain end
### motif4: diphenylpyreneCNEO with CC=C chain end and N protonated
motifA_smiles = lambda nEO: "OCC"*nEO + "Oc6ccc(c2ccc3ccc5c(c1ccc(O" + "CCO"*nEO + ")cc1)ccc4ccc2c3c45)cc6"
motifB_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "Oc6ccc(c2ccc3ccc5c(c1ccc(O" + "CCO"*nEO + "CC=C)cc1)ccc4ccc2c3c45)cc6"
motifC_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "NCc6ccc(c2ccc3ccc5c(c1ccc(CN" + "CCO"*nEO + "CC=C)cc1)ccc4ccc2c3c45)cc6"
motifD_smiles = lambda nEO: "C=CC" + "OCC"*nEO + "[N+]Cc6ccc(c2ccc3ccc5c(c1ccc(C[N+]" + "CCO"*nEO + "CC=C)cc1)ccc4ccc2c3c45)cc6"

diphenylpyrene = "c6ccc(c2ccc3ccc5c(c1ccccc1)ccc4ccc2c3c45)cc6"
core = util.gen_mol_from_smiles(diphenylpyrene)

d = dock(sf_name='ad4', receptor_name='GCDOH')
d.load_results(files('silc.data.tutorial').joinpath('GCDOH_2_ligand_ad4.pdbqt'), n_ligand=2)
ligand0, ligand1 = d.ligand_mol(n_ligand=2, pose_id=0)    # parrallel
#ligand0, ligand1 = d.ligand_mol(n_ligand=2, pose_id=19)    # X-shape

expanded_core0 = util.expand_substructure(ligand0, core)
expanded_core1 = util.expand_substructure(ligand1, core)

nEO = 9
stretch = True
motif = util.gen_mol_from_smiles(motifB_smiles(nEO))
motif0 = util.align_to_substructure(motif, expanded_core0, stretch=stretch)
motif = util.gen_mol_from_smiles(motifB_smiles(nEO))
motif1 = util.align_to_substructure(motif, expanded_core1, stretch=stretch)

writer = AllChem.PDBWriter("target0.pdb")
writer.write(ligand0)
writer = AllChem.PDBWriter("target1.pdb")
writer.write(ligand1)

writer = AllChem.PDBWriter("motif0.pdb")
writer.write(motif0)
writer = AllChem.PDBWriter("motif1.pdb")
writer.write(motif1)


receptor = AllChem.MolFromPDBFile(str(files('silc.data.receptor').joinpath('GCDOH.pdbqt')))
complex = util.optimize_complex(receptor, [motif0, motif1], [expanded_core0, expanded_core1])
writer = AllChem.PDBWriter("complex.pdb")
writer.write(complex)
