import os
from collections import OrderedDict

from rdkit.Chem import AllChem

def clean_PDBQT(pdbqt_file):
    '''
    Remove the charge and element in PDBQT file written by AutoDock
    '''
    new_l = []
    new_file = pdbqt_file[:-4].split('/')[-1] + "_ADcleaned.pdb"

    with open(pdbqt_file, "r") as f:
        for l in f:
            new_l.append(l[:66] + "\n")

    with open(new_file, "w") as f:
        f.writelines(new_l)

    return new_file


def gen_mol_from_pdb_and_smiles(pdb_file, smiles, hydrogen=False):
    '''
    Generate molecule with coordinate from PDB and correct bond order and aromaticity from SMILES.
    https://github.com/rdkit/rdkit/issues/1031

    This is useful when converting the PDBQT generated by AutoDock to sanitized molecule (with correct
    aromaticity)

    return rdkit molecule
    '''
    if hydrogen:
        mol_template = AllChem.MolFromSmiles(smiles)
        mol_template = AllChem.AddHs(mol_template)
        mol = AllChem.MolFromPDBFile(pdb_file, removeHs=False)
        mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
        return mol
    else:
        mol_template = AllChem.MolFromSmiles(smiles)
        mol = AllChem.MolFromPDBFile(pdb_file, removeHs=True)
        mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol)
        mol_with_hydrogens = AllChem.AddHs(mol)
        mol_with_hydrogens = AllChem.ConstrainedEmbed(mol_with_hydrogens, mol)
        return mol_with_hydrogens
        

def gen_mol_from_smiles(smiles, hydrogen=False):
    '''
    Generate a molecule with given smiles, and optimize the 3D structure.
    Hydrogen can be added optionally.

    Return rdkit molecule
    '''
    mol = AllChem.MolFromSmiles(smiles)
    if hydrogen:
        mol = AllChem.AddHs(mol)
    rtn = AllChem.EmbedMolecule(mol)
    if rtn == -1:    # fail to generate conformer by AllChem.EmbedMolecule
        AllChem.Compute2DCoords(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def add_Hs(mol):
    '''
    add hydrogens and optimize the structure
    '''
    mol_with_hydrogens = AllChem.AddHs(mol)
    mol_with_hydrogens = AllChem.ConstrainedEmbed(mol_with_hydrogens, mol)
    return mol_with_hydrogens


def expand_substructure(mol, core, expand_iteration=1, only_include=None):
    '''
    Find the substructure from mol that matches core, then expand it to n bonded neighbors.
    n is set by expand_iteration
    only_include means only include a list of element. For example, when only_include = None, all elements are included;
    when only_include = ['H', 'C'], only hydrogen and carbon atoms are consider valid neighbor atoms.

    Return rdkit molecule
    '''
    if expand_iteration == -1:
        return AllChem.Mol(core)

    if expand_iteration < -1:
        raise RuntimeError("cannot have negative iteration number")

    if not mol.HasSubstructMatch(core):
        raise RuntimeError("No match found")

    match = mol.GetSubstructMatches(core)
    if len(match) > 1:
        raise RuntimeError("More than one substructure matched. Total %d mathces." % len(match))

    # Add the atoms from the matched substructure and their 1st neighbors
    core_atom_indices = list(match[0])
    keep_atom_indices = list(match[0])
    if expand_iteration > 0:
        for atom_idx in core_atom_indices:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in keep_atom_indices:
                    neigh_atom = mol.GetAtomWithIdx(neighbor_idx)
                    if only_include is None or AllChem.GetPeriodicTable().GetElementSymbol(neigh_atom.GetAtomicNum()) in only_include:
                        keep_atom_indices.append(neighbor_idx)

    # Initialize a new empty molecule to store the expanded substructure
    expanded_core = AllChem.RWMol(mol)
    expanded_core.BeginBatchEdit()
    for atom in expanded_core.GetAtoms():
        atom_idx = atom.GetIdx()
        if atom_idx not in keep_atom_indices:
            expanded_core.RemoveAtom(atom_idx)
    expanded_core.CommitBatchEdit()

    expand_iteration -= 1
    return expand_substructure(mol, expanded_core, expand_iteration)


def align_to_substructure(mol, core, stretch=False, stretch_ratio_min=0.5, stretch_ratio_max=0.7):
    '''
    Align a molecule to a core substructure.
    Optionally stretch the molecule between chain ends to obtain extended conformation.

    core must have at least one conformer
    '''
    m = AllChem.Mol(mol)
    match = m.GetSubstructMatch(core)
    m = AllChem.ConstrainedEmbed(m, core, randomseed=123, maxAttempts=100000)
    mp = AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94')
    ff = AllChem.MMFFGetMoleculeForceField(m, mp)

    if stretch:
        # Find chain end
        chain_end_atom_idx = []

        # Iterate through the atoms in the molecule
        for atom in m.GetAtoms():
            # Check if the atom is connected to only one other atom (chain end)
            if atom.GetDegree() == 1:
                chain_end_atom_idx.append(atom.GetIdx())

        # calculate the topological distance matrix, to get the number of bonds between chain ends
        dm = AllChem.GetDistanceMatrix(m)
        # for each pair of chain ends, add a distance constrain
        for atomi in chain_end_atom_idx:
            for atomj in chain_end_atom_idx:
                if atomi > atomj:
                    distance = dm[atomi,atomj]
                    ### force chain ends to be far away from each other
                    ff.MMFFAddDistanceConstraint(atomi, atomj, stretch, distance*stretch_ratio_min, distance*stretch_ratio_max, 0.5)
        # constrain the motion of the core part
        for atom_idx in match:
            ff.MMFFAddPositionConstraint(atom_idx, 0.01,200)

        ff.Minimize(maxIts=10000)

    return m


def generate_multiple_conformers(mol, num_confs, prune_rms_thresh=1.0):
    '''
    Generate multiple conformers for a molecule
    '''
    ps = AllChem.ETKDG()
    ps.pruneRmsThresh = prune_rms_thresh
    #ps.optimizerForceTol = 0.01
    ps.useRandomCoords = True
    #ps.maxIterations = 1000
    m = AllChem.Mol(mol)
    cids = AllChem.EmbedMultipleConfs(m, numConfs=num_confs, params=ps)    # generate more confs in case some will be pruned
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(m, confId=cid)
    return m


def read_conformer_from_pdb(mol, pdb_file):
    '''
    Read the coordination of a conformer for a molecule from a pdb file
    '''
    mol_tmp = AllChem.MolFromPDBFile(pdb_file, sanitize=False, removeHs=False)
    mol.AddConformer(mol_tmp.GetConformer(), assignId=True)
    return mol


def optimize_complex(receptor, ligand_list, ligand_core_list, max_its=10000):
    '''
    Use MMFF to optimize the complex structure, while keeping the receptor and the core of ligand rigid.
    '''
    mols = [receptor] + ligand_list
    complex = AllChem.CombineMols(mols[0], mols[1])
    for i in range(2,len(mols)):
        complex = AllChem.CombineMols(complex, mols[i])
    mp = AllChem.MMFFGetMoleculeProperties(complex, mmffVariant='MMFF94')
    ff = AllChem.MMFFGetMoleculeForceField(complex, mp)

    # fix the position of receptor
    match = complex.GetSubstructMatch(receptor)
    for atom_idx in match:
        ff.MMFFAddPositionConstraint(atom_idx, 0.01, 200)
    # fix the position of ligand core
    for i in range(len(ligand_core_list)):
        matches = complex.GetSubstructMatches(ligand_core_list[i])
        for match in matches:
            for atom_idx in match:
                ff.MMFFAddPositionConstraint(atom_idx, 0.01, 200)
    ff.Minimize(maxIts=max_its)

    frags = AllChem.GetMolFrags(complex, asMols=True)
    receptor_new = frags[0]
    if len(frags) > 2:
        ligand_list_new = frags[1:]
    elif len(frags) == 2:
        ligand_list_new = [frags[1]]
    else:
        raise RuntimeError("Something is wrong")

    ### check if the fragments are assigned correctly
    receptor_check = receptor_new.HasSubstructMatch(receptor)
    ligand_check = []
    for ligand_new, ligand in zip(ligand_list_new, ligand_list):
        ligand_check.append(ligand_new.HasSubstructMatch(ligand))

    if not receptor_check or not all(ligand_check):
        raise RuntimeError("Something is wrong")

    return receptor_new, ligand_list_new


def set_resname_by_substructure(mol, core, res_name, res_num, res_num_mode='='):
    '''
    Find matching substructure (core) in a given molecule (mol), and rename its residue name to res_name and set new residue number by res_num.

    res_name should has two letters, the third letter is save for index (A, B, C, ...) in case there are multiple matched substructure.
    res_name will be converted to uppercase if lowercase is provided

    res_num should be an integer >= 1, it will be the initial residue number for the first matched substructure, and then it will increase
    by 1 each time another matched substrcture is found

    res_num_mode:
        '=': assign res_num as residue number
        '+':    add res_num to the current residue number
    '''
    m = AllChem.MolFromPDBBlock(AllChem.MolToPDBBlock(mol), removeHs=False)    # generate residue name using MolToPDBBlock
    c = AllChem.Mol(core)
    c = AllChem.RemoveAllHs(c)
    core_h = expand_substructure(m, c, expand_iteration=1, only_include=['H'])

    if not m.HasSubstructMatch(core_h):
        raise RuntimeError("No match found")

    if len(res_name) != 2:
        raise RuntimeError('resname should has two letters.')

    match = m.GetSubstructMatches(core_h)
    next_alpha = lambda s: chr((ord(s.upper())+1 - 65) % 26 + 65)    # return the next uppercase alphabet, e.g. next_alpha('A') = 'Z', next_alpha('Z') = 'A'
    res_idx = 'A'
    for i in range(len(match)):
        atom_indices = list(match[i])
        for atom_idx in atom_indices:
            atom = m.GetAtomWithIdx(atom_idx)
            mi = atom.GetMonomerInfo()
            mi.SetResidueName(res_name[0].upper() + res_name[1].upper() + res_idx)
            if res_num_mode == '=':
                mi.SetResidueNumber(res_num)
            elif res_num_mode == '+':
                tmp = mi.GetResidueNumber()
                mi.SetResidueNumber(res_num + tmp)
            else:
                raise RuntimeError("invalid res_num_mode")
        res_idx = next_alpha(res_idx)
        res_num += 1
    return m


def save_pdb(mol, name, conf_id=-1, bond=True):
    '''
    Save molecule to pdb file.
    '''
    if len(name) < 5 or name[-4:] != ".pdb":
        name = "%s.pdb" % name
    if bond:
        AllChem.MolToPDBFile(mol, name, confId=conf_id, flavor=4)
    else:
        ### flavor=2 still generate a few 'CONECT' lines
        lines = AllChem.MolToPDBBlock(mol, confId=conf_id).split('\n')
        with open(name, "w") as f:
            for line in lines:
                if not line.startswith("CONECT"):
                    f.write(line + '\n')


def read_mol2_charge(mol2file):
    '''
    Read mol2 file generated by antechamber and get the partial charge.
    '''
    partial_charge = OrderedDict()
    with open(mol2file, "r") as f:
        readflag = False
        for line in f:
            if line.startswith("@<TRIPOS>BOND"):
                break
            elif readflag:
                ws = line.split()
                partial_charge[(int(ws[0]), ws[1])] = float(ws[8])
            elif line.startswith("@<TRIPOS>ATOM"):
                readflag = True
    return partial_charge


def read_amber_mol2(mol2file):
    '''
    Read AMBER mol2 file and generate rdkit molecule.
    '''
    with open(mol2file, "r") as f:
        block = f.read()
    block_list = list(block)
    for i in range(len(block_list)):
        if block_list[i] == "h" and block_list[i+1] != '':
            block_list[i] = "H"
            block_list[i+1] = ""
        elif block_list[i] == "c" and block_list[i+1] != '':
            block_list[i] = "C"
            block_list[i+1] = ""
        elif block_list[i] == "n" and block_list[i+1] != '':
            block_list[i] = "N"
            block_list[i+1] = ""
        elif block_list[i] == "o" and block_list[i+1] != '':
            block_list[i] = "O"
            block_list[i+1] = ""
        elif block_list[i] == "f" and block_list[i+1] != '':
            block_list[i] = "F"
            block_list[i+1] = ""
    newblock = "".join(block_list)
    mol = AllChem.MolFromMol2Block(newblock, removeHs=False)
    return mol


def replace_dummy(smiles, new=None, replace_mass_label=False):
    '''
    Replace dummmy atoms in smiles. At most two dummy atoms allowed in a smiles string.
    Dummy atoms should be denoted as "[1*]" and "[2*]".
    new should either be None (no change) or a list of character (including the empty one: "") with the same number of dummy atoms.
    '''
    if smiles.count("*") > 3:
            raise RuntimeError("Too many dummy atoms. Should be less than 3.")
    if new is None:
        return smiles
    elif not replace_mass_label:
        if smiles.count("*") == 1:
            return smiles.replace("*", new[0])
        elif smiles.count("*") == 2:
            loc = smiles.find("[1*]")
            smiles = smiles[:loc+2] + new[0] + smiles[loc+3:]
            loc = smiles.find("[2*]")
            smiles = smiles[:loc+2] + new[1] + smiles[loc+3:]
            return smiles
    elif replace_mass_label:
        if smiles.count("*") == 1:
            loc = smiles.find("[1*]")
            smiles = smiles[:loc] + new[0] + smiles[loc+4:]
            return smiles
        elif smiles.count("*") == 2:
            ### replace the first dummy atom
            loc = smiles.find("[1*]")
            smiles = smiles[:loc] + new[0] + smiles[loc+4:]
            ### replace the second dummy atom
            loc = smiles.find("[2*]")
            smiles = smiles[:loc] + new[1] + smiles[loc+4:]
            return smiles


def check_leap_log(leap_log):
    '''
    Check the status of amber LEaP run
    
    Return: A list of [total number of errors, total number of warnings, total number of notes].
    '''
    head = "Exiting LEaP:"
    status = []
    with open(leap_log, "r") as f:
        for line in f:
            if line.startswith(head):
                part = line[len(head):]
                part = part.replace('.', ';')
                w = part.split(';')
                num_error = int(w[0].split()[2])
                num_warning = int(w[1].split()[2])
                num_note = int(w[2].split()[2])
                status.append([num_error, num_warning, num_note])
    total = [0, 0, 0]
    for s in status:
        for i in range(3):
            total[i] += s[i]
    return total