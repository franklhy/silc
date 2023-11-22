import os
import shutil
import subprocess
from collections import OrderedDict
from importlib_resources import files

import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem

from .. import util

class charge:
    '''
    Handle partial charge for molecules.
    '''

    def __init__(self, method='bcc', read_only=False, antechamber_status=1, remove_antechamber_intermediate_files=True):
        '''
        read_only: read only mode, only read partial charge from work_path where other partial charge calculation results are saved
        '''
        self.method = method
        self.read_only = read_only
        self.antechamber_status = antechamber_status
        if remove_antechamber_intermediate_files:
            self.raif = 'y'
        else:
            self.raif = 'n'
        self.mol = None
        self.mol_charge = None
        self.work_path = None
        self.num_confs = 1
        self.charge = OrderedDict()    # key: (atom index, atom name)
        self.charge_avgstd = OrderedDict()    # key: (atom index, atom name)


    def set_work_path(self, work_path="amber_charge"):
        work_path = str(work_path)
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def set_molecule(self, mol):
        self.mol = AllChem.Mol(mol)
        self.mol = util.add_Hs(self.mol)


    def set_molecule_from_smile(self, smiles):
        self.mol = util.gen_mol_from_smiles(smiles, hydrogen=True)


    def set_num_confs(self, num_confs):
        self.num_confs = num_confs


    def run(self):
        '''
        Either read the conformers' atomic coordinates and charges generated previously (read only mode),
        or generate new conformers and calculate partial charge form them
        '''
        if self.read_only:
            self.read_multiconformer_crd()
            self.read_multiconformer_charge()
            self.symmetrize_charge()
        else:
            self.calculate_partial_charge()
            self.read_multiconformer_charge()
            self.symmetrize_charge()


    def calculate_partial_charge(self):
        '''
        Use antechamber to calculate partial charge for atoms of a molecule, average over several gas phase conformations. For larger molecules, it will take very long.
        Right now, only consider GAFF2 force field and AM1-BCC charge.
        '''
        if self.read_only:
            raise RuntimeError("Cannot calculate partial charge in read only mode.")
        cwd = os.getcwd()

        # generate molecule and conformers
        prune_rms_thresh = 1.0
        while self.mol.GetNumConformers() != self.num_confs:
            self.mol = util.generate_multiple_conformers(self.mol, self.num_confs, prune_rms_thresh=prune_rms_thresh)
            prune_rms_thresh *= 0.5
        print("Generated %d conformers. prune_rms_thresh = %f" % (self.num_confs, prune_rms_thresh))

        # run antechamber to calculate partial charge
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)
        os.makedirs(self.work_path)
        os.chdir(self.work_path)
        for cid in range(self.num_confs):
            util.save_pdb(self.mol, "confId%d.pdb" % cid, conf_id=cid, bond=False)
            subprocess.run(["antechamber", "-i", "confId%d.pdb" % cid, "-fi", "pdb", "-o", "confId%d.mol2" % cid, "-fo", "mol2", "-c", "bcc", "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
        os.chdir(cwd)


    def read_multiconformer_crd(self):
        '''
        Read the coordinate of multiple conformers
        '''
        self.mol.RemoveAllConformers()
        for cid in range(self.num_confs):
            util.read_conformer_from_pdb(self.mol, os.path.join(self.work_path, "confId%d.pdb" % cid))


    def read_multiconformer_charge(self):
        '''
        Read atomic partial charge of multiple conformers
        '''
        # gather partial charge from all conformers
        self.charge = util.read_mol2_charge(os.path.join(self.work_path, "confId0.mol2"))
        for key in self.charge.keys():
            self.charge[key] = [self.charge[key]]
        for cid in range(1, self.num_confs):
            charge_tmp = util.read_mol2_charge(os.path.join(self.work_path, "confId%d.mol2" % cid))
            for key in charge_tmp.keys():
                self.charge[key].append(charge_tmp[key])

        # average partial charge among different conformers
        charge_sum = 0
        for key in self.charge.keys():
            self.charge_avgstd[key] = [np.mean(self.charge[key]), np.std(self.charge[key])]
            charge_sum += self.charge_avgstd[key][0]
        for key in self.charge.keys():
            self.charge_avgstd[key][0] -= charge_sum / len(self.charge_avgstd)


    def symmetrize_charge(self):
        '''
        Check and make sure that partial charges on equivalent atoms are the same
        Note: the '-eq' option of antechamber should have already handeled it.
        '''
        sym_chg = {}
        new_sym_chg = {}
        sym = list(AllChem.CanonicalRankAtoms(self.mol, breakTies=False))  # get symmetry class for each atom
        for i,(k,v) in enumerate(self.charge_avgstd.items()):
            if sym[i] in sym_chg:
                sym_chg[sym[i]].append(v[0])
            else:
                sym_chg[sym[i]] = [v[0]]
        for k,v in sym_chg.items():
            if len(v) > 1:
                uniq = np.unique(v)
                if len(uniq) == 1:
                    continue
                else:
                    print("Warning: Partial charges of equivalent atoms are different. Enforce charge symmetrization.")
                    avg = np.mean(v)
                    new_sym_chg[k] = avg
                    dev = (np.array(v) - avg) / avg
                    if np.max(dev) > 0.01:
                        print("Warning: Partial charges of equivalent atoms differs more by 1%%.")
        for i,(k,v) in enumerate(self.charge_avgstd.items()):
            if sym[i] in new_sym_chg:
                self.charge_avgstd[k] = new_sym_chg[sym[i]]


    def write_charge(self, file_name="partial_charge.txt"):
        # write out the partial charge, using the same format by antechamber
        count = 0
        with open(file_name, "w") as f:
            for key, value in self.charge_avgstd.items():
                f.write("%10.6lf" % value[0])
                count += 1
                if count == 8:
                    count = 0
                    f.write("\n")
        file_abspath = os.path.join(os.path.abspath(os.getcwd()), file_name)
        return file_abspath


    def plot_charge(self):
        fig,ax = plt.subplots(1)
        colors = {'H': 'grey', 'C': 'black', 'O': 'red', 'N': 'blue'}
        for key in self.charge.keys():
            if key[1][0] in list(colors.keys()):
                color = colors[key[1][0]]
            else:
                color = 'green'
            ax.plot([np.mean(self.charge[key])] * len(self.charge[key]), self.charge[key], ls='none', marker='o', color=color, fillstyle='none')
        minq = 100
        maxq = -100
        for k,v in self.charge.items():
            if minq > np.min(v):
                minq = np.min(v)
            if maxq < np.max(v):
                maxq = np.max(v)
        xs = np.linspace(minq-0.1, maxq+0.1, 10)
        ys = xs * 1.0
        ax.plot(xs, ys, 'k--')

        ax.set_xlabel('average atomic charge')
        ax.set_ylabel('conformer atomic charge')

        from matplotlib.lines import Line2D
        legend_elements = []
        for key, value in colors.items():
            legend_elements.append(Line2D([0], [0], color=value, ls='none', marker='o', fillstyle='none', label=key))
        legend_elements.append(Line2D([0], [0], color='green', ls='none', marker='o', fillstyle='none', label='Other'))
        ax.legend(handles=legend_elements)


    def draw_mol_with_charge(self, image_size=(800, 800)):
        self.mol_charge = AllChem.MolFromPDBBlock(AllChem.MolToPDBBlock(self.mol), removeHs=False)
        AllChem.Compute2DCoords(self.mol_charge)
        for atom in self.mol_charge.GetAtoms():
            aname = atom.GetMonomerInfo().GetName().strip()
            aidx = atom.GetIdx() + 1
            lbl = "%f" % self.charge_avgstd[(aidx,aname)][0]
            atom.SetProp('atomNote', lbl)
            atom.SetProp("_TriposPartialCharge", lbl)
        img = AllChem.Draw.MolsToGridImage([self.mol_charge], molsPerRow=1, subImgSize=image_size)
        return img


class setup:
    '''
    Use antechamber and other amber tools to setup GAFF2 force fields as well as amber input topology and coordinate files.
    '''
    def __init__(self, antechamber_status=1, remove_antechamber_intermediate_files=True):
        self.n_motif = 0
        self.antechamber_status = antechamber_status
        if remove_antechamber_intermediate_files:
            self.raif = 'y'
        else:
            self.raif = 'n'
        self.work_path = None
        self.motif_pdb = []
        self.motif_charge_file = []
        self.receptor = None


    def set_work_path(self, work_path="amber_forcefield"):
        work_path = str(work_path)
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def add_motif(self, pdb_file, charge_file):
        pdb_abs = os.path.abspath(pdb_file)
        charge_abs = os.path.abspath(charge_file)
        self.motif_pdb.append(pdb_abs)
        self.motif_charge_file.append(charge_abs)
        self.n_motif += 1


    def add_receptor(self, receptor_name="GCDOH"):
        self.receptor = receptor_name


    def run(self):
        cwd = os.getcwd()
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)
        os.makedirs(self.work_path)
        os.chdir(self.work_path)

        for i in range(self.n_motif):
            subprocess.run(["antechamber", "-i", self.motif_pdb[i], "-fi", "pdb", "-o", "motif%d.mol2" % i, "-fo", "mol2", "-c", "rc", "-cf", self.motif_charge_file[i], "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
            subprocess.run(["parmchk2", "-i", "motif%d.mol2" % i, "-f", "mol2", "-o", "motif%d.frcmod" % i, "-s", "gaff2"])

        with open("tleap.in", "w") as f:
            f.write("source leaprc.gaff2\n\n")

            for i in range(self.n_motif):
                f.write("motif%d = loadmol2 motif%d.mol2\n"  % (i, i))
                f.write("check motif%d\n" % i)
                f.write("loadamberparams motif%d.frcmod\n\n" % i)
                #f.write("saveoff motif%d motif%.lib\n" % (i, i))

            if self.receptor is not None:
                shutil.copy2(files('silc.data.receptor').joinpath("glycam04.dat"), ".")
                shutil.copy2(files('silc.data.receptor').joinpath("frcmod.q4md"), ".")
                shutil.copy2(files('silc.data.receptor').joinpath("%s.pdb" % self.receptor), ".")
                shutil.copy2(files('silc.data.receptor').joinpath("q4md-CD.off"), ".")

                f.write("loadamberparams glycam04.dat\n")
                f.write("loadamberparams frcmod.q4md\n")
                f.write("%s = loadPDB %s.pdb\n" % (self.receptor, self.receptor))
                f.write("loadoff q4md-CD.off\n\n")

                f.write("complex = combine {%s" % self.receptor)
                for i in range(self.n_motif):
                    f.write(" motif%d" % i)
                f.write("}\n")
                f.write("check complex\n\n")

                f.write("saveamberparm complex complex.prmtop complex.rst7\n")
                f.write("savepdb complex complex.pdb\n\n")

                f.write("loadoff solvents.lib\n")
                f.write("solvateOct complex TIP3PBOX 14.0\n")
                f.write("saveamberparm complex complex_solv.prmtop complex_solv.rst7\n")
                f.write("savepdb complex complex_solv.pdb\n\n")

            f.write("quit\n")
        subprocess.run(["tleap", "-f", "tleap.in"])

        os.chdir(cwd)


class residue:
    def __init__(self, charge_method='bcc', read_only=False, antechamber_status=1, remove_antechamber_intermediate_files=True):
        '''
        read_only: read only mode, only read partial charge from work_path where other partial charge calculation results are saved
        '''
        self.charge_method = charge_method
        self.read_only = read_only
        self.antechamber_status = antechamber_status
        if remove_antechamber_intermediate_files:
            self.raif = 'y'
        else:
            self.raif = 'n'
        self.resname = None
        self.resname_head = None    # residue name for head capping residue
        self.resname_tail = None    # residue name for tail capping residue
        self.replace_dummy_with = None
        self.smiles_with_dummy = None
        self.smiles = None
        self.work_path = None
        self.num_confs = 5

    def set_resname(self, resname):
        '''
        Set residue name
        '''
        if len(resname) != 3 or not isinstance(resname, str):
            raise RuntimeError("Need three letters for residue name.")
        self.resname = resname.upper()


    def set_resname_head(self, resname):
        '''
        Set residue name for head capping residue
        '''
        if len(resname) != 3 or not isinstance(resname, str):
            raise RuntimeError("Need three letters for residue name.")
        self.resname_head = resname.upper()


    def set_resname_tail(self, resname):
        '''
        Set residue name for tail capping residue
        '''
        if len(resname) != 3 or not isinstance(resname, str):
            raise RuntimeError("Need three letters for residue name.")
        self.resname_tail = resname.upper()


    def set_dummy_replacement(self, dummy_replacement):
        if isinstance(dummy_replacement, str) and len(dummy_replacement) == 1:
            self.replace_dummy_with = dummy_replacement
        else:
            raise RuntimeError("Invalide dummy replacement atom.")


    def set_smiles(self, smiles):
        if self.replace_dummy_with is None:
            raise RuntimeError("Please set dummy replacement first using set_dummy_replacement().")
        if smiles.count("*") == 0:
            raise RuntimeError("Need dummy atoms to perform reaction.")
        self.smiles_with_dummy = smiles
        self.smiles = smiles.replace("*", self.replace_dummy_with)


    def set_work_path(self, work_path):
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def set_num_confs(self, num_confs):
        '''
        Number of conformers for partial charge calculation.
        '''
        self.num_confs = num_confs


    def run(self):
        cwd = os.getcwd()
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)
        os.makedirs(self.work_path)
        os.chdir(self.work_path)

        with open("smiles.txt", "w") as f:
            f.write("%s\n" % self.smiles_with_dummy)
            f.write("*:%s\n" % self.replace_dummy_with)

        # calculate partial charge
        chg = charge(method=self.charge_method)
        chg.set_work_path(os.path.join(self.work_path, "amber_charge"))
        chg.set_molecule_from_smile(self.smiles)
        chg.set_num_confs(self.num_confs)
        chg.run()
        chf = chg.write_charge(file_name=os.path.join("amber_charge", "partial_charge.txt"))

        # prepare mainchain file
        mol = AllChem.MolFromSmiles(self.smiles_with_dummy)
        dummy = []
        for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetMass() == 0:
                dummy.append(i)

        mol = AllChem.MolFromPDBBlock(AllChem.MolToPDBBlock(chg.mol), removeHs=False)
        remove = []
        head = []
        for atom_idx in dummy:
            atom = mol.GetAtomWithIdx(atom_idx)
            remove.append(atom.GetMonomerInfo().GetName().strip())
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                neigh_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neigh_atom.GetAtomicNum() == 1:
                    remove.append(neigh_atom.GetMonomerInfo().GetName().strip())
                elif neigh_atom.GetAtomicNum() in [6,7,8]:
                    head.append(neigh_atom.GetMonomerInfo().GetName().strip())

        fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.smiles_with_dummy))    # formal charge

        if len(head) == 2:
            with open("mainchain.mol", "w") as f:
                f.write("HEAD_NAME %s\n" % head[0])
                f.write("TAIL_NAME %s\n" % head[1])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
        elif len(head) == 1:
            ### head capping residue
            with open("mainchain_head.mol", "w") as f:
                f.write("TAIL_NAME %s\n" % head[0])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
            ### tail capping residue
            with open("mainchain_tail.mol", "w") as f:
                f.write("HEAD_NAME %s\n" % head[0])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
        else:
            raise RuntimeError("cannot create a residue for the molecule given its SMILES: %s" % self.smiles_with_dummy)


        # prepare residues with amber tools
        subprocess.run(["antechamber", "-i", os.path.join("amber_charge", "confId0.mol2"), "-fi", "mol2", "-o", "molecule.ac", "-fo", "ac", "-c", "rc", "-cf", chf, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
        if len(head) == 2:
            if self.resname is None:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule.prepi", "-f", "prepi", "-m", "mainchain.mol", "-rn", "RES", "-rf", "molecule.res"])
            else:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule.prepi", "-f", "prepi", "-m", "mainchain.mol", "-rn", self.resname, "-rf", "molecule.res"])
        elif len(head) == 1:
            if self.resname_head is None:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_head.prepi", "-f", "prepi", "-m", "mainchain_head.mol", "-rn", "RSH", "-rf", "molecule_head.res"])
            else:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_head.prepi", "-f", "prepi", "-m", "mainchain_head.mol", "-rn", self.resname_head, "-rf", "molecule_head.res"])

            if self.resname_tail is None:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_tail.prepi", "-f", "prepi", "-m", "mainchain_tail.mol", "-rn", "RST", "-rf", "molecule_tail.res"])
            else:
                subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_tail.prepi", "-f", "prepi", "-m", "mainchain_tail.mol", "-rn", self.resname_tail, "-rf", "molecule_tail.res"])

        os.chdir(cwd)
