import os
import shutil
import subprocess
from collections import OrderedDict
from importlib_resources import files

import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from .. import util

class charge:
    '''
    Handle partial charge for molecules.
    '''

    def __init__(self, charge_method='bcc', nproc=1, read_only=False, antechamber_status=1, remove_antechamber_intermediate_files=True):
        '''
        read_only: read only mode, only read partial charge from work_path where other partial charge calculation results are saved
        '''
        self.charge_method = charge_method
        self.nproc = nproc
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
        self.formal_charge = 0
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
        self.formal_charge = AllChem.GetFormalCharge(AllChem.MolFromSmiles(smiles))


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
            if self.charge_method == "bcc":
                subprocess.run(["antechamber", "-i", "confId%d.pdb" % cid, "-fi", "pdb", "-o", "confId%d.mol2" % cid, "-fo", "mol2", "-c", "bcc", "-nc", "%d" % self.formal_charge, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
            elif self.charge_method == "resp":
                if shutil.which("g16") is not None:
                    subprocess.run(["antechamber", "-i", "confId%d.pdb" % cid, "-fi", "pdb", "-o", "confId%d.gjf" % cid, "-fo", "gcrt", "-nc", "%d" % self.formal_charge, "-gn", "%%NProcShared=%d" % self.nproc, "-s", "%d" % self.antechamber_status, "-pf", self.raif])
                    print("Running Gaussian 16 to calcualte ESP charge.")
                    subprocess.run(["g16","confId%d.gjf" % cid])
                #elif shutil.which("g09") is not None:
                    # NOTE: During testing, Gaussian 09 does not generate the gesp file as desired, and the resp charge fitting will fail when running the "antechamber -c resp"
                    #subprocess.run(["antechamber", "-i", "confId%d.pdb" % cid, "-fi", "pdb", "-o", "confId%d.gjf" % cid, "-fo", "gcrt", "-nc", "%d" % self.formal_charge, "-gn", "%%NProcShared=%d" % self.nproc, "-gv", "1", "-ge", "confId%d.gesp" % cid, "-s", "%d" % self.antechamber_status, "-pf", self.raif])
                    #print("Running Gaussian 09 to calcualte ESP charge.")
                    #subprocess.run(["g09","confId%d.gjf" % cid])
                else:
                    raise RuntimeError("Gaussian is not available. Fail to calculate ESP charge.")
                subprocess.run(["antechamber", "-i", "confId%d.log" % cid, "-fi", "gout", "-o", "confId%d.mol2" % cid, "-fo", "mol2", "-c", "resp", "-nc", "%d" % self.formal_charge, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
            else:
                raise RuntimeError("Invalid charge method: %s" % self.charge_method)
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
                    std = np.std(v)
                    new_sym_chg[k] = [avg, std]
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


class residue:
    def __init__(self, charge_method='bcc', nproc=1, antechamber_status=1, remove_antechamber_intermediate_files=True):
        '''
        Create Amber residue
        '''
        self.charge_method = charge_method
        self.nproc = nproc
        self.antechamber_status = antechamber_status
        if remove_antechamber_intermediate_files:
            self.raif = 'y'
        else:
            self.raif = 'n'
        self.restype = None
        self.resname_head = None    # residue name for head capping residue
        self.resname_tail = None    # residue name for tail capping residue
        self.smiles_with_dummy = None
        self.smiles = None
        self.num_dummy = None
        self.replace_dummy_with = []
        self.work_path = None
        self.database = None
        self.num_confs = 5

    def set_restype(self, restype):
        '''
        Set residue type. Choose from "tail", "core", "bridge"
        '''
        if restype in ["tail", "core", "bridge"]:
            self.restype = restype
        else:
            raise RuntimeError("Residue type should be chosen from [\"tail\", \"core\", \"bridge\"]")

        if restype == "tail":
            self.resname_head = "TLA"
            self.resname_tail = "TLB"
        elif restype == "core":
            self.resname_head = "CRA"
            self.resname_tail = "CRB"
        elif restype == "bridge":
            ### bridge should be a symmetric molecule, so there will be only one residue name
            self.resname_head = "BRD"
            self.resname_tail = "BRD"


    def set_smiles(self, smiles):
        smiles = AllChem.MolToSmiles(AllChem.MolFromSmiles(smiles))    # get canonical smiles
        self.num_dummy = smiles.count("*")
        if self.num_dummy == 0:
            raise RuntimeError("Need dummy atoms to perform reaction.")
        elif self.num_dummy == 1:
            if smiles.find("[1*]") == -1:
                raise RuntimeError("Dummy atom should be label as \"[1*]\".")
        elif self.num_dummy == 2:
            if smiles.find("[1*]") == -1 or smiles.find("[2*]") == -1:
                raise RuntimeError("Dummy atom should be label as \"[1*]\" (as head dummy atom) and \"[2*]\" (as tail dummy atom).")
        else:
            raise RuntimeError("Too many dummy atoms. Should be less than three.")
        self.smiles_with_dummy = smiles
        

    def set_dummy_replacement(self, dummy_replacement):
        if self.smiles_with_dummy is None:
            raise RuntimeError("Please set dummy replacement first using set_dummy_replacement().")
        if isinstance(dummy_replacement, list) and len(dummy_replacement) == self.num_dummy:
            for dummy in dummy_replacement:
                if len(dummy) == 1 or dummy == "[NH3+]":
                    self.replace_dummy_with.append(dummy)
                else:
                    raise RuntimeError("Invalid dummy replacement atom.")
        else:
            raise RuntimeError("Invalid dummy replacement atom. It should be a list a atom character, with the same number of input dummy atoms.")

        if self.num_dummy == 1:
            self.smiles = self.smiles_with_dummy
            loc = self.smiles.find("[1*]")
            self.smiles = self.smiles[:loc] + self.replace_dummy_with[0] + self.smiles[loc+4:]
        elif self.num_dummy == 2:
            self.smiles = self.smiles_with_dummy
            ### replace the first dummy atom
            loc = self.smiles.find("[1*]")
            self.smiles = self.smiles[:loc] + self.replace_dummy_with[0] + self.smiles[loc+4:]
            ### replace the second dummy atom
            loc = self.smiles.find("[2*]")
            self.smiles = self.smiles[:loc] + self.replace_dummy_with[1] + self.smiles[loc+4:]


    def set_work_path(self, work_path):
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def set_database(self, database=None):
        '''
        default path is the silc/data/residue
        '''
        if database is None:
            self.database = str(files('silc.data').joinpath('residue'))
        elif database[0] == '/' or database == '~':
            self.database = database
        else:
            self.database = os.path.join(os.path.abspath(os.getcwd()), database)


    def set_num_confs(self, num_confs):
        '''
        Number of conformers for partial charge calculation.
        '''
        self.num_confs = num_confs


    def run(self):
        # check if the residue already exists in the work path
        # if so, then don't need to generate the residue again
        if self._check_path(self.work_path):
            print("Residue exists in work path at %s" % self.work_path)
            return
        
        # check if the residue already exists in the database
        # if so, copy the residue from the database
        res_datapath = self._check_database()
        if res_datapath:
            shutil.copytree(res_datapath, self.work_path)
            return

        # if the residue is not in the current work path, nor in the database
        # then start to prepare the residue
        cwd = os.getcwd()
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)
        os.makedirs(self.work_path)
        os.chdir(self.work_path)
        
        with open("smiles.txt", "w") as f:
            f.write("%s\n" % self.smiles_with_dummy)
            for i in range(self.num_dummy):
                f.write("[%d*]:%s\n" % (i+1, self.replace_dummy_with[i]))

        img = Draw.MolsToGridImage([AllChem.MolFromSmiles(self.smiles_with_dummy),
                                    AllChem.MolFromSmiles(self.smiles)], 
                                    molsPerRow=1, subImgSize=(800, 400), useSVG=True)
        with open("residue.svg", "w") as svg:
            if isinstance(img, str):
                svg.write(img)
            else:
                svg.write(img.data)

        # calculate partial charge
        chg = charge(charge_method=self.charge_method, nproc=self.nproc)
        chg.set_work_path(os.path.join(self.work_path, "amber_charge"))
        chg.set_molecule_from_smile(self.smiles)
        chg.set_num_confs(self.num_confs)
        chg.run()
        chf = chg.write_charge(file_name=os.path.join("amber_charge", "partial_charge.txt"))
        
        # identify dummy atoms
        mol = AllChem.MolFromSmiles(self.smiles_with_dummy)
        dummy = []
        for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetMass() == 0:
                dummy.append(i)
        if len(dummy) != self.num_dummy:
            raise RuntimeError("Something is wrong")
        # identify head ([1*]) and tail ([2*]) dummy atoms
        # atom index follows the sequence in the input smiles string, and this can be use to identify which dummy atom is head or tail
        # dummy[0] is the head dummy atom, and dummy[1] is the tail dummy atom
        if len(dummy) == 2:
            loc1 = self.smiles_with_dummy.find("[1*]")
            loc2 = self.smiles_with_dummy.find("[2*]")
            if loc1 > loc2:
                dummy_tmp = dummy[0]
                dummy[0] = dummy[1]
                dummy[1] = dummy_tmp

        mol = AllChem.MolFromPDBBlock(AllChem.MolToPDBBlock(chg.mol), removeHs=False)
        remove = []
        head = []   # head[0] is the head atom, and head[1] is the tail atom
        for i in range(len(dummy)):
            atom = mol.GetAtomWithIdx(dummy[i])
            remove.append(atom.GetMonomerInfo().GetName().strip())
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                neigh_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neigh_atom.GetAtomicNum() == 1:
                    remove.append(neigh_atom.GetMonomerInfo().GetName().strip())
                elif neigh_atom.GetAtomicNum() in [6,7,8]:
                    head.append(neigh_atom.GetMonomerInfo().GetName().strip())

        fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.smiles))    # formal charge

        # prepare mainchain file
        if len(head) == 2:
            # head residue
            with open("mainchain_head.mc", "w") as f:
                f.write("HEAD_NAME %s\n" % head[0])
                f.write("TAIL_NAME %s\n" % head[1])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
            # tail residue (reverse the head and tail atoms)
            with open("mainchain_tail.mc", "w") as f:
                f.write("HEAD_NAME %s\n" % head[1])
                f.write("TAIL_NAME %s\n" % head[0])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
        elif len(head) == 1:
            # with only one dummy atom, we need to find another chain end
            chain_end_atom_idx = []
            mol = AllChem.MolFromPDBBlock(AllChem.MolToPDBBlock(chg.mol))
            for atom in mol.GetAtoms():
                ## Check if the atom is connected to only one other atom (chain end)
                if atom.GetDegree() == 1:
                    chain_end_atom_idx.append(atom.GetIdx())
            ## calculate the topological distance matrix, to get the number of bonds between chain ends
            dm = AllChem.GetDistanceMatrix(mol)
            ## for all chain ends, find the one farthest from the dummy atom
            distance = 0
            end = dummy[0]
            for atomi in chain_end_atom_idx:
                if dm[atomi,dummy[0]] > distance:
                    distance = dm[atomi,dummy[0]]
                    end = atomi
            end = mol.GetAtomWithIdx(end).GetMonomerInfo().GetName().strip()

            # head capping residue
            with open("mainchain_head.mc", "w") as f:
                f.write("HEAD_NAME %s\n" % end)
                f.write("TAIL_NAME %s\n" % head[0])
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
            # tail capping residue
            with open("mainchain_tail.mc", "w") as f:
                f.write("HEAD_NAME %s\n" % head[0])
                f.write("TAIL_NAME %s\n" % end)
                for i in range(len(remove)):
                    f.write("OMIT_NAME %s\n" % remove[i])
                f.write("CHARGE %g\n" % fc)
        else:
            raise RuntimeError("cannot create a residue for the molecule given its SMILES: %s" % self.smiles_with_dummy)

        # prepare residues with amber tools
        subprocess.run(["antechamber", "-i", os.path.join("amber_charge", "confId0.mol2"), "-fi", "mol2", "-o", "molecule.ac", "-fo", "ac", "-nc", "%d" % fc, "-c", "rc", "-cf", chf, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
        subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_head.prepi", "-f", "prepi", "-m", "mainchain_head.mc", "-rn", self.resname_head, "-rf", "molecule_head.res"])
        subprocess.run(["prepgen", "-i", "molecule.ac", "-o", "molecule_tail.prepi", "-f", "prepi", "-m", "mainchain_tail.mc", "-rn", self.resname_tail, "-rf", "molecule_tail.res"])
        os.remove("ATOMTYPE.INF")
        os.remove("NEWPDB.PDB")
        os.remove("PREP.INF")
        os.chdir(cwd)


    def _check_path(self, path):
        try:
            with open(os.path.join(path, "smiles.txt"), "r") as f:
                # check smiles string (with dummy atoms)
                s = f.readline().strip('\n')
                if s != self.smiles_with_dummy:
                    return False
                # check dummy atoms
                dummy_flag = True
                for i in range(self.num_dummy):
                    dummy = f.readline().strip('\n')[5:]
                    if dummy != self.replace_dummy_with[i]:
                        dummy_flag = False
                if not dummy_flag:
                    return False
            # check if there are enough conformer
            for i in range(self.num_confs):
                with open(os.path.join(path, "amber_charge", "confId%i.mol2" % i), "r") as f:
                    lines = f.readlines()
                    if len(lines) == 0:
                        return False
            # check if prepi file exists:
            if not os.path.isfile(os.path.join(path,"molecule_head.prepi")) or not os.path.isfile(os.path.join(path,"molecule_tail.prepi")):
                return False
            return True  
        except:
            return False
        

    def _check_database(self):
        '''
        Check if the amber residue is already created and included in the database path

        return None is not found, else return the path of the residue from the database
        '''
        res_path = os.path.join(self.database, self.restype)
        if not os.path.exists(res_path):
            reture None
        sub_path = [f.path for f in os.scandir(res_path) if f.is_dir()]
        for sp in sub_path:
            if self._check_path(sp):
                print("Find residue in database at %s" % sp)
                return sp
        return None
