import os
import shutil
import subprocess
from importlib_resources import files

import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Geometry import Point3D

from . import util
from .force_field import gaff2


class binding_molecule:
    '''
    Build a ditopic molecule or a binding motif.
    '''
    def __init__(self, charge_method='bcc', nproc=1):
        self.charge_method = charge_method
        self.nproc = nproc    # only used when charge_method is "resp"
        self.core_smiles = None
        self.tail_smiles = None
        self.bridge_smiles = None
        self.core_dummy_replacement = None
        self.tail_dummy_replacement = None
        self.bridge_dummy_replacement = None
        self.core_num_confs_for_charge = 5
        self.tail_num_confs_for_charge = 5
        self.bridge_num_confs_for_charge = 5
        self.mol = None
        self.work_path = None
        self.database_path = None
        self.ditopic_pdb = None
        self.ditopic_mol2 = None
        self.ditopic_mol2_Tripos = None
        self.ditopic_charge = None
        self.ditopic_smiles = None
        self.motif_pdb = None
        self.motif_mol2 = None
        self.motif_mol2_Tripos = None
        self.motif_charge = None
        self.motif_smiles = None


    def set_core_smiles(self, smiles, dummy_replacement):
        self._check_smiles(smiles)
        if smiles.count("*") != 2:
            raise RuntimeError("Need 2 dummy atoms for core SMILES to perform reaction.")
        if not isinstance(dummy_replacement, list) or len(dummy_replacement) != 2:
            raise RuntimeError("Need 2 dummy atoms replacements for core.")
        self._check_smiles(smiles)
        self.core_smiles = smiles
        self.core_dummy_replacement = dummy_replacement


    def set_tail_smiles(self, smiles, dummy_replacement):
        self._check_smiles(smiles)
        if smiles.count("*") != 1:
            raise RuntimeError("Need 1 dummy atom for tail SMILES to perform reaction.")
        if not isinstance(dummy_replacement, list) or len(dummy_replacement) != 1:
            raise RuntimeError("Need 1 dummy atoms replacements for tail.")
        self.tail_smiles = smiles
        self.tail_dummy_replacement = dummy_replacement


    def set_bridge_smiles(self, smiles, dummy_replacement):
        self._check_smiles(smiles)
        if smiles.count("*") != 2:
            raise RuntimeError("Need 2 dummy atoms for bridge SMILES to perform reaction.")
        if not isinstance(dummy_replacement, list) or len(dummy_replacement) != 2:
            raise RuntimeError("Need 2 dummy atoms replacements for bridge.")
        self.bridge_smiles = smiles
        self.bridge_dummy_replacement = dummy_replacement


    def set_core_num_confs_for_charge(self, num_confs):
        self.core_num_confs_for_charge = num_confs


    def set_tail_num_confs_for_charge(self, num_confs):
        self.tail_num_confs_for_charge = num_confs


    def set_bridge_num_confs_for_charge(self, num_confs):
        self.bridge_num_confs_for_charge = num_confs


    def set_work_path(self, work_path="ditopic_molecule"):
        work_path = str(work_path)
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def set_database_path(self, database_path=None):
        self.database_path = database_path


    def clear_work_path(self):
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)


    def create_ditopic_molecule(self, solvate=False, counter_anion="Cl-", counter_cation="Na+", nmol=1, translate=[0.0, 0.0, 0.0]):
        '''
        Create a ditopic molecule with the following structure: tail-core-bridge-core-tail

        nmol: create nmol copies of the binding motif
        translate: should be a list of three numbers, used to translate copied molecules to prevent overlapping
        '''
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        # create a smile string for the ditopic molecule by rdkit reaction
        A = AllChem.MolFromSmiles(util.replace_dummy(self.bridge_smiles, ["K", "I"]))
        Ba = AllChem.MolFromSmiles(util.replace_dummy(self.core_smiles, ["I", "K"]))
        Bb = AllChem.MolFromSmiles(util.replace_dummy(self.core_smiles, ["K", "I"]))
        Ca = AllChem.MolFromSmiles(util.replace_dummy(self.tail_smiles, ["K"]))
        Cb = AllChem.MolFromSmiles(util.replace_dummy(self.tail_smiles, ["I"]))
        prod = self._reaction([Ca, Ba, A, Bb, Cb])
        self.ditopic_smiles = AllChem.MolToSmiles(prod)
        img = Draw.MolsToGridImage([AllChem.MolFromSmiles(self.ditopic_smiles),],
                                    molsPerRow=1, subImgSize=(1200, 800), useSVG=True)
        with open("ditopic.svg", "w") as svg:
            if isinstance(img, str):
                svg.write(img)
            else:
                svg.write(img.data)

        # prepare residues
        self._prepare_amber_residue("core", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge)
        self._prepare_amber_residue("tail", self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge)
        self._prepare_amber_residue("bridge", self.bridge_smiles, self.bridge_dummy_replacement, self.bridge_num_confs_for_charge)

        # combine residues
        if nmol == 1:
            appendix = ""
        else:
            appendix = "_%dmol" % nmol
        with open("tleap_ditopic%s_step1.in" % appendix, "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_ditopic%s_step1.log\n" % appendix)
            #f.write("loadamberparams %s\n\n" % files('silc.data.residue').joinpath('frcmod.silc'))
            f.write("loadamberprep core/molecule_head.prepi\n")
            f.write("loadamberprep core/molecule_tail.prepi\n")
            f.write("loadamberprep tail/molecule_head.prepi\n")
            f.write("loadamberprep tail/molecule_tail.prepi\n")
            f.write("loadamberprep bridge/molecule_head.prepi\n")
            f.write("loadamberprep bridge/molecule_tail.prepi\n")
            f.write("mol = sequence {TLA CRA BRD CRB TLB}\n")
            f.write("savepdb mol ditopic.pdb\n")
            f.write("savemol2 mol ditopic_Tripos.mol2 0\n")
            f.write("savemol2 mol ditopic.mol2 1\n")
            if nmol > 1:
                for i in range(nmol):
                    f.write("mol%d = sequence {TLA CRA BRD CRB TLB}\n" % i)
                    f.write("translate mol%d {%f %f %f}\n" % (i, translate[0]*i, translate[1]*i, translate[2]*i))
                f.write("mol_comb = combine {")
                for i in range(nmol):
                    f.write("mol%d " % i)
                f.write("}\n")
                f.write("savepdb mol_comb ditopic%s.pdb\n" % appendix)
                f.write("savemol2 mol_comb ditopic_Tripos%s.mol2 0\n" % appendix)
                f.write("savemol2 mol_comb ditopic%s.mol2 1\n" % appendix)
            f.write("quit\n")
        result = subprocess.run(["tleap", "-f", "tleap_ditopic%s_step1.in" % appendix])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_ditopic%s_step1.log" % appendix)
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])
        
        result = subprocess.run(["parmchk2", "-i", "ditopic%s.mol2" % appendix, "-f", "mol2", "-o", "ditopic.frcmod", "-s", "gaff2"])
        if result.returncode != 0:
            raise RuntimeError("parmchk2 run error.")

        with open("tleap_ditopic%s_step2.in" % appendix, "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_ditopic%s_step2.log\n" % appendix)
            f.write("loadamberparams ditopic.frcmod\n\n")
            f.write("mol_comb = loadmol2 ditopic%s.mol2\n" % appendix)
            f.write("saveamberparm mol_comb ditopic%s.prmtop ditopic%s.rst7\n" % (appendix, appendix))
            if solvate:
                f.write("source leaprc.water.tip3p\n")
                f.write("solvateBox mol_comb TIP3PBOX 14.0 iso\n")
                fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.ditopic_smiles))    # formal charge
                if fc > 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 mol_comb %s 0\n" % counter_anion)
                elif fc < 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 mol_comb %s 0\n" % counter_cation)
                f.write("savepdb mol_comb ditopic%s_solv.pdb\n" % appendix)
                f.write("saveamberparm mol_comb ditopic%s_solv.prmtop ditopic%s_solv.rst7\n" % (appendix, appendix))
            f.write("quit\n")
        result = subprocess.run(["tleap", "-f", "tleap_ditopic%s_step2.in" % appendix])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_ditopic%s_step2.log" % appendix)
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])

        self.ditopic_pdb = os.path.join(os.path.abspath(os.getcwd()), "ditopic.pdb")
        self.ditopic_mol2 = os.path.join(os.path.abspath(os.getcwd()), "ditopic.mol2")
        self.ditopic_mol2_Tripos = os.path.join(os.path.abspath(os.getcwd()), "ditopic_Tripos.mol2")
        self.ditopic_charge = util.read_mol2_charge(self.ditopic_mol2)
    
        os.chdir(cwd)


    def create_binding_motif(self, solvate=False, counter_anion="Cl-", counter_cation="Na+", nmol=1, translate=[0.0, 0.0, 0.0]):
        '''
        Create a binding motif with the following structure: tail-core-tail

        nmol: create nmol copies of the binding motif
        translate: should be a list of three numbers, used to translate copied molecules to prevent overlapping
        '''
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        # create a smile string for the binding motif by rdkit reaction
        B = AllChem.MolFromSmiles(util.replace_dummy(self.core_smiles, ["K", "I"]))
        Ca = AllChem.MolFromSmiles(util.replace_dummy(self.tail_smiles, ["K"]))
        Cb = AllChem.MolFromSmiles(util.replace_dummy(self.tail_smiles, ["I"]))
        prod = self._reaction([Ca, B, Cb])
        self.motif_smiles = AllChem.MolToSmiles(prod)
        img = Draw.MolsToGridImage([AllChem.MolFromSmiles(self.motif_smiles),],
                                    molsPerRow=1, subImgSize=(1200, 800), useSVG=True)
        with open("motif.svg", "w") as svg:
            if isinstance(img, str):
                svg.write(img)
            else:
                svg.write(img.data)

        # prepare residues
        self._prepare_amber_residue("core", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge)
        self._prepare_amber_residue("tail", self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge)

        # combine residues
        if nmol == 1:
            appendix = ""
        else:
            appendix = "_%dmol" % nmol
        with open("tleap_motif%s_step1.in" % appendix, "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_motif%s_step1.log\n" % appendix)
            f.write("loadamberprep core/molecule_head.prepi\n")
            f.write("loadamberprep core/molecule_tail.prepi\n")
            f.write("loadamberprep tail/molecule_head.prepi\n")
            f.write("loadamberprep tail/molecule_tail.prepi\n")
            f.write("mol = sequence {TLA CRA TLB}\n")
            f.write("savepdb mol motif.pdb\n")
            f.write("savemol2 mol motif_Tripos.mol2 0\n")
            f.write("savemol2 mol motif.mol2 1\n")
            if nmol > 1:
                for i in range(nmol):
                    f.write("mol%d = sequence {TLA CRA TLB}\n" % i)
                    f.write("translate mol%d {%f %f %f}\n" % (i, translate[0]*i, translate[1]*i, translate[2]*i))
                f.write("mol_comb = combine {")
                for i in range(nmol):
                    f.write("mol%d " % i)
                f.write("}\n")
                f.write("savepdb mol_comb motif%s.pdb\n" % appendix)
                f.write("savemol2 mol_comb motif_Tripos%s.mol2 0\n" % appendix)
                f.write("savemol2 mol_comb motif%s.mol2 1\n" % appendix)
            f.write("quit\n")
        result = subprocess.run(["tleap", "-f", "tleap_motif%s_step1.in" % appendix])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_motif%s_step1.log" % appendix)
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])
        
        result = subprocess.run(["parmchk2", "-i", "motif.mol2", "-f", "mol2", "-o", "motif.frcmod", "-s", "gaff2"])
        if result.returncode != 0:
            raise RuntimeError("parmchk2 run error.")

        with open("tleap_motif%s_step2.in" % appendix, "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_motif%s_step2.log\n" % appendix)
            f.write("loadamberparams motif.frcmod\n\n")
            f.write("mol_comb = loadmol2 motif%s.mol2\n" % appendix)
            f.write("saveamberparm mol_comb motif%s.prmtop motif%s.rst7\n" % (appendix, appendix))
            if solvate:
                f.write("source leaprc.water.tip3p\n")
                f.write("solvateBox mol_comb TIP3PBOX 14.0 iso\n")
                fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.motif_smiles))    # formal charge
                if fc > 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 mol_comb %s 0\n" % counter_anion)
                elif fc < 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 mol_comb %s 0\n" % counter_cation)
                f.write("savepdb mol_comb motif%s_solv.pdb\n" % appendix)
                f.write("saveamberparm mol_comb motif%s_solv.prmtop motif%s_solv.rst7\n" % (appendix, appendix))
            f.write("quit\n")
        result = subprocess.run(["tleap", "-f", "tleap_motif%s_step2.in" % appendix])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_motif%s_step2.log" % appendix)
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])

        self.motif_pdb = os.path.join(os.path.abspath(os.getcwd()), "motif.pdb")
        self.motif_mol2 = os.path.join(os.path.abspath(os.getcwd()), "motif.mol2")
        self.motif_mol2_Tripos = os.path.join(os.path.abspath(os.getcwd()), "motif_Tripos.mol2")
        self.motif_charge = util.read_mol2_charge(self.motif_mol2)

        os.chdir(cwd)


    def write_ditopic_molecule_charge_file(self, file_name):
        '''
        write the partial charge into a file, using the same format by antechamber
        '''
        count = 0
        with open(file_name, "w") as f:
            for key, value in self.ditopic_charge.items():
                f.write("%10.6lf" % value)
                count += 1
                if count == 8:
                    count = 0
                    f.write("\n")
        file_abspath = os.path.join(os.path.abspath(os.getcwd()), file_name)
        return file_abspath


    def write_binding_motif_charge_file(self, file_name):
        '''
        write the partial charge into a file, using the same format by antechamber
        '''
        count = 0
        with open(file_name, "w") as f:
            for key, value in self.motif_charge.items():
                f.write("%10.6lf" % value)
                count += 1
                if count == 8:
                    count = 0
                    f.write("\n")
        file_abspath = os.path.join(os.path.abspath(os.getcwd()), file_name)
        return file_abspath


    def _prepare_amber_residue(self, restype, smiles, dummy_replacement, num_confs):
        '''
        restype: residue type, choose from ["tail", "core", "bridge"]
        smiles: SMILES string, with one or two dummy atom (*)
        dummy_replacement: when generate the residue, replace all dummy atoms with an atom given by dummy_replacement
                           this replacement atom should be chosen careful to reflect the chemical environment after new bond formation
        num_confs: number of conformers for charge calculation
        '''
        res = gaff2.residue(charge_method=self.charge_method, nproc=self.nproc)
        res.set_restype(restype)
        res.set_smiles(smiles)
        res.set_dummy_replacement(dummy_replacement)
        res.set_num_confs(num_confs)
        res.set_work_path(restype)
        res.set_database(self.database_path)
        res.run()


    def _check_smiles(self, smiles):
        num_dummy = smiles.count("*")
        if num_dummy == 0:
            raise RuntimeError("Need dummy atoms to perform reaction.")
        elif num_dummy == 1:
            if smiles.find("[1*]") == -1:
                raise RuntimeError("Dummy atom should be label as \"[1*]\".")
        elif num_dummy == 2:
            if smiles.find("[1*]") == -1 or smiles.find("[2*]") == -1:
                raise RuntimeError("Dummy atom should be label as \"[1*]\" (as head dummy atom) and \"[2*]\" (as tail dummy atom).")
        else:
            raise RuntimeError("Too many dummy atoms. Should be less than three.")


    def _reaction(self, reactants):
        '''
        Carry out the following reaction and return the final product (Z):
        reactants[0] + reactants[1] -> A + reactants[2] -> B + ... -> Y + reactants[-1] -> Z
        
        reactants: A list of rdkit molecule, each should have a [1K] and [2I] atom representing the head and tail dummy atoms.
                   For the first reactant, it should have a [2I] atom as the tail dummy atom.
                   For the last reactant, it should have a [1K] atom as the head dummy atom.
        '''
        ### define reaction
        ### K is used as the tail dummy atom in the head molecule, I is used as the head dummy atom in the tail molecule
        rxn = rdChemReactions.ReactionFromSmarts("[C,c,O,N:1][K:2].[I:3][C,c,O,N:4]>>[C,c,O,N:1][C,c,O,N:4]")
        prod = reactants[0]
        if len(reactants) >= 2:
            for i in range(len(reactants)-1):
                reacts = (prod, reactants[i+1])
                products = rxn.RunReactants(reacts)
                prod = products[0][0]
                AllChem.SanitizeMol(prod)
        else:
            raise RuntimeError("Not enough reactants.")
        return prod


class complex():
    def __init__(self, antechamber_status=1, remove_antechamber_intermediate_files=True):
        self.antechamber_status = antechamber_status
        if remove_antechamber_intermediate_files:
            self.raif = 'y'
        else:
            self.raif = 'n'
        self.work_path = None
        self.binding_molecule = None
        self.dock = None


    def set_work_path(self, work_path="complex"):
        work_path = str(work_path)
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def set_binding_molecule(self, binding_molecule):
        self.binding_molecule = binding_molecule


    def set_dock(self, dock):
        self.dock = dock
        

    def clear_work_path(self):
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)


    def create_receptor_motif_complex(self, n_motif, dock_pose_id, solvate=False, counter_anion="Cl-", counter_cation="Na+"):
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        ligands = self.dock.ligand_mol(n_ligand=n_motif, pose_id=dock_pose_id)
        core = AllChem.MolFromSmiles(util.replace_dummy(self.binding_molecule.core_smiles, new=["[H]", "[H]"], replace_mass_label=True))
        motifs = []
        expanded_cores = []
        for i in range(n_motif):
            expanded_corei = util.expand_substructure(ligands[i], core, expand_iteration=1)
            expanded_cores.append(expanded_corei)
            motif_template = AllChem.MolFromSmiles(self.binding_molecule.motif_smiles)
            motif_template = AllChem.AddHs(motif_template)
            motif = AllChem.MolFromPDBFile(self.binding_molecule.motif_pdb, removeHs=False)
            motif = AllChem.Mol(motif)
            motif = AllChem.AssignBondOrdersFromTemplate(motif_template, motif)
            motif = util.align_to_substructure(motif, expanded_corei, stretch=True)
            motifs.append(motif)
        self.dock.receptor, motifs = util.optimize_complex(self.dock.receptor, motifs, expanded_cores)

        util.save_pdb(self.dock.receptor, "receptor.pdb", bond=False)
        for i in range(n_motif):
            util.save_pdb(motifs[i], "motif%d_nobond.pdb" % i, bond=False)
            self.binding_molecule.write_binding_motif_charge_file("motif%d_charge.txt" % i)
            fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.binding_molecule.motif_smiles))    # formal charge
            #subprocess.run(["antechamber", "-i", "motif%d_nobond.pdb" % i, "-fi", "pdb", "-o", "motif%d.mol2" % i, "-fo", "mol2", "-nc", "%d" % fc, "-c", "rc", "-cf", "motif%d_charge.txt" % i, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
            #subprocess.run(["parmchk2", "-i", "motif%d.mol2" % i, "-f", "mol2", "-o", "motif%d.frcmod" % i, "-s", "gaff2"])

        with open("tleap_motif_complex.in", "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_motif_complex.log\n")
            f.write("source %s\n\n" % files('silc.data.receptor.q4md-CD').joinpath('script1.ff'))
            f.write("loadamberparams %s/motif.frcmod\n\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            for i in range(n_motif):
                f.write("motif%d = sequence {TLA CRA TLB}\n" % i)
                f.write("motif%d = loadpdb motif%d_nobond.pdb\n" % (i, i))
            f.write("\n")

            f.write("receptor = loadPDB receptor.pdb\n")
            f.write("bond receptor.1.C4 receptor.8.O1\n\n")

            f.write("complex = combine {receptor")
            for i in range(n_motif):
                f.write(" motif%d" % i)
            f.write("}\n")
            f.write("charge complex\n")
            f.write("check complex\n\n")

            f.write("saveamberparm complex motif_complex.prmtop motif_complex.rst7\n")
            f.write("savemol2 complex motif_complex_Tripos.mol2 0\n")
            f.write("savemol2 complex motif_complex.mol2 1\n")
            f.write("savepdb complex motif_complex.pdb\n\n")

            if solvate:
                f.write("loadoff solvents.lib\n")
                f.write("solvateBox complex TIP3PBOX 14.0 iso\n")
                if fc > 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_anion)
                elif fc < 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_cation)
                f.write("saveamberparm complex motif_complex_solv.prmtop motif_complex_solv.rst7\n")
                f.write("savepdb complex motif_complex_solv.pdb\n\n")

            f.write("quit\n")
        result = subprocess.run(["tleap", "-I", "%s" % files('silc.data.receptor').joinpath('q4md-CD'), "-f", "tleap_motif_complex.in"])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_motif_complex.log")
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])

        os.chdir(cwd)


    def create_receptor_motif_2to2_complex(self, dock_pose_id, receptor1_translation=[0.,0.,0.], receptor2_translation=[0.,0.,0.], solvate=False, counter_anion="Cl-", counter_cation="Na+"):
        '''
        receptor1_translation & receptor2_translation:
            translate receptors along their three principal axes (obtained from the eigendecomposition of inertia matrix)
            receptor1: the receptor that already bind to the ligands in dock
            receptor2: the extra receptor that is added to the complexes
        '''
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        ligands = self.dock.ligand_mol(n_ligand=2, pose_id=dock_pose_id)
        core = AllChem.MolFromSmiles(util.replace_dummy(self.binding_molecule.core_smiles, new=["[H]", "[H]"], replace_mass_label=True))
        motifs = []
        expanded_cores = []
        for i in range(2):
            expanded_corei = util.expand_substructure(ligands[i], core, expand_iteration=1)
            expanded_cores.append(expanded_corei)
            motif_template = AllChem.MolFromSmiles(self.binding_molecule.motif_smiles)
            motif_template = AllChem.AddHs(motif_template)
            motif = AllChem.MolFromPDBFile(self.binding_molecule.motif_pdb, removeHs=False)
            motif = AllChem.Mol(motif)
            motif = AllChem.AssignBondOrdersFromTemplate(motif_template, motif)
            motif = util.align_to_substructure(motif, expanded_corei, stretch=True)
            motifs.append(motif)
        self.dock.receptor, motifs = util.optimize_complex(self.dock.receptor, motifs, expanded_cores)

        ### receptor's coordinate and principle axis of inertia
        receptor = AllChem.Mol(self.dock.receptor)
        conf = receptor.GetConformer()
        receptor_coord = []
        for i in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            receptor_coord.append((pos.x, pos.y, pos.z))
        receptor_coord = np.array(receptor_coord)
        cog = np.mean(receptor_coord, axis=0)
        receptor_coord -= cog
        ### inertia matrix
        inertia = np.zeros((3,3))
        for i in range(len(receptor_coord)):
            p = receptor_coord[i]
            inertia += np.dot(p,p)*np.identity(3)-np.outer(p,p)
        ### principle axis of inetria
        _, evec = np.linalg.eigh(inertia)

        ### translocate the 1st receptor
        receptor1 = AllChem.Mol(self.dock.receptor)
        conf = receptor1.GetConformer()
        ### assume the receptor is plate-like, and the main axis is perpendicular to the plane of the plate,
        ### this axis is evec[:,-1] with the largest eigen value (principal moments of inertia)
        delta = Point3D(*(evec[:,-1] * receptor1_translation[0] + evec[:,0] * receptor1_translation[1] + evec[:,1] * receptor1_translation[2]))
        for i in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, pos + delta)

        ### prepare the 2nd receptor
        receptor2 = AllChem.Mol(self.dock.receptor)
        conf = receptor2.GetConformer()
        ### inversion of the 2nd receptor w.r.t. center of gemoretry (cog)
        for i in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            delta = Point3D(*(-2*receptor_coord[i]))
            conf.SetAtomPosition(i, pos + delta)
        ### translocate the 2nd receptor along
        ### assume the receptor is plate-like, and the main axis is perpendicular to the plane of the plate,
        ### this axis is evec[:,-1] with the largest eigen value (principal moments of inertia)
        delta = Point3D(*(evec[:,-1] * receptor2_translation[0] + evec[:,0] * receptor2_translation[1] + evec[:,1] * receptor2_translation[2]))
        for i in range(conf.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, pos + delta)
    
        util.save_pdb(receptor1, "receptor0.pdb", bond=False)
        util.save_pdb(receptor2, "receptor1.pdb", bond=False)
        for i in range(2):
            util.save_pdb(motifs[i], "motif%d_nobond.pdb" % i, bond=False)
            self.binding_molecule.write_binding_motif_charge_file("motif%d_charge.txt" % i)
            fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.binding_molecule.motif_smiles))    # formal charge
            #subprocess.run(["antechamber", "-i", "motif%d_nobond.pdb" % i, "-fi", "pdb", "-o", "motif%d.mol2" % i, "-fo", "mol2", "-nc", "%d" % fc, "-c", "rc", "-cf", "motif%d_charge.txt" % i, "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
            #subprocess.run(["parmchk2", "-i", "motif%d.mol2" % i, "-f", "mol2", "-o", "motif%d.frcmod" % i, "-s", "gaff2"])

        with open("tleap_motif_complex.in", "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_motif_complex.log\n")
            f.write("source %s\n\n" % files('silc.data.receptor.q4md-CD').joinpath('script1.ff'))
            f.write("loadamberparams %s/motif.frcmod\n\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            for i in range(2):
                f.write("motif%d = sequence {TLA CRA TLB}\n" % i)
                f.write("motif%d = loadpdb motif%d_nobond.pdb\n" % (i, i))
            f.write("\n")
            for i in range(2):
                f.write("receptor%d = loadPDB receptor%d.pdb\n" % (i, i))
                f.write("bond receptor%d.1.C4 receptor%d.8.O1\n\n" % (i, i))

            f.write("complex = combine {receptor0 receptor1 motif0 motif1}\n")
            f.write("charge complex\n")
            f.write("check complex\n\n")

            f.write("saveamberparm complex motif_complex.prmtop motif_complex.rst7\n")
            f.write("savemol2 complex motif_complex_Tripos.mol2 0\n")
            f.write("savemol2 complex motif_complex.mol2 1\n")
            f.write("savepdb complex motif_complex.pdb\n\n")

            if solvate:
                f.write("loadoff solvents.lib\n")
                f.write("solvateBox complex TIP3PBOX 14.0 iso\n")
                if fc > 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_anion)
                elif fc < 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_cation)
                f.write("saveamberparm complex motif_complex_solv.prmtop motif_complex_solv.rst7\n")
                f.write("savepdb complex motif_complex_solv.pdb\n\n")

            f.write("quit\n")
        result = subprocess.run(["tleap", "-I", "%s" % files('silc.data.receptor').joinpath('q4md-CD'), "-f", "tleap_motif_complex.in"])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_motif_complex.log")
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])

        os.chdir(cwd)


    def create_receptor_ditopic_complex(self, solvate=False, counter_anion="Cl-", counter_cation="Na+"):
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        ligands = self.dock.ligand_mol(n_ligand=2, pose_id=0)   
        core = AllChem.MolFromSmiles(util.replace_dummy(self.binding_molecule.core_smiles, new=["", ""], replace_mass_label=True))
        expanded_core = util.expand_substructure(ligands[0], core, expand_iteration=0)
        ditopic_template = AllChem.MolFromSmiles(self.binding_molecule.ditopic_smiles)
        ditopic_template = AllChem.AddHs(ditopic_template)
        ditopic = AllChem.MolFromPDBFile(self.binding_molecule.ditopic_pdb, removeHs=False)
        ditopic = AllChem.Mol(ditopic)
        ditopic = AllChem.AssignBondOrdersFromTemplate(ditopic_template, ditopic)
        match = ditopic.GetSubstructMatches(expanded_core)
        pair1 = [(match[0][i], i) for i in range(len(match[0]))]
        pair2 = [(match[1][i], i) for i in range(len(match[1]))]
        AllChem.AlignMol(ditopic, expanded_core, atomMap=pair1)
        comb = AllChem.CombineMols(ditopic, self.dock.receptor)
        AllChem.AlignMol(comb, expanded_core, atomMap=pair2)
        comb = AllChem.CombineMols(comb, self.dock.receptor)

        frags = AllChem.GetMolFrags(comb, asMols=True)
        ditopic = frags[0]
        receptor0 = frags[1]
        receptor1 = frags[2]

        util.save_pdb(ditopic, "ditopic_nobond.pdb", bond=False)
        util.save_pdb(receptor0, "receptor0_nobond.pdb", bond=False)
        util.save_pdb(receptor1, "receptor1_nobond.pdb", bond=False)
        self.binding_molecule.write_binding_motif_charge_file("ditopic_charge.txt")
        fc = AllChem.GetFormalCharge(AllChem.MolFromSmiles(self.binding_molecule.motif_smiles))    # formal charge
        #subprocess.run(["antechamber", "-i", "ditopic_nobond.pdb", "-fi", "pdb", "-o", "ditopic.mol2", "-fo", "mol2", "-nc", "%d" % fc, "-c", "rc", "-cf", "ditopic_charge.txt", "-at", "gaff2", "-s", "%d" % self.antechamber_status, "-pf", self.raif])
        #subprocess.run(["parmchk2", "-i", "ditopic.mol2", "-f", "mol2", "-o", "ditopic.frcmod", "-s", "gaff2"])

        with open("tleap_ditopic_complex.in", "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("logfile leap_ditopic_complex.log\n")
            f.write("source %s\n\n" % files('silc.data.receptor.q4md-CD').joinpath('script1.ff'))
            f.write("loadamberparams %s/ditopic.frcmod\n\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/core/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/tail/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/bridge/molecule_head.prepi\n" % self.binding_molecule.work_path)
            f.write("loadamberprep %s/bridge/molecule_tail.prepi\n" % self.binding_molecule.work_path)
            f.write("ditopic = sequence {TLA CRA BRD CRB TLB}\n")
            f.write("ditopic = loadpdb ditopic_nobond.pdb\n")
            f.write("receptor0 = loadPDB receptor0_nobond.pdb\n")
            f.write("bond receptor0.1.C4 receptor0.8.O1\n\n")
            f.write("receptor1 = loadPDB receptor1_nobond.pdb\n")
            f.write("bond receptor1.1.C4 receptor1.8.O1\n\n")
            f.write("\n")

            f.write("complex = combine {ditopic receptor0 receptor1}\n")
            f.write("charge complex\n")
            f.write("check complex\n\n")

            f.write("saveamberparm complex ditopic_complex.prmtop ditopic_complex.rst7\n")
            f.write("savemol2 complex ditopic_complex_Tripos.mol2 0\n")
            f.write("savemol2 complex ditopic_complex.mol2 1\n")
            f.write("savepdb complex ditopic_complex.pdb\n\n")

            if solvate:
                f.write("loadoff solvents.lib\n")
                f.write("solvateBox complex TIP3PBOX 14.0 iso\n")
                if fc > 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_anion)
                elif fc < 0:
                    f.write("source leaprc.water.tip3p\n")
                    f.write("addIons2 complex %s 0\n" % counter_cation)
                f.write("saveamberparm complex ditopic_complex_solv.prmtop ditopic_complex_solv.rst7\n")
                f.write("savepdb complex ditopic_complex_solv.pdb\n\n")

            f.write("quit\n")
        result = subprocess.run(["tleap", "-I", "%s" % files('silc.data.receptor').joinpath('q4md-CD'), "-f", "tleap_ditopic_complex.in"])
        if result.returncode != 0:
            raise RuntimeError("tleap run error.")
        leap_status = util.check_leap_log("leap_ditopic_complex.log")
        if leap_status[0] != 0:
            raise RuntimeError("tleap run error. Total number of errors = %d" % leap_status[0])

        os.chdir(cwd)
