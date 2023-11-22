import os
import shutil
import subprocess
from importlib_resources import files

import numpy as np
import matplotlib.pyplot as plt
from rdkit.Chem import AllChem

from . import util
from .force_field import gaff2


class binding_molecule:
    '''
    Build a ditopic motif molecule
    '''
    def __init__(self, charge_method='bcc'):
        self.charge_method = charge_method
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
        self.resource_path = None


    def set_core_smiles(self, smiles, dummy_replacement):
        if smiles.count("*") != 2:
            raise RuntimeError("Need 2 dummy atoms for core SMILES to perform reaction.")
        self.core_smiles = smiles
        self.core_dummy_replacement = dummy_replacement


    def set_tail_smiles(self, smiles, dummy_replacement):
        if smiles.count("*") != 1:
            raise RuntimeError("Need 1 dummy atom for tail SMILES to perform reaction.")
        self.tail_smiles = smiles
        self.tail_dummy_replacement = dummy_replacement


    def set_bridge_smiles(self, smiles, dummy_replacement):
        if smiles.count("*") != 2:
            raise RuntimeError("Need 2 dummy atoms for bridge SMILES to perform reaction.")
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


    def set_resource_path(self, resource_path=None):
        if not resource_path:
            self.resource_path = files('silc.data.tutorial').joinpath("amber_charge")
        else:
            self.resource_path = resource_path


    def clear_work_path(self):
        if os.path.exists(self.work_path):
            shutil.rmtree(self.work_path)


    def create_ditopic_molecule(self, solvate=False):
        '''
        Create a ditopic molecule with the following structure: tail-core-bridge-core-tail
        '''
        cwd = os.getcwd()

        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        # check if residues are already prepared, if not, prepare them.
        if not self.check_amber_residues("COR", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge):
            self.prepare_amber_residue("COR", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge)
        if not self.check_amber_residues(["TLA", "TLB"], self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge):
            self.prepare_amber_residue(["TLA", "TLB"], self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge)
        if not self.check_amber_residues("BRD", self.bridge_smiles, self.bridge_dummy_replacement, self.bridge_num_confs_for_charge):
            self.prepare_amber_residue("BRD", self.bridge_smiles, self.bridge_dummy_replacement, self.bridge_num_confs_for_charge)

        # combine residues
        with open("tleap_ditopic.in", "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("loadamberprep COR/molecule.prepi\n")
            f.write("loadamberprep TLA_TLB/molecule_head.prepi\n")
            f.write("loadamberprep TLA_TLB/molecule_tail.prepi\n")
            f.write("loadamberprep BRD/molecule.prepi\n")
            f.write("mol = sequence {TLA COR BRD COR TLB}\n")
            f.write("savepdb mol ditopic.pdb\n")
            f.write("savemol2 mol ditopic.mol2 1\n")
            f.write("saveamberparm mol ditopic.prmtop ditopic.rst7\n")
            if solvate:
                f.write("source leaprc.water.tip3p\n")
                f.write("solvateOct mol TIP3PBOX 14.0\n")
                f.write("savepdb mol ditopic_solv.pdb\n")
                f.write("saveamberparm mol ditopic_solv.prmtop ditopic_solv.rst7\n")
            f.write("quit\n")
        subprocess.run(["tleap", "-f", "tleap_ditopic.in"])
        os.rename("leap.log", "leap_ditopic.log")

        os.chdir(cwd)


    def create_binding_motif(self, solvate=False):
        '''
        Create a binding motif with the following structure: tail-core-tail
        '''
        cwd = os.getcwd()

        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        # check if residues are already prepared, if not, prepare them.
        if not self.check_amber_residues("COR", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge):
            self.prepare_amber_residue("COR", self.core_smiles, self.core_dummy_replacement, self.core_num_confs_for_charge)
        if not self.check_amber_residues(["TLA", "TLB"], self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge):
            self.prepare_amber_residue(["TLA", "TLB"], self.tail_smiles, self.tail_dummy_replacement, self.tail_num_confs_for_charge)

        # combine residues
        with open("tleap_motif.in", "w") as f:
            f.write("source leaprc.gaff2\n")
            f.write("loadamberprep COR/molecule.prepi\n")
            f.write("loadamberprep TLA_TLB/molecule_head.prepi\n")
            f.write("loadamberprep TLA_TLB/molecule_tail.prepi\n")
            f.write("mol = sequence {TLA COR TLB}\n")
            f.write("savepdb mol motif.pdb\n")
            f.write("savemol2 mol motif.mol2 1\n")
            f.write("saveamberparm mol motif.prmtop motif.rst7\n")
            if solvate:
                f.write("source leaprc.water.tip3p\n")
                f.write("solvateOct mol TIP3PBOX 14.0\n")
                f.write("savepdb mol motif_solv.pdb\n")
                f.write("saveamberparm mol motif_solv.prmtop motif_solv.rst7\n")
            f.write("quit\n")
        subprocess.run(["tleap", "-f", "tleap_motif.in"])
        os.rename("leap.log", "leap_motif.log")

        os.chdir(cwd)


    def check_amber_residues(self, resname, smiles, dummy_replacement, num_confs_for_charge):
        '''
        Check if the amber residue is created.
        '''
        if smiles.count("*") == 2 and isinstance(resname, str):
            res_path = resname
        elif smiles.count("*") == 1 and isinstance(resname, list) and len(resname) == 2:
            res_path = "%s_%s" % (resname[0], resname[1])
        else:
            raise RuntimeError("invalid input parameters")
        res_path = os.path.join(self.work_path, res_path)

        try:
            with open(os.path.join(res_path, "smiles.txt"), "r") as f:
                s = f.readline().strip('\n')
                if s != smiles:
                    return False
                dummy = f.readline().strip('\n')[2]
                if dummy != dummy_replacement:
                    return False
            if smiles.count("*") == 2:
                with open(os.path.join(res_path, "molecule.prepi"), "r") as f:
                    for i, line in enumerate(f):
                        if i == 4:
                            if line.split()[0] != resname:
                                return False
            elif smiles.count("*") == 1:
                with open(os.path.join(res_path, "molecule_head.prepi"), "r") as f:
                    for i, line in enumerate(f):
                        if i == 4:
                            if line.split()[0] != resname[0]:
                                return False
                with open(os.path.join(res_path, "molecule_tail.prepi"), "r") as f:
                    for i, line in enumerate(f):
                        if i == 4:
                            if line.split()[0] != resname[1]:
                                return False
            for i in range(num_confs_for_charge):
                with open(os.path.join(res_path, "amber_charge", "confId%i.mol2" % i), "r") as f:
                    lines = f.readlines()
                    if len(lines) == 0:
                        return False
        except:
            return False
        return True


    def prepare_all_residues(self):
        self.prepare_amber_residue("COR", self.core_smiles, self.core_num_confs_for_charge)
        self.prepare_amber_residue(["TLA", "TLB"], self.tail_smiles, self.tail_num_confs_for_charge)
        self.prepare_amber_residue("BRD", self.bridge_smiles, self.bridge_num_confs_for_charge)


    def prepare_amber_residue(self, resname, smiles, dummy_replacement, num_confs, only_head=True):
        '''
        resname: three letter residue name.
                 If the residue contains only one dummy atom (which can "react" with dummy atoms in other residues), then
                 it should be a list of two residue name, for the residue to appear at the head and tail of a chain.
        smiles: SMILES string, with one or two dummy atom (*)
        dummy_replacement: when generate the residue, replace all dummy atoms with an atom given by dummy_replacement
                           this replacement atom should be chosen careful to reflect the chemical environment after new bond formation
        num_confs: number of conformers for charge calculation
        only_head: if the residue contains only one dummy atom (which can "react" with dummy atoms in other residues), then
                   True means the dummy atom is head atom (define by HEAD_NAME in amber main chain file) and connects to the previous residue;
                   False means the dummy atom is tail atom (define by TAIL_NAME in amber main chain file) and connects to the next residue;
        '''
        res = gaff2.residue(charge_method=self.charge_method)
        res.set_dummy_replacement(dummy_replacement)
        res.set_smiles(smiles)
        res.set_num_confs(num_confs)

        if smiles.count("*") == 2 and isinstance(resname, str):
            res.set_resname(resname)
            res.set_work_path(resname)
        elif smiles.count("*") == 1 and isinstance(resname, list) and len(resname) == 2:
            res.set_resname_head(resname[0])
            res.set_resname_tail(resname[1])
            res.set_work_path("%s_%s" % (resname[0], resname[1]))
        res.run()
