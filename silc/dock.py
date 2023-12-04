import os
from rdkit.Chem import AllChem
import meeko
from vina import Vina
import py3Dmol

from importlib_resources import files

from . import util

class dock:
    def __init__(self, sf_name, receptor_name, receptor_pdb_path=None, receptor_map_path=None):
        self.sf_name = sf_name
        self.receptor_name = receptor_name
        self.receptor_pdb_path = receptor_pdb_path
        self.receptor_map_path = receptor_map_path
        self.use_receptor_database = False
        self.receptor_pdb = None
        self.receptor = None
        self.ligand = None
        self.work_path = None
        self.result_pdbqt = {} # result pdbqt of n_ligand
        self.result_file = {} # result file of n_ligand

        self.v = Vina(sf_name=self.sf_name)

        # try to find the receptor in the database (i.e. from data/receptor)
        if self._find_receptor(self.receptor_name) and self.receptor_pdb_path is None:
            self.receptor_pdb = str(files("silc.data.receptor").joinpath("%s.pdb" % self.receptor_name))
            print("Found receptor %s in database at %s" % (self.receptor_name, self.receptor_pdb))
            self.use_receptor_database = True
        elif self._find_receptor(self.receptor_name) and self.receptor_pdb_path is not None:
            self.receptor_pdb = os.path.join(self.receptor_pdb_path, "%s.pdb" % self.receptor_name)
        else:
            raise RuntimeError("Cannot find receptor.")
        self.receptor = AllChem.MolFromPDBFile(self.receptor_pdb, removeHs=False)

        ### load affinity precalculated affinity map
        if self.use_receptor_database:
            self.v.load_maps(str(files("silc.data.receptor").joinpath("%s_%s" % (self.receptor_name, self.sf_name))))
        else:
            self.v.load_maps(os.path.join(self.receptor_map_path, "%s_%s" % (self.receptor_name, self.sf_name)))


    def _find_receptor(self, receptor_name):
        return files("silc.data.receptor").joinpath("%s.pdb" % self.receptor_name).is_file()


    def set_work_path(self, work_path="docking"):
        work_path = str(work_path)
        if work_path[0] == '/' or work_path == '~':
            self.work_path = work_path
        else:
            self.work_path = os.path.join(os.path.abspath(os.getcwd()), work_path)


    def prepare_ligand_from_smiles(self, smiles):
        '''
        Use meeko to prepare ligand
        '''
        ligand = util.gen_mol_from_smiles(smiles, hydrogen=True)
        meeko_prep = meeko.MoleculePreparation()
        meeko_prep.prepare(ligand)
        self.ligand_pdbqt = meeko_prep.write_pdbqt_string()


    def run(self, n_ligand=1, n_poses=20, exhaustiveness=32, save_name=None):
        '''
        run docking
        '''
        cwd = os.getcwd()
        if not os.path.exists(self.work_path):
            os.makedirs(self.work_path)
        os.chdir(self.work_path)

        self.v.set_ligand_from_string([self.ligand_pdbqt] * n_ligand)
        self.v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        if save_name is None:
            save_name = "%s_%d_ligand_%s" % (self.receptor_name, n_ligand, self.sf_name)
        self.v.write_poses("%s.pdbqt" % save_name, n_poses=n_poses, overwrite=True)
        self.result_file[n_ligand] = "%s.pdbqt" % save_name
        self.result_pdbqt[n_ligand] = self.v.poses(n_poses=n_poses)

        os.chdir(cwd)


    def visualize(self, n_ligand, pose_id):
        '''
        Visualize the complex structure calculated/loaded.

        n_ligand: number of ligands (the same ligand)
        pose_id:  the id of pose generated by vina

        Return nothing.
        '''
        receptor_block = AllChem.MolToMolBlock(self.receptor)
        viewer = py3Dmol.view(width=500, height=500)
        viewer.addModel(receptor_block, 'mol')
        viewer.setStyle({'model':-1,}, {'stick': {'color':'red'}})

        pdbqt_mol = meeko.PDBQTMolecule(self.result_pdbqt[n_ligand])
        mol = meeko.RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if n_ligand == 1:
            mol = [mol]
        for i in range(len(mol)):
            ligand = AllChem.Mol(mol[i], confId=pose_id)
            ligand_block = AllChem.MolToMolBlock(ligand)
            viewer.addModel(ligand_block, 'mol')
            viewer.setStyle({'model':-1,}, {'stick': {'color':'green'}})

        viewer.zoomTo()
        viewer.show()


    def ligand_mol(self, n_ligand, pose_id, hydrogen=False):
        '''
        Return a list of ligand molecules in the complex, with their conformer determined by pose_id.
        '''
        pdbqt_mol = meeko.PDBQTMolecule(self.result_pdbqt[n_ligand])
        mol = meeko.RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        ligand = []
        for i in range(len(mol)):
            moli = AllChem.Mol(mol[i], confId=pose_id)
            if not hydrogen:
                moli = AllChem.RemoveHs(moli)
            ligand.append(moli)
        return ligand


    def load_results(self, pdbqt_name, n_ligand):
        '''
        Load the vina result from PDBQT file.
        '''
        f = open(pdbqt_name, "r")
        self.result_pdbqt[n_ligand] = f.read()
        f.close()
