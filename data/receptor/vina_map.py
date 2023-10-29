from rdkit.Chem import AllChem
import meeko
from vina import Vina

receptor_name = 'GCDOH'
receptor = AllChem.MolFromPDBFile("%s.pdb" % receptor_name)
centroid = AllChem.ComputeCentroid(receptor.GetConformer())
v = Vina(sf_name="vina")
v.set_receptor("%s.pdbqt" % receptor_name)
v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=[15, 15, 15])
v.write_maps(map_prefix_filename="%s_vina" % receptor_name)
