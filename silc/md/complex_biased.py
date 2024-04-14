#!/usr/bin/env python
'''
Run a biased simulation that calculate the binding free energy between one motif and the gamma-CD.
It applies the Funnel ABF method, using the distance between the motif core and the 'bottom' of gamma-CD as the collective variable.
'''

import os
import dill as pickle
import numpy as np

# ParmEd Imports
from parmed import unit as u

# pysages imports
import pysages
from pysages import Grid
from pysages.colvars import Projection_on_Axis_mobile
from pysages.methods import CVRestraints, Funnel_ABF, Funnel_Logger, get_funnel_force

# silc imports
from silc.md.util import generate_simulation, prepare_restart

win_min = 0
win_max = 1

restart = True
timesteps = 10000000
restart_name = "state.pickle"
input_files = ['../../complexes_one_motif/motif_complex_solv.prmtop', '../../complexes_one_motif/motif_complex_solv.rst7']
output_name = "complex_solv"
T = 298 * u.kelvin
NPT_steps = 100000
minimize_steps = 500

# find restart folder and restart files
restart_file, output_path = prepare_restart(restart_name)

if restart and restart_file is not None:
    input_files[1] = "final_NPT.xml"    # use the box size from unbiased NPT simulation.
    print(restart_file)
    with restart_file.open("rb") as rf:
        state = pickle.load(rf)
    state = pysages.run(state, generate_simulation, timesteps,
                        context_args={"input_files":input_files, "output_name":output_name, "output_path":output_path,
                                      "T":T, "minimize_steps":0, "NPT_steps":0, "NVT_steps":0})
else:
    sim = generate_simulation(input_files, minimize_steps=0, NPT_steps=0, NVT_steps=0)
    core_residue_names = ["CRA",]
    CD_residue_names = ["MGO",]
    CD_upper_atom_names = ['O2', 'O3']
    CD_lower_atom_names = ['C6', 'O6']

    core_atom_indexes = []
    CD_atom_indexes = []
    CD_upper_atom_indexes = []
    CD_lower_atom_indexes = []
    for atom in sim.topology.atoms():
        if atom.residue.name in core_residue_names:
            core_atom_indexes.append(atom.index)
        elif atom.residue.name in CD_residue_names:
            CD_atom_indexes.append(atom.index)
            if atom.name in CD_upper_atom_names:
                CD_upper_atom_indexes.append(atom.index)
            elif atom.name in CD_lower_atom_names:
                CD_lower_atom_indexes.append(atom.index)
    print(core_atom_indexes)
    print(CD_atom_indexes)
    print(CD_upper_atom_indexes)
    print(CD_lower_atom_indexes)

    box_vectors = sim.context.getSystem().getDefaultPeriodicBoxVectors()
    a = box_vectors[0].value_in_unit(u.nanometer)
    b = box_vectors[1].value_in_unit(u.nanometer)
    c = box_vectors[2].value_in_unit(u.nanometer)
    box = [a[0], b[1], c[2]]

    state = sim.context.getState(getPositions=True)
    pos = state.getPositions(asNumpy=True)
    ref = pos.value_in_unit(u.nanometer)[CD_atom_indexes]
    A = np.mean(pos.value_in_unit(u.nanometer)[CD_lower_atom_indexes], axis=0)
    B = np.mean(pos.value_in_unit(u.nanometer)[CD_upper_atom_indexes], axis=0)
    Z_0 = np.min(np.linalg.norm(pos.value_in_unit(u.nanometer)[CD_lower_atom_indexes] - A, axis=1))
    Z_0 -= 0.1
    print(Z_0)
    Zcc = 2.0
    R_cyl = 0.2
    k_cone = 10000.0
    k_cv = 10000.0
    cv_min = win_min
    cv_max = win_max
    cv_buffer = 0.05

    cvs = (
        Projection_on_Axis_mobile(
            [core_atom_indexes, CD_atom_indexes, [CD_atom_indexes[0]]],
            references=ref,
            weights_lig=None,
            weights_prot=None,
            A=A,
            B=B,
            box=box,
        ),
    )
    grid = Grid(lower=(cv_min, ), upper=(cv_max,), shape=(128,))
    #cv_restraints = CVRestraints(lower=(cv_min - cv_buffer, ), upper=(cv_max + cv_buffer,), ku=10000, kl=10000)

    funnel_force = get_funnel_force(
                        [core_atom_indexes, CD_atom_indexes, [CD_atom_indexes[0]]], ref, A, B,
                        Zcc, Z_0, R_cyl, k_cone, k_cv, cv_min, cv_max, cv_buffer, box
                   )

    sampling_method = Funnel_ABF(cvs, grid=grid, N=1000, ext_force=funnel_force)
    funnel_file = "funnel.dat"
    callback = Funnel_Logger(funnel_file, 1000)
    state = pysages.run(sampling_method, generate_simulation, timesteps, callback,
                            context_args={"input_files":input_files, "output_name":output_name, "output_path":output_path,
                                    "T":T, "minimize_steps":minimize_steps, "NPT_steps":NPT_steps, "NVT_steps":0})

with open(os.path.join(output_path, restart_name), "wb") as rf:
    # TODO: Make it easier to save different states if their running time differs
    pickle.dump(state, rf)


result = pysages.analyze(state)
energy = np.asarray(result["free_energy_corrected"])
forces = np.asarray(result["force_corrected"])
hist = np.asarray(result["histogram"])
grid = np.asarray(result["mesh"])
np.savetxt(os.path.join(output_path, "FES.txt"), np.hstack([grid, energy.reshape(-1, 1)]))
np.savetxt(os.path.join(output_path, "Forces.txt"), np.hstack([grid, forces.reshape(-1, grid.shape[1])]))
np.savetxt(os.path.join(output_path, "Histogram.txt"), np.hstack([grid, hist.reshape(-1, 1)]))
