#!/usr/bin/env python
'''
Run a biased simulation that calculate the stacking free energy between two motif molecules
It applies the ABF method, using the distance between two motif cores as the collective variable.
'''


import os
import dill as pickle
import numpy as np

# ParmEd Imports
from parmed import unit as u

# pysages imports
import pysages
from pysages import Grid
from pysages.methods import ABF, CVRestraints

# silc imports
from silc.md.util import generate_simulation, prepare_restart
from silc.md.logger import ABFLogger
from silc.md.collective_variables import DistancePBC, DistancesProduct


restart = True
timesteps = 20000000
restart_name = "state.pickle"
input_files = ['../molecules/motif_2mol_solv.prmtop', '../molecules/motif_2mol_solv.rst7']
output_name = "motif_2mol_solv"
T = 298 * u.kelvin
NPT_steps = 100000
minimize_steps = 500


# find restart folder and restart files
restart_file, output_path = prepare_restart(restart_name)

if restart and restart_file is not None:
    input_files[1] = "final_NPT.xml"    # use the box size from unbiased NPT simulation.
    with restart_file.open("rb") as rf:
        state = pickle.load(rf)
    state = pysages.run(state, generate_simulation, timesteps,
                        context_args={"input_files":input_files, "output_name":output_name, "output_path":output_path,
                                      "T":T, "minimize_steps":0, "NPT_steps":0, "NVT_steps":0})
else:
    # Extract atom indexes from residue names
    sim = generate_simulation(input_files, minimize_steps=0, NPT_steps=0, NVT_steps=0)
    residue_indexes = []
    for residue in sim.topology.residues():
        if residue.name == 'BRD':
            residue_indexes.append(residue.index)
    for residue in sim.topology.residues():
        if residue.name == 'TLA':
            residue_indexes.append(residue.index)   
    for residue in sim.topology.residuges():
        if residue.name == 'TLB':
            residue_indexes.append(residue.index)

    print(residue_indexes)

    atom_indexes=[[],[],[],[],[],[]]
    for atom in sim.topology.atoms():
        if atom.residue.index in residue_indexes:
            i = atom.residue.index
            atom_indexes[residue_indexes.index(i)].append(atom.index)
    print(atom_indexes)

    box_vectors = sim.context.getSystem().getDefaultPeriodicBoxVectors()
    a = box_vectors[0].value_in_unit(u.nanometer)
    b = box_vectors[1].value_in_unit(u.nanometer)
    c = box_vectors[2].value_in_unit(u.nanometer)
    box = [a[0], b[1], c[2]]

    # Define collective variables and sampling methods.
    indices_0 = atom_indexes[0]
    indices_1 = atom_indexes[1]
    indices_2 = atom_indexes[2]
    indices_3 = atom_indexes[3]
    indices_4 = atom_indexes[4]
    indices_5 = atom_indexes[5]

    # CV1: distance between two cores;
    cv = [DistancePBC([indices_0, indices_1], box), DistancesProduct([indices_2, indices_4, indices_3, indices_5])]
    grid = Grid(lower=(0.3,0.1), upper=(4.,9.), shape=(64,64))
    cv_restraints = CVRestraints(lower=(0.3,0.1), upper=(4.,9.), ku=100, kl=100)
    sampling_method = ABF(cv, grid, restraints=cv_restraints)
    os.mkdir("ABF_logger")
    callback = ABFLogger("logger", output_path="ABF_logger", period_hist_force=timesteps//10, period_cv=timesteps//1000)
    state = pysages.run(sampling_method, generate_simulation, timesteps, callback,
                        context_args={"input_files":input_files, "output_name":output_name, "output_path":output_path,
                                      "T":T, "minimize_steps":minimize_steps, "NPT_steps":NPT_steps, "NVT_steps":0})

with open(os.path.join(output_path, restart_name), "wb") as rf:
    # TODO: Make it easier to save different states if their running time differs
    pickle.dump(state, rf)

# Analyze the results
result = pysages.analyze(state)
energy = np.asarray(result["free_energy"])
forces = np.asarray(result["mean_force"])
grid = np.asarray(result["mesh"])
np.savetxt(os.path.join(output_path, "FES.txt"), np.hstack([grid, energy.reshape(-1, 1)]))
np.savetxt(os.path.join(output_path, "Forces.txt"), np.hstack([grid, forces.reshape(-1, grid.shape[1])]))
