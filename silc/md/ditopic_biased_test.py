#!/usr/bin/env python
#from __future__ import division, print_function

import os
import dill as pickle
import numpy as np

# ParmEd Imports
from parmed import unit as u

# pysages imports
import pysages
from pysages import Grid
from pysages.colvars import Distance
from pysages.methods import ABF, CVRestraints

# silc imports
from silc.md.util import generate_simulation, prepare_restart
from silc.md.logger import ABFLogger
from silc.md.collective_variables import DistancesSum2,DistancesSum3


restart = True
timesteps = 10000000
restart_name = "state.pickle"
input_files = ['../ditopic_solv.prmtop', '../ditopic_solv.rst7']
output_name = "ditopic_solv"
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
    bridge_residue_name = "BRD"
    core_residue_name = ['CRA', 'CRB']
    residue_indexes = []
    for residue in sim.topology.residues():
        if residue.name == bridge_residue_name:
            residue_indexes.append(residue.index)
    for residue in sim.topology.residues():
        if residue.name in core_residue_name:
            residue_indexes.append(residue.index)

    atom_indexes=[[],[],[]]
    for atom in sim.topology.atoms():
        if atom.residue.index in residue_indexes:
            i = atom.residue.index
            atom_indexes[residue_indexes.index(i)].append(atom.index)
    print(atom_indexes)

    # Define collective variables and sampling methods.
    indices_0 = atom_indexes[0]
    indices_1 = atom_indexes[1]
    indices_2 = atom_indexes[2]
    # CV1: distance between two cores;
    # CV2: sum of distances between cores and the bridge
    # Note that CV2 is always larger than CV1, so the lower limit of CV2 should be larger than the upper limit of CV1,
    # otherwise part of the 2D sample space will be unphysical and could leads to 'nan' in the final forces and free energy.
    cv = [DistancesSum2([indices_0, indices_1, indices_2]), DistancesSum3([indices_0, indices_1, indices_2])]
    grid = Grid(lower=(1.32, 0.01), upper=(4., 2.), shape=(32, 32))
    cv_restraints = CVRestraints(lower=(1.32, 0.01), upper=(4., 2.), ku=100, kl=100)
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
