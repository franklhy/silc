#!/usr/bin/env python
#from __future__ import division, print_function

import sys

import numpy as np

# OpenMM Imports
#import simtk.openmm as mm
#import simtk.openmm.app as app
import openmm as mm
import openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
from parmed.openmm import StateDataReporter, NetCDFReporter

# pysages imports
import pysages
from pysages import Grid
from pysages.colvars import Distance
from pysages.methods import ABF, CVRestraints

# silc imports
from silc.md.logger import ABFLogger
from silc.md.collective_variables import DistancesSum


biased = True
timesteps = 10000000
input_files = ('../ditopic_solv.prmtop', '../ditopic_solv.rst7')
T = 298 * u.kelvin
NPT_steps = 100000
minimize_steps = 500

def generate_simulation(input_files=input_files, T=T, NPT_steps=NPT_steps, minimize_steps=minimize_steps, platform_name="CUDA"):
    # Load the Amber files
    print('Loading AMBER files...')
    complex_solv = load_file(*input_files)

    # Create the OpenMM system
    print('Creating OpenMM System')
    system = complex_solv.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=8.0*u.angstroms,
                                    constraints=app.HBonds,
    )

    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
                            T,       # Temperature of heat bath
                            1.0/u.picoseconds,  # Friction coefficient
                            1.0*u.femtoseconds, # Time step
    )

    # barostat
    barostat = mm.MonteCarloBarostat(1.0*u.atmosphere, T)
    barostat_id = system.addForce(barostat)

    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
    # the platform to use the default (fastest) platform
    platform = mm.Platform.getPlatformByName(platform_name)
    prop = dict(CudaPrecision='mixed') # Use mixed single/double precision

    # Create the Simulation object
    sim = app.Simulation(complex_solv.topology, system, integrator, platform, prop)

    # Set the particle positions
    sim.context.setPositions(complex_solv.positions)

    # Minimize the energy
    if minimize_steps > 0:
        print('Minimizing energy')
        sim.minimizeEnergy(maxIterations=minimize_steps)

    # Set up the reporters to report energies and coordinates
    sim.reporters.append(
            StateDataReporter("log.ditopic_solv", 10000, step=True, potentialEnergy=True,
                              kineticEnergy=True, temperature=True, volume=True,
                              density=True)
    )
    sim.reporters.append(
            StateDataReporter(sys.stdout, 10000, step=True, potentialEnergy=True,
                              kineticEnergy=True, temperature=True, volume=True,
                              density=True)
    )
    sim.reporters.append(NetCDFReporter('ditopic_solv.nc', 100000, crds=True))

    # Run NPT dynamics
    if NPT_steps > 0:
        print('Running NPT dynamics')
        sim.step(NPT_steps)
    system.removeForce(barostat_id)
    sim.context.reinitialize(preserveState=True)

    return sim


if not biased:
    sim = generate_simulation(input_files)

    # Run NVT dynamics
    print('Running NVT dynamics')
    sim.step(timesteps)
    sim.saveState('final.xml')
else:
    # Extract atom indexes from residue names
    sim = generate_simulation(input_files, NPT_steps=0, minimize_steps=0)
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
    cv = [Distance([indices_1, indices_2]), DistancesSum([indices_0, indices_1, indices_2])]
    grid = Grid(lower=(0.33, 1.0), upper=(1.93, 4.2), shape=(32, 32))
    cv_restraints = CVRestraints(lower=(0.33, 1.0), upper=(1.93, 4.2), ku=10, kl=10)
    sampling_method = ABF(cv, grid, restraints=cv_restraints)

    # Run biased dynamics and save results
    callback = ABFLogger("logger", period_hist_force=timesteps//10, period_cv=timesteps//1000)
    state = pysages.run(sampling_method, generate_simulation, timesteps, callback)
    result = pysages.analyze(state)
    energy = np.asarray(result["free_energy"])
    forces = np.asarray(result["mean_force"])
    grid = np.asarray(result["mesh"])
    np.savetxt("FES.txt", np.hstack([grid, energy.reshape(-1, 1)]))
    np.savetxt("Forces.txt", np.hstack([grid, forces.reshape(-1, grid.shape[1])]))
