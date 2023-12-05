#!/usr/bin/env python
#from __future__ import division, print_function

import sys

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


biased = False
timesteps = 10000000


input_files = ('../force_field/complex_solv.prmtop', '../force_field/complex_solv.rst7')
T = 298 * u.kelvin


def generate_simulation(input_files, T=T, platform_name="CUDA"):
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
    print('Minimizing energy')
    sim.minimizeEnergy(maxIterations=500)
    
    # Set up the reporters to report energies and coordinates
    sim.reporters.append(
            StateDataReporter("log.comlex_solv", 10000, step=True, potentialEnergy=True,
                              kineticEnergy=True, temperature=True, volume=True,
                              density=True)
    )
    sim.reporters.append(
            StateDataReporter(sys.stdout, 10000, step=True, potentialEnergy=True,
                              kineticEnergy=True, temperature=True, volume=True,
                              density=True)
    )
    sim.reporters.append(NetCDFReporter('complex_solv.nc', 100000, crds=True))

    # Run NPT dynamics
    print('Running NPT dynamics')
    sim = generate_simulation(input_files)
    sim.step(100000)
    system.removeForce(barostat_id)
    sim.context.reinitialize(preserveState=True)

    return sim



if not biased:
    # Run NVT dynamics
    print('Running NVT dynamics')
    sim.step(timesteps)
    sim.saveState('final.xml')
else:
    # add logic to extract indices from res names
    indices_1 = [1]
    indices_2 = [2]
    cv = Distance([indices_1, indices_2])
    grid = Grid(lower=(0.9,), upper=(3.0,), size=(32,))
    cv_restraints = CVRestraints(lower=(0.9,), upper=(3.0,), ku=10, kl=10)
    sampling_method = ABF(cv, grid, restraints=cv_restraints)
    state = pysages.run(sampling_method, generate_simulation, timesteps)
    result = pysages.analyze(state)
    fe = result["free_energy"]
    mesh = result["mesh"]
    hist = result["histogram"]
    # plt.plot(mesh, fe)
    # plt.plot(mesh, hist)