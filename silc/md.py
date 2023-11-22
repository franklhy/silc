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

# Load the Amber files
print('Loading AMBER files...')
complex_solv = load_file('../force_field/complex_solv.prmtop', '../force_field/complex_solv.rst7')

# Create the OpenMM system
print('Creating OpenMM System')
system = complex_solv.createSystem(nonbondedMethod=app.PME,
                                nonbondedCutoff=8.0*u.angstroms,
                                constraints=app.HBonds,
)

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        298*u.kelvin,       # Temperature of heat bath
                        1.0/u.picoseconds,  # Friction coefficient
                        1.0*u.femtoseconds, # Time step
)

# barostat
barostat = mm.MonteCarloBarostat(1.0*u.atmosphere, 298*u.kelvin)
barostat_id = system.addForce(barostat)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CUDA')
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
sim.step(100000)

# Run NVT dynamics
system.removeForce(barostat_id)
sim.context.reinitialize(preserveState=True)
print('Running NVT dynamics')
sim.step(10000000)
sim.saveState('final.xml')
