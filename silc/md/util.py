#!/usr/bin/env python
#from __future__ import division, print_function

import os
import sys
from glob import glob
from pathlib import Path

# OpenMM Imports
#import simtk.openmm as mm
#import simtk.openmm.app as app
import openmm as mm
import openmm.app as app

# ParmEd Imports
from parmed import load_file, unit as u
#from parmed.openmm import StateDataReporter


def generate_simulation(
        input_files, output_name=None, output_path=None, T=0, 
        minimize_steps=0, NPT_steps=0, NVT_steps=0, state_freq=10000, traj_freq=100000,
        platform_name="CUDA"):
    # Load the Amber files
    print('Loading AMBER files...')
    configuration = load_file(*input_files)

    # Create the OpenMM system
    print('Creating OpenMM System')
    system = configuration.createSystem(nonbondedMethod=app.PME,
                                    nonbondedCutoff=8.0*u.angstroms,
                                    constraints=app.HBonds,
    )

    # Create the integrator to do Langevin dynamics
    integrator = mm.LangevinIntegrator(
                            T,       # Temperature of heat bath
                            1.0/u.picoseconds,  # Friction coefficient
                            1.0*u.femtoseconds, # Time step
    )

    # Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
    # the platform to use the default (fastest) platform
    platform = mm.Platform.getPlatformByName(platform_name)
    prop = dict(CudaPrecision='mixed') # Use mixed single/double precision

    # Create the Simulation object
    sim = app.Simulation(configuration.topology, system, integrator, platform, prop)

    # Set the particle positions
    sim.context.setPositions(configuration.positions)

    # Minimize the energy
    if minimize_steps > 0:
        print('Minimizing energy')
        sim.minimizeEnergy(maxIterations=minimize_steps)

    if output_name is not None and output_path is not None:
        # Set up the reporters to report energies and coordinates
        sim.reporters.append(
            app.StateDataReporter(os.path.join(output_path, "log.%s" % output_name),
                            state_freq, step=True, potentialEnergy=True, kineticEnergy=True,
                            temperature=True, volume=True, density=True)
        )
        sim.reporters.append(
            app.StateDataReporter(sys.stdout, state_freq, step=True, potentialEnergy=True, kineticEnergy=True,
                            temperature=True, volume=True, density=True)
        )
        #sim.reporters.append(
        #    NetCDFReporter(os.path.join(output_path, '%s.nc' % output_name),
        #                   traj_freq, crds=True)
        #)
        sim.reporters.append(
            app.DCDReporter(os.path.join(output_path, '%s.dcd' % output_name),
                            traj_freq, enforcePeriodicBox=False)
        )

    # Run NPT dynamics
    if NPT_steps > 0:
        print('Running NPT dynamics')
        barostat = mm.MonteCarloBarostat(1.0*u.atmosphere, T)
        barostat_id = system.addForce(barostat)
        sim.step(NPT_steps)
        sim.saveState('final_NPT.xml')
        system.removeForce(barostat_id)
        sim.context.reinitialize(preserveState=True)
        print('Done with NPT dynamics')

    # Run NVT dynamics
    if NVT_steps > 0:
        print('Running NVT dynamics')
        sim.step(NVT_steps)
        sim.saveState('final_NVT.xml')
        print('Done with NVT dynamics')

    return sim

def prepare_restart(restart_name):
    output_path = "."
    restart_file = None
    restarts = sorted([int(p.split('_')[1]) for p in glob("restart*")])
    if len(restarts) == 0 and Path(restart_name).is_file():
        restart_file = Path(restart_name).resolve()
        output_path = "restart_1"
    elif len(restarts) > 0 and Path(os.path.join("restart_%d" % restarts[-1], restart_name)).is_file():
        restart_file = Path(os.path.join("restart_%d" % restarts[-1], restart_name)).resolve()
        output_path = "restart_%d" % (restarts[-1] + 1)

    if output_path != ".":
        os.mkdir(output_path)

    return restart_file, output_path