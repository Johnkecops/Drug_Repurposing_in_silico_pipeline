#!/usr/bin/env python3
"""
Script 4: Molecular Dynamics Simulation with OpenMM

Purpose: Run MD simulations to validate docking poses and assess binding stability.
Covers: System preparation, energy minimisation, equilibration, production MD.

Usage:
    python 04_run_md_simulation.py --pdb complex.pdb --output_dir md_output \\
                                   --production_time_ns 200
"""

import argparse
import os
from openmm.app import *
from openmm import *
from openmm.unit import *

def run_md_simulation(protein_pdb, output_dir='md_output', 
                      production_time_ns=200, verbose=True):
    """
    Run molecular dynamics simulation of protein-ligand complex.
    
    Args:
        protein_pdb: Protein structure (prepared PDB)
        output_dir: Output directory for trajectory and logs
        production_time_ns: Production MD duration in nanoseconds
        verbose: Print progress
    
    Returns:
        Trajectory filename, energy log filename
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        print(f"Setting up MD simulation ({production_time_ns} ns production)...")
    
    # Load structure
    pdb = PDBFile(protein_pdb)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    if verbose:
        print(f"Loaded {pdb.topology.getNumAtoms()} atoms, "
              f"{pdb.topology.getNumResidues()} residues")
    
    # Add hydrogens and solvate
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, model='tip3p', padding=10*angstroms,
                        ionicStrength=0.15*molar)
    
    if verbose:
        print(f"Solvated system: {modeller.topology.getNumAtoms()} atoms total")
    
    # Create system
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds,
        rigidWater=True
    )
    
    # Add barostat and thermostat
    system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Energy minimisation
    if verbose:
        print("Energy minimisation...")
    simulation.minimizeEnergy(maxIterations=1000)
    state = simulation.context.getState(getEnergy=True)
    if verbose:
        print(f"  Final energy: {state.getPotentialEnergy()/kilojoule_per_mole:.2f} kJ/mol")
    
    # NVT equilibration (50 ps)
    if verbose:
        print("NVT equilibration (50 ps)...")
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(25000)  # 50 ps at 2 fs/step
    
    # NPT equilibration (100 ps)
    if verbose:
        print("NPT equilibration (100 ps)...")
    simulation.step(50000)  # 100 ps
    
    # Production MD
    if verbose:
        print(f"Production MD ({production_time_ns} ns)...")
    
    trajectory_file = os.path.join(output_dir, 'trajectory.dcd')
    energy_log = os.path.join(output_dir, 'energy.log')
    
    simulation.reporters.append(
        DCDReporter(trajectory_file, 1000)  # Save every 2 ps
    )
    simulation.reporters.append(
        StateDataReporter(energy_log, 1000, step=True, time=True,
                         potentialEnergy=True, kineticEnergy=True,
                         temperature=True, density=True, progress=True)
    )
    
    # Calculate number of steps
    n_steps = int((production_time_ns * nanosecond) / (2 * femtosecond))
    simulation.step(n_steps)
    
    if verbose:
        print(f"\n✓ MD simulation complete!")
        print(f"  Trajectory: {trajectory_file}")
        print(f"  Energy log: {energy_log}")
    
    return trajectory_file, energy_log

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run MD simulation with OpenMM')
    parser.add_argument('--pdb', required=True, help='Input PDB file (protein+ligand)')
    parser.add_argument('--output_dir', default='md_output', help='Output directory')
    parser.add_argument('--production_time_ns', type=int, default=200,
                       help='Production MD duration (ns)')
    
    args = parser.parse_args()
    
    traj, energy = run_md_simulation(args.pdb, args.output_dir, 
                                     args.production_time_ns)
