from openmm.app import PDBFile, ForceField, Topology, Simulation, NoCutoff, Modeller
from openmm import NonbondedForce, VerletIntegrator
from openmm.unit import kilocalorie_per_mole, kilojoule_per_mole, nanometer, elementary_charge
from Bio.PDB import PDBIO  # ✅ Used to write Biopython structure to temp PDB
from openmm.app import PDBFile
import numpy as np
import tempfile  # ✅ Creates a temporary PDB file for OpenMM conversion

from openmm import Vec3
from openmm.unit import nanometer
from scipy.spatial.transform import Rotation as R  # ✅ For applying rotations
import os

from openmm.app import PDBFile

def save_pdb_in_angstroms(topology, positions, output_pdb):
    """Save an OpenMM PDB file with coordinates in Ångstroms instead of nanometers."""
    
    # ✅ Convert positions from nanometers (`nm`) to Ångstroms (`Å`)
    angstrom_positions = [Vec3(pos[0] * 10, pos[1] * 10, pos[2] * 10) for pos in positions]

    # ✅ Save the PDB file with converted coordinates
    with open(output_pdb, "w") as pdb_file:
        PDBFile.writeFile(topology, angstrom_positions, pdb_file)
    
    print(f"📂 Saved PDB with Ångstrom coordinates: {output_pdb}")

# Example usage:
# save_pdb_in_angstroms(modeller.topology, modeller.positions, "output_angstroms.pdb")

class EnergyCalculator:
    def __init__(self, biopython_structure):

        """Initialize EnergyCalculator by converting a Biopython Structure to OpenMM, adding hydrogens, and computing energies."""
        self.structure = biopython_structure  # ✅ Store Biopython structure

        # ✅ Step 1: Convert Biopython Structure to OpenMM
        self.topology, self.pdbfile = self._convert_to_openmm(biopython_structure)

        # ✅ Step 2: Add Missing Hydrogens
        self.modeller = Modeller(self.topology, self.pdbfile.positions)
        self.forcefield = ForceField('charmm36.xml')
        self.modeller.addHydrogens(self.forcefield)  # ✅ Adds hydrogens automatically

        # ✅ Step 3: Load OpenMM System with Hydrogens
        self.system = self.forcefield.createSystem(self.modeller.topology, nonbondedMethod=NoCutoff)

        # ✅ Step 4: Create OpenMM Simulation
        self.simulation = Simulation(self.modeller.topology, self.system, VerletIntegrator(0.001))
        self.simulation.context.setPositions(self.modeller.positions)

        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions()

    def _convert_to_openmm(self, biopython_structure):
        """Convert a Biopython Structure to an OpenMM Topology & PDBFile."""
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb:
            # ✅ Step 1: Write Biopython structure to a temporary PDB file
            io = PDBIO()
            io.set_structure(biopython_structure)
            io.save(temp_pdb.name)

            # ✅ Step 2: Load PDB in OpenMM
            openmm_pdb = PDBFile(temp_pdb.name)
            topology = openmm_pdb.topology

            # ✅ Extract atom positions (BEFORE setting up the system)
            initial_positions = openmm_pdb.positions  # ✅ OpenMM stores positions

        return topology, openmm_pdb

    def compute_pairwise_nonbonded_energy(self, intermolecular_only=False):
        """Compute non-bonded energy (Lennard-Jones + Electrostatics) **per residue pair**."""
        print(f"🔍 Computing {'intermolecular' if intermolecular_only else 'total'} pairwise non-bonded energies...")

        # ✅ Extract OpenMM NonbondedForce
        nonbonded_force = None
        for i in range(self.system.getNumForces()):
            force = self.system.getForce(i)
            if isinstance(force, NonbondedForce):
                nonbonded_force = force
                break

        if nonbonded_force is None:
            raise ValueError("❌ No NonbondedForce found in the system.")

        # ✅ Extract energy state from OpenMM
        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)

        # ✅ Extract all residues from all chains into a single list
        all_residues = [res for chain in self.modeller.topology.chains() for res in chain.residues()]
        
        residue_pairs = {}
        # ✅ Loop through all residue pairs, including different chains
        for i, res1 in enumerate(all_residues):
            for j, res2 in enumerate(all_residues):
                ####
                ##continue
                if i >= j:
                    continue  # ✅ Avoid double-counting and self-interactions
        
                if intermolecular_only and res1.chain.id == res2.chain.id:
                    continue  # ✅ Skip interactions within the same chain if intermolecular_only is True
        
                # ✅ Compute energy for this residue pair
                pair_energy = self.compute_residue_energy(res1, res2, positions, nonbonded_force)
                residue_pairs[((res1.name, res1.id), (res2.name, res2.id))] = pair_energy

        return residue_pairs

    def compute_residue_energy(self, res1, res2, positions, nonbonded_force):
        """Compute energy contribution between two residues using OpenMM."""
        pair_energy = 0.0

        for atom1 in res1.atoms():
            for atom2 in res2.atoms():
                q1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(atom1.index)
                q2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(atom2.index)

                # ✅ Convert to correct units
                q1 = q1.value_in_unit(elementary_charge)
                q2 = q2.value_in_unit(elementary_charge)
                sigma1 = sigma1.value_in_unit(nanometer)
                sigma2 = sigma2.value_in_unit(nanometer)
                epsilon1 = epsilon1.value_in_unit(kilojoule_per_mole)
                epsilon2 = epsilon2.value_in_unit(kilojoule_per_mole)

                # ✅ Compute distance
                pos1 = positions[atom1.index]._value
                pos2 = positions[atom2.index]._value
                r_ij = np.linalg.norm(pos1 - pos2)

                if r_ij == 0:
                    continue  # ✅ Skip self interactions

                r_ij *= 10

                # ✅ Electrostatic Energy
                k_elec = 138.935485  # Electrostatic constant in kJ/mol·nm·e²
                E_elec = (k_elec * q1 * q2) / (4 * r_ij * r_ij)

                # ✅ Lennard-Jones Energy
                sigma_ij = (sigma1 + sigma2) / 2
                epsilon_ij = np.sqrt(epsilon1 * epsilon2)
                E_LJ = 4 * epsilon_ij * ((sigma_ij / r_ij)**12 - (sigma_ij / r_ij)**6)

                # ✅ Sum contributions and convert to kcal/mol
                pair_energy += (E_elec + E_LJ) / 4.184  

        return pair_energy
    
    def monte_carlo_sampling(self, moving_chain_id, num_steps=4000, max_translation=0.5, max_rotation=10, temp_initial=3000, decay_constant=1000, output_dir="data"):
        """Perform Monte Carlo sampling and save PDB files at each step."""
        print(f"🔄 Starting Monte Carlo: Moving Chain {moving_chain_id}, {num_steps} steps...")
    
        kB = 0.0019872041  # Boltzmann constant in kcal/(mol*K)
    
        # ✅ Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
    
        # ✅ Get Initial Energy
        current_energy = sum(self.compute_pairwise_nonbonded_energy(intermolecular_only=True).values())
    
        state = self.simulation.context.getState(getPositions=True)
        positions = state.getPositions()  # ✅ OpenMM Vec3 format
    
        # ✅ Convert OpenMM positions to a NumPy array (NO UNITS)
        positions_array = np.array([[pos[0].value_in_unit(nanometer),  
                                     pos[1].value_in_unit(nanometer),  
                                     pos[2].value_in_unit(nanometer)] for pos in positions])
    
        # ✅ Get Atom Indices for the Moving Chain
        moving_chain_atoms = [atom.index for chain in self.modeller.topology.chains() if chain.id == moving_chain_id for atom in chain.atoms()]
    
        # ✅ Compute Initial Center of Mass (COM)
        com = np.mean(positions_array[moving_chain_atoms], axis=0)
    
        for step in range(num_steps):
            # ✅ Generate random translation
            translation = np.random.uniform(-max_translation, max_translation, 3)
    
            # ✅ Generate random rotation matrix
            rotation_angles = np.random.uniform(-max_rotation, max_rotation, 3)  # ✅ Rotation in degrees
            rotation_matrix = R.from_euler('xyz', rotation_angles, degrees=True).as_matrix()
    
            # ✅ Apply Transformations **Only to Moving Chain**
            perturbed_positions = positions_array.copy()
            for atom_index in moving_chain_atoms:
                original_pos = positions_array[atom_index]
    
                # ✅ Move to COM Frame (Correct Rotation)
                relative_pos = original_pos - com  
    
                # ✅ Apply Rotation **without distortion**
                rotated_pos = rotation_matrix @ relative_pos  
    
                # ✅ Apply Translation **after** rotation
                new_pos = rotated_pos + com + translation  
    
                # ✅ Debugging: Print Movement of First Atom
                if step % 500 == 0 and atom_index == moving_chain_atoms[0]:  
                    print(f"🔬 Step {step}: Atom {atom_index} moved from {original_pos} to {new_pos}")
    
                perturbed_positions[atom_index] = new_pos  
    
            # ✅ Step 0 Check: Detect & Fix Shrinking
            if step == 0:
                initial_com = np.mean(positions_array[moving_chain_atoms], axis=0)
                new_com = np.mean(perturbed_positions[moving_chain_atoms], axis=0)
                print(f"✅ Initial COM: {initial_com}")
                print(f"✅ New COM after first move: {new_com}")
    
                # 🚨 If COM shifts too much, reset positions
                if np.linalg.norm(new_com - initial_com) > max_translation * 5:  
                    print("⚠️ WARNING: Step 0 transformation is too large! Resetting...")
                    perturbed_positions = positions_array.copy()  # ✅ Restore original positions
    
            # ✅ Convert back to OpenMM Vec3 before setting positions
            perturbed_positions_openmm = [Vec3(*pos) * nanometer for pos in perturbed_positions]
            self.simulation.context.setPositions(perturbed_positions_openmm)
    
            # ✅ Compute New Energy
            new_energy = sum(self.compute_pairwise_nonbonded_energy(intermolecular_only=True).values())
    
            # ✅ Compute Acceptance Probability
            temperature = temp_initial * np.exp(-step / decay_constant)
            delta_E = new_energy - current_energy
            acceptance_prob = min(1, np.exp(-delta_E / (kB * temperature)))
    
            # ✅ Accept or Reject Move Based on Energy Score
            if np.random.rand() < acceptance_prob:
                positions_array = perturbed_positions  
                current_energy = new_energy  
                print(f"✅ Step {step}: Accepted (ΔE = {delta_E:.3f} kcal/mol, New Energy = {new_energy:.3f})")
    
                # ✅ Convert positions to unitless before writing PDB
                pdb_filename = os.path.join(output_dir, f"step_{step:04d}.pdb")
                save_pdb_in_angstroms(self.modeller.topology, perturbed_positions, pdb_filename)
                print(f"📂 Saved PDB: {pdb_filename}")
    
            else:
                self.simulation.context.setPositions([Vec3(*pos) * nanometer for pos in positions_array])  
                print(f"❌ Step {step}: Rejected (ΔE = {delta_E:.3f} kcal/mol, Energy Unchanged)")
    
        print("🏁 Monte Carlo simulation complete!")
    
