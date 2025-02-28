from openmm.app import ForceField, Simulation
from openmm import VerletIntegrator, unit
from openmm.app import NoCutoff
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from src.protein import Protein  # Import our Protein class

class EnergyCalculator:
    def __init__(self, protein: Protein, forcefield="amber14-all.xml"):
        """
        Initialize the CHARMM-based energy calculator.

        Args:
            protein (Protein): Instance of the Protein class.
            forcefield (str): OpenMM force field XML file.
        """
        self.pdb = protein.get_openmm_pdb()  # Get the (fixed) structure
        self.forcefield = ForceField(forcefield)
        self.system = self.forcefield.createSystem(self.pdb.topology, nonbondedMethod=unit.NoCutoff)
        self.integrator = VerletIntegrator(0.001 * unit.picoseconds)
        self.simulation = Simulation(self.pdb.topology, self.system, self.integrator)

    def compute_energy(self):
        """
        Compute the potential energy of the system.

        Returns:
            float: Potential energy in kJ/mol.
        """
        self.simulation.context.setPositions(self.pdb.positions)
        state = self.simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        return energy.value_in_unit(unit.kilojoule_per_mole)

# Example Usage
if __name__ == "__main__":
    pdb_path = "./data/5t35_AB.pdb"  # Update with actual PDB file
    protein = Protein(pdb_path)

    # OPTIONALLY fix missing atoms
    user_choice = input("Do you want to fix missing atoms? (y/n): ").strip().lower()
    if user_choice == "y":
        protein.fix_missing_atoms()

    protein.print_summary()  # Show structure details

    calculator = EnergyCalculator(protein)
    energy = calculator.compute_energy()
    print(f"⚡ Potential Energy: {energy:.2f} kJ/mol")

