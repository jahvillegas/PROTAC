import sys
import os
# Ensure Python can find src/ directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from protein import *

# Step 1: Fix missing side chains in-memory
protein = Protein("data/incomplete.pdb", "data/aa_template.pdb")
fixed_structure = protein.copy_side_chains()  # Returns modified structure

# Step 2: Save the fixed structure to a PDB file
protein.save_pdb(fixed_structure, "data/fixed_protein.pdb")

# Step 3: Compute CHARMM energy using OpenMM (directly from structure)
calculator = EnergyCalculator(fixed_structure)
energy = calculator.compute_energy()

print(f"✅ Final Potential Energy: {energy}")
