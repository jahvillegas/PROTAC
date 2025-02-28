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

# Define translation vector [dx, dy, dz]
translation_vector = np.array([10.0, 0.0, 0.0])  # Move chain B by 10 Å along X-axis

# Define rotation matrix (90-degree rotation around Z-axis)
theta = np.radians(90)
rotation_matrix = np.array([
    [np.cos(theta), -np.sin(theta), 0],
    [np.sin(theta),  np.cos(theta), 0],
    [0, 0, 1]
])

# Move only chain D
modified_structure = protein.displace_chain(protein.structure, chain_id="D", translation=translation_vector, rotation_matrix=rotation_matrix)

# ✅ Use modified_structure when saving
protein.save_pdb(modified_structure, "data/modified.pdb")  # ✅ Now compatible with your original function

calculator = EnergyCalculator(modified_structure)
energy = calculator.compute_energy()

print(f"✅ Final Potential Energy: {energy}")
