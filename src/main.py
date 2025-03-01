import sys
import os
# Ensure Python can find src/ directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from protein import *
from energy  import *

# Step 1: Fix missing side chains in-memory
protein = Protein("data/incomplete.pdb", "data/aa_template.pdb")
fixed_structure = protein.copy_side_chains()  # Returns modified structure


# Step 2: Save the fixed structure to a PDB file
protein.save_pdb(fixed_structure, "data/fixed_protein.pdb")

# ✅ Pass Biopython Structure to EnergyCalculator
energy_calculator = EnergyCalculator(fixed_structure)

# ✅ Compute Non-Bonded Energies Using OpenMM
#pairwise_energies = energy_calculator.compute_pairwise_nonbonded_energy(intermolecular_only=True)

#print("\n🔍 Sorted Pairwise Non-Bonded Energies:")
#sorted_energies = sorted(pairwise_energies.items(), key=lambda x: abs(x[1]), reverse=True)

#for (res1, res2), energy in sorted_energies[:20]:  # Show top 20 strongest interactions
#    print(f"{res1} ↔ {res2}: {energy:.3f} kcal/mol")

# ✅ Select Chain to Move (e.g., "A")
moving_chain = "D"

# ✅ Run Monte Carlo with Pairwise Non-Bonded Energy Scoring
energy_calculator.monte_carlo_sampling(moving_chain)

