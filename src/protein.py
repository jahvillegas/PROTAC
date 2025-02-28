from openmm.app import PDBFile, Modeller, ForceField
from Bio.PDB import PDBParser
from pdbfixer import PDBFixer

class Protein:
    def __init__(self, pdb_file: str):
        """
        Load a protein structure from a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
        """
        self.pdb_file = pdb_file
        self.openmm_pdb = PDBFile(pdb_file)  # Load structure
        self.fixed_pdb_file = None  # Will be set if fix_missing_atoms() is called

    def fix_missing_atoms(self, forcefield="amber14-all.xml", heavy_atoms=True):
        """
        Fix missing atoms in the PDB file using PDBFixer (for heavy atoms) and OpenMM Modeller (for hydrogens).

        Args:
            forcefield (str): OpenMM force field XML file for hydrogen addition.
            heavy_atoms (bool): If True, fixes missing side chains and heavy atoms.
        """
        print("🔧 Fixing missing atoms...")

        # Use PDBFixer to rebuild heavy atoms
        fixer = PDBFixer(self.pdb_file)
        fixer.findMissingResidues()  # Detect missing residues
        fixer.findMissingAtoms()  # Detect missing heavy atoms
        fixer.addMissingAtoms()  # Add missing heavy atoms

        # Save the fixed PDB file with heavy atoms
        self.fixed_pdb_file = self.pdb_file.replace(".pdb", "_fixed.pdb")
        with open(self.fixed_pdb_file, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"✅ Fixed PDB saved as: {self.fixed_pdb_file}")

        # Use OpenMM Modeller to add missing hydrogens
        forcefield = ForceField(forcefield, "amber14/tip3p.xml")
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.addHydrogens(forcefield)

        # Save the final fixed PDB file with hydrogens
        self.fixed_pdb_file = self.fixed_pdb_file.replace(".pdb", "_h.pdb")
        with open(self.fixed_pdb_file, "w") as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

        print(f"✅ Added missing hydrogens: {self.fixed_pdb_file}")

    def get_openmm_pdb(self):
        """
        Returns the OpenMM PDB object. Uses the fixed PDB if available.

        Returns:
            PDBFile: OpenMM-compatible PDB object.
        """
        if self.fixed_pdb_file:
            return PDBFile(self.fixed_pdb_file)
        return self.openmm_pdb

    def print_summary(self):
        """
        Prints a summary of the structure.
        """
        pdb_file_to_parse = self.fixed_pdb_file if self.fixed_pdb_file else self.pdb_file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file_to_parse)

        num_atoms = len([atom for atom in structure.get_atoms()])
        num_residues = len([res for res in structure.get_residues()])
        chains = [chain.get_id() for chain in structure.get_chains()]

        print(f"🔬 Number of atoms: {num_atoms}")
        print(f"🔬 Number of residues: {num_residues}")
        print(f"🔬 Chains: {', '.join(chains)}")

