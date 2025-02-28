from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB import Atom as PDBAtom
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np

forcefield_path = "/Users/josevillegas/Library/Python/3.9/lib/python/site-packages/openmm/app/data/"
forcefield = ForceField(forcefield_path + 'charmm36.xml', forcefield_path + 'charmm36/water.xml')

class Protein:
    def __init__(self, pdb_file, template_pdb):
        """Load a protein structure and a template PDB containing all 20 amino acids."""
        self.pdb_file = pdb_file
        self.template_pdb = template_pdb
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("protein", pdb_file)
        self.template_structure = self.parser.get_structure("template", template_pdb)

    def copy_side_chains(self):
        """Copy missing side chains and add OXT to the C-terminus."""
        print("🔄 Copying and aligning side chains from template PDB...")

        sup = Superimposer()

        for model in self.structure:
            for chain in model:
                residues = list(chain)  # Convert to a list for indexing
                last_residue = residues[-1]  # Get the C-terminal residue
                residue_name = last_residue.get_resname()
                residue_id = last_residue.get_id()[1]

                # Add OXT if missing
                if "OXT" not in [atom.get_name() for atom in last_residue]:
                    print(f"➕ Adding OXT to {residue_name} at {residue_id}")

                    # Get C-terminal C and O coordinates
                    c_atom = last_residue["C"].get_coord()
                    o_atom = last_residue["O"].get_coord()

                    # Estimate OXT position (project from C=O bond)
                    oxt_position = 2 * o_atom - c_atom  # Approximate OXT placement

                    # Create OXT atom
                    oxt_atom = PDBAtom.Atom("OXT", oxt_position, 1.0, 1.0, " ", "OXT", "O")

                    # Add to residue
                    last_residue.add(oxt_atom)

                for residue in chain:
                    residue_name = residue.get_resname()
                    residue_id = residue.get_id()[1]

                    # Find corresponding template residue
                    template_residue = None
                    for t_model in self.template_structure:
                        for t_chain in t_model:
                            for t_residue in t_chain:
                                if t_residue.get_resname() == residue_name:
                                    template_residue = t_residue
                                    break
                            if template_residue:
                                break
                        if template_residue:
                            break

                    if not template_residue:
                        print(f"⚠️ No template found for {residue_name} at {residue_id}, skipping.")
                        continue

                    # Extract only backbone atoms (N, CA, C) for alignment
                    backbone_atoms = ["N", "CA", "C"]
                    target_backbone = [atom for atom in residue if atom.get_name() in backbone_atoms]
                    template_backbone = [atom for atom in template_residue if atom.get_name() in backbone_atoms]

                    if len(target_backbone) != len(template_backbone):
                        print(f"❌ Backbone mismatch for {residue_name} {residue_id}, skipping.")
                        continue

                    # Align backbone atoms and apply transformation
                    sup.set_atoms(target_backbone, template_backbone)
                    sup.apply(template_residue)

                    # Copy missing side-chain atoms from the template
                    existing_atoms = {atom.get_name() for atom in residue}
                    for atom in template_residue:
                        if atom.get_name() not in backbone_atoms and atom.get_name() not in existing_atoms:
                            new_atom = atom.copy()
                            residue.add(new_atom)
                            print(f"✅ Copied {new_atom.get_name()} for {residue_name} at {residue_id}")

        print(f"✅ Side chains copied and aligned correctly.")
        return self.structure  # Return structure in memory

    def save_pdb(self, structure, output_pdb):
        """Save the modified structure to a PDB file."""
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        print(f"📂 Saved fixed structure as: {output_pdb}")

    def displace_chain(self, structure, chain_id, translation=None, rotation_matrix=None):
        """
        Mov e one chain while keeping the others fixed.

        Parameters:
        - structure: The PDB structure to modify
        - chain_id (str): Chain identifier to move (e.g., "B").
        - translation (numpy array): A vector [dx, dy, dz] for translation.
        - rotation_matrix (numpy array): A 3x3 rotation matrix.
    
        Returns:
        - structure: The modified PDB structure
        """
        print(f"🔄 Displacing chain {chain_id}...")

        model = next(iter(structure))  # Get the first model
        if chain_id not in model:
            print(f"❌ Chain {chain_id} not found! Skipping transformation.")
            return structure

        chain = model[chain_id]  # Get the chain to move

        for residue in chain:
            for atom in residue:
                old_coord = atom.get_coord()
                new_coord = old_coord

                # Apply translation if provided
                if translation is not None:
                    new_coord += translation  # Apply translation

                # Apply rotation if provided
                if rotation_matrix is not None:
                    new_coord = np.dot(rotation_matrix, new_coord)  # Apply rotation

                atom.set_coord(new_coord)

        print(f"✅ Chain {chain_id} displaced successfully!")
        return structure  # ✅ Return modified structure

class EnergyCalculator:
    def __init__(self, structure):
        """Initialize OpenMM and load the structure object directly."""
        self.structure = structure  # Store structure in memory
        self.pdb = self.structure_to_pdb(self.structure)  # Convert to PDB object

        # Load CHARMM force field
        self.forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')

        # Create a system
        self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
        self.modeller.addHydrogens(self.forcefield)

        self.system = self.forcefield.createSystem(self.modeller.topology, 
                                                   nonbondedMethod=NoCutoff, 
                                                   constraints=HBonds)

        self.integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        self.simulation = Simulation(self.modeller.topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.modeller.positions)

    def structure_to_pdb(self, structure):
        """Convert a Biopython structure object to an OpenMM PDBFile object (in memory)."""
        io = PDBIO()
        io.set_structure(structure)
        import io as sys_io
        buffer = sys_io.StringIO()
        io.save(buffer)
        buffer.seek(0)
        return PDBFile(buffer)

    def compute_energy(self):
        """Compute the potential energy of the system."""
        state = self.simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        print(f"⚡ Potential Energy: {energy}")
        return energy

