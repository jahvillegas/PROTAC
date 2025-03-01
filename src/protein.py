from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB import Atom as PDBAtom
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
from openmm import Platform
import random
platform = Platform.getPlatformByName('CPU')
platform.setPropertyDefaultValue('DeterministicForces', 'true')

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

