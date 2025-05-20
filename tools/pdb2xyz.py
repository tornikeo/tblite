#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def pdb_to_xyz(pdb_file, xyz_file):
    """
    Convert a PDB file to an XYZ file.
    
    Parameters:
    pdb_file (str): Path to the input PDB file.
    xyz_file (str): Path to the output XYZ file.
    """
    # Read the PDB file
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    
    if mol is None:
        print(f"Error: Could not read PDB file {pdb_file}.")
        return
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(mol)
    if result != 0:  # Embedding failed
        print(f"Error: 3D coordinate generation failed for {pdb_file}.")
        return
    
    # Write to XYZ file
    with open(xyz_file, 'w') as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write("Generated from PDB\n")
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x} {pos.y} {pos.z}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python pdb2xyz.py <input.pdb> <output.xyz>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    xyz_file = sys.argv[2]
    
    pdb_to_xyz(pdb_file, xyz_file)
    print(f"Converted {pdb_file} to {xyz_file}.")