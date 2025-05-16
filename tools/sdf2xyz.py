#!/usr/bin/env python3
import sys
from rdkit import Chem

def sdf_to_xyz(sdf_file, xyz_file):
  suppl = Chem.SDMolSupplier(sdf_file)
  with open(xyz_file, 'w') as xyz:
    for mol in suppl:
      if mol is None:
        continue
      atoms = mol.GetAtoms()
      coords = mol.GetConformers()[0].GetPositions()
      xyz.write(f"{len(atoms)}\n\n")
      for atom, coord in zip(atoms, coords):
        symbol = atom.GetSymbol()
        x, y, z = coord
        xyz.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")

def main():
  if len(sys.argv) != 3:
    print("Usage: sdf2xyz.py <input.sdf> <output.xyz>")
    sys.exit(1)

  sdf_file = sys.argv[1]
  xyz_file = sys.argv[2]

  try:
    sdf_to_xyz(sdf_file, xyz_file)
    print(f"Converted {sdf_file} to {xyz_file}")
  except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

if __name__ == "__main__":
  main()