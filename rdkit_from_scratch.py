from rdkit import Chem

# Create an empty molecule or modify an existing one
rw_mol = Chem.RWMol()


# Add atoms (Hydrogen and Oxygen in this case)
h_atom = Chem.Atom(1)  # Hydrogen
o_atom = Chem.Atom(8)  # Oxygen

# Get the indices of the atoms added
h_idx = rw_mol.AddAtom(h_atom)
o_idx = rw_mol.AddAtom(o_atom)

# Add a bond between Hydrogen and Oxygen (single bond)
rw_mol.AddBond(h_idx, o_idx, Chem.BondType.SINGLE)

# Convert the RW molecule back to a standard RDKit molecule
mol = rw_mol.GetMol()

# Print the SMILES representation of the molecule
print(Chem.MolToSmiles(mol))
