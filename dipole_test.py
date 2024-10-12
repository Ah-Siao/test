import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
import math
import numpy as np

with open('smiles_127_conformers.pkl', 'rb') as f:
    data = pickle.load(f)

data = data[:10]
dipole_list = []
total = np.array([])
for mol in data:
    xyz = Chem.MolToXYZBlock(mol)
    psi4.core.set_output_file('output.dat', False)
    psi4.geometry(xyz)
    psi4.set_options({'basis': "6-31G*"})
    # Perform the calculation and get the dipole moment
    energy, wfn = psi4.energy('hf', return_wfn=True)
    dipole_moment = wfn.variable('SCF DIPOLE')
    print(f"Dipole moment :{dipole_moment} Debye")
    x, y, z = dipole_moment
    dipole = math.sqrt(x**2+y**2+z**2)
    print(x, y, z, dipole)
    new_info = np.array([x, y, z, dipole])
    if not total.any():
        total = new_info
    else:
        total = np.vstack([total, new_info])
    np.save('dipoles.npy', total)
