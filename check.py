import pickle
from rdkit import Chem

with open('test_Drugs.pkl', 'rb') as f:
    data = pickle.load(f)


new_data = [data[0]]
first_smiles = Chem.MolToSmiles(data[0])
for i in range(1, len(data)):
    smiles = Chem.MolToSmiles(data[i])
    if smiles == first_smiles:
        new_data.append(data[i])
    else:
        break

with open(f'smiles_{len(new_data)}_conformers.pkl', 'wb') as f:
    pickle.dump(new_data, f)
