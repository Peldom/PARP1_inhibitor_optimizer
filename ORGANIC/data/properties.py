from scopy.druglikeness import molproperty
from scopy import ScoConfig
from rdkit.Chem import Descriptors
from rdkit import Chem
import os
mol = Chem.MolFromMolFile('test.sdf')
# props =
# print(props)

# print(mol)
with open('test.sdf','r')as f:
    for s in f.readlines:
        print(s)