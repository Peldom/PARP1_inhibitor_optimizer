import sys

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors as rdescriptors
import numpy as np
def num_hydrogen_bond_acceptors(mol):
    return rdescriptors.CalcNumLipinskiHBA(mol)

def num_hydrogen_bond_donors(mol):
    return rdescriptors.CalcNumLipinskiHBD(mol)

def MW(mol):
    return Descriptors.MolWt(mol)

def logP(mol):
    return Descriptors.MolLogP(mol)

def TPSA(mol):
    return Descriptors.TPSA(mol)

def num_rotatable_bonds(mol):
    return Descriptors.NumRotatableBonds(mol)

def num_heavy_atoms(mol):
    return mol.GetNumHeavyAtoms()

### Lipinski checks, return True/False
def lipinski_HBA(mol):
    return num_hydrogen_bond_acceptors(mol) <= 10

def lipinski_HBD(mol):
    return num_hydrogen_bond_acceptors(mol) <= 5

def lipinski_MW(mol):
    return MW(mol) < 500

def lipinski_logP(mol):
    return logP(mol) <= 5

def lipinski_TPSA(mol):
    return TPSA(mol) < 140

def lipinski_rotatable_bonds(mol):
    return num_rotatable_bonds(mol) < 10

def lipinski_heavy_atoms(mol):
    return 20 < num_heavy_atoms(mol) < 70
### Lipinski rule violation counter
def num_lipinski_violations(mol, *args, ruleset='basic'):
    rulesets = {'basic': (lipinski_HBA, lipinski_HBD, lipinski_logP, lipinski_MW),
                'extended': (lipinski_HBA, lipinski_HBD, lipinski_logP, lipinski_MW, lipinski_TPSA, lipinski_heavy_atoms, lipinski_rotatable_bonds)}
    rulescount = len(rulesets[ruleset])
    # count of all rules - count of passed rules = count of failed rules
    # print(rulescount,sum([f(mol) for f in rulesets[ruleset]]))
    return rulescount - sum([f(mol) for f in rulesets[ruleset]])
def smi2mol(smipath):
    # smipath='drugs_sub.smi'
    with open(smipath,'r') as f:
        for mol in f:
            if mol:
                # yield Chem.MolFromSmiles(mol)
                try:
                     yield Chem.MolFromSmiles(mol)
                except Exception:pass




# violations = [num_lipinski_violations(m, ruleset='extended') for m in smi2mol('lipinski_Valueki_0.smi')]
# print(float(sum(violations)) / len(violations))
# print(float(sum(violations)))
# violations
v=[]
for m in smi2mol('metrics_0.smi'):
    try:
        v.append(num_lipinski_violations(m, ruleset='extended'))

    except:pass
violations=v
print(float(sum(violations)) / len(violations))
def ls_result(smipath):
    mols=[]
    with open(smipath,'r') as f:
        for mol in f:
            mols.append(mol)
    print(len(mols))
    b=dict(zip(mols,violations))
    with open('metrics_results.smi','w') as f:
        for k in list(b.keys()):
            if b[k]==0:
                f.write(k)
            # elif b[k]==1:pass
            # else:print(k)

ls_result('metrics_0.smi')












