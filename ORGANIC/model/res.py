import sys

import matplotlib.pyplot as plt
from rdkit.Chem import Draw
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
def get_violations(files):
    v=[]
    for m in smi2mol(files):
        try:
            v.append(num_lipinski_violations(m, ruleset='extended'))

        except:pass
    violations=v
    return violations
# print(float(sum(violations)) / len(violations))
def ls_result(smipath,violations):
    mols=[]
    lip=[]
    with open(smipath,'r') as f:
        for mol in f:
            mols.append(mol)
    # print(len(mols))
    b=dict(zip(mols,violations))
    for k in list(b.keys()):
        if b[k] == 0:
            lip.append(k)
    return lip
    # with open('metrics_results.smi','w') as f:
    #     for k in list(b.keys()):
    #         if b[k]==0:
    #             f.write(k)
            # elif b[k]==1:pass
            # else:print(k)

def res_(files,list):
    with open(files,mode='w',encoding='utf-8') as g:
        # with open(r'metrics_results.smi',mode='r') as f:
        for line in list:
            # print(line)
            try:
                mol1 = Chem.MolFromSmiles(line)
                mols = mol1
                img = Draw.MolsToImage([mols], subImgSize=(512, 512))
                # img.save("{}.png".format(i))
                # i+=1
                g.write(Chem.MolToSmiles(mols)+'\n')
            except(ValueError):pass



def res(files,list):
    with open(files,mode='a+',encoding='utf-8') as g:
        # with open(r'metrics_results.smi',mode='r') as f:
        for line in list:
            # print(line)
            try:
                mol1 = Chem.MolFromSmiles(line)
                mols = mol1
                img = Draw.MolsToImage([mols], subImgSize=(512, 512))
                # img.save("{}.png".format(i))
                # i+=1
                g.write(Chem.MolToSmiles(mols)+'\n')
            except(ValueError):pass
import os
def list_set(path):
    # path='test1.txt'
    l=[]
    with open(path,'r')as f:
        for line in f:
            l.append(line.strip())

    l=list(set(l))
    # print(l)

    with open(path,'w')as f:
        for k in l:
            f.write(k+'\n')
if __name__ == '__main__':

    # os.system("python3 organic.py")
    # print('++++++++++++++++第{}次循环++++++++++++++++'.format(j))
    for i in range(10):

        violations = get_violations('./epoch_data/metrics_{}.smi'.format(i))
        lip = ls_result('./epoch_data/metrics_{}.smi'.format(i), violations)
        # res('./epoch_data/metrics_{}.smi'.format(i),lip)
        res('success_4_.smi', lip)
        res_('success.smi',lip)
        list_set('success_4_.smi')
        list_set('success.smi')
        with open('../data/trainingsets/Value_Ki.smi',mode='a+',encoding='utf-8')as h:
            with open('success.smi', mode='r', encoding='utf-8')as j:
                for line in j:
                    h.write(line)

        # list_set('../data/trainingsets/Value_Ki.smi')
        os.remove('./epoch_data/metrics_{}.smi'.format(i))
        print('################第{}次xiao循环####################'.format(i))
    violations = get_violations('../data/trainingsets/Value_Ki.smi'.format(i))
    lip = ls_result('../data/trainingsets/Value_Ki.smi'.format(i), violations)
    res('../data/trainingsets/Value_Ki.smi', lip)

# if __name__ == '__main__':
#     violations=get_violations('model/epoch_data/lipinski_Valueki_0.smi')
#     lip=ls_result('model/epoch_data/lipinski_Valueki_0.smi',violations)
#     res('model/epoch_data/lipinski_Valueki_0.smi',lip)
#     res('success_4_.smi',lip)