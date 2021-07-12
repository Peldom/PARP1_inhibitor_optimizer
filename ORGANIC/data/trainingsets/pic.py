from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
# import sys
# sys.path.append('../model')

#from nn_metrics import KerasNN
i=0
results=[]
with open(r'parp_4.smi',mode='r') as f:
    for line in f:

        # print(line)
        try:
            mol1 = Chem.MolFromSmiles(line)
            mols = mol1
            img = Draw.MolsToImage([mols], subImgSize=(512, 512))
            img.save("./parp_4/{}.png".format(i))
            i+=1
            results.append(Chem.MolToSmiles(mols))

        except(ValueError):pass
print(results)
#
# mol1 = Chem.MolFromSmiles('Fc1c(C(=O)N2CCN(C(=O)c3cnccc3)CC2)cc(CC2=NNC(=O)C3=C2CCCC3)cc1')
# mols = mol1
# img = Draw.MolsToImage([mols], subImgSize=(500, 500))
# img.save("{}.png".format('1'))
#
