from rdkit import Chem
from rdkit.Chem import Draw
import time
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
# import sys
# sys.path.append('../model')

#from nn_metrics import KerasNN
# i=0
# for m in range(4):
#     with open(r'J:\GAN\新数据\ORGANIC-1(复件)\model\epoch_data\metrics-logP_opv_sub.smi_{}.smi'.format(m),mode='r') as f:
#         for line in f:
#
#             # print(line)
#             try:
#                 mol1 = Chem.MolFromSmiles(line)
#                 mols = mol1
#                 img = Draw.MolsToImage([mols], subImgSize=(512, 512))
#                 # img.save("{}.png".format(i))
#                 i+=1
#                 # print(i)
#             except(ValueError):pass
#     print(i)
#     print('这是第{}个文件'.format(m+1))
#     time.sleep(5)



i=0
with open(r'D:\360Downloads\ORGANIC-1(复件)\model\epoch_data\lipinski_Valueki_0.smi',mode='r') as f:
    for line in f:

        # print(line)
        try:
            mol1 = Chem.MolFromSmiles(line)
            mols = mol1
            img = Draw.MolsToImage([mols], subImgSize=(512, 512))
            # img.save("{}.png".format(i))
            i+=1
            # print(i)
        except(ValueError):pass
print(i)
# time.sleep(5)






# mol1 = Chem.MolFromSmiles('CCC1NC3CN(C#C)C1')
# mols = mol1
# img = Draw.MolsToImage([mols], subImgSize=(200, 200))
# img.save("{}.png".format('line'))