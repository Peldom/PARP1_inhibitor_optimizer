from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem

import re

file_object=open(r"D:\360Downloads\操作页面\ORGANIC\data\trainingsets\4_parp.smi")
w=Chem.SDWriter("parp_4.sdf")
for line in file_object.readlines() :
    # [smi,zid,pid]=line.split()
    m=Chem.MolFromSmiles(line)
    #print(Chem.MolToMolBlock(m)) 
    # m.SetProp("_Name",zid)
    AllChem.Compute2DCoords(m)
   
    w.write(m)
    