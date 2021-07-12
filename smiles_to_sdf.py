from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem

import re

file_object=open("decoys.smi")
w=Chem.SDWriter("test.sdf")
for line in file_object.readlines() :
    [smi,zid,pid]=line.split()
    m=Chem.MolFromSmiles(smi)
    #print(Chem.MolToMolBlock(m)) 
    m.SetProp("_Name",zid)
    AllChem.Compute2DCoords(m)
   
    w.write(m)
    