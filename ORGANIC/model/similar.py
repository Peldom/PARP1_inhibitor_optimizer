from rdkit import Chem
from rdkit.Chem import AllChem as ch
from rdkit.Chem import Draw
from rdkit import DataStructs


with open('dgl.smi','r') as f:
    mols = [Chem.MolFromSmiles(x) for x in f if x is not None]
len(mols)
def inquire(mol):
    #读入查询分子，计算指纹
    nicotine = ch.MolFromSmiles(mol)
    nicotine_fingerprint = ch.GetMorganFingerprint(nicotine, 2)
    return nicotine_fingerprint
l=['c1ccc2c(c1)c(n[nH]c2=O)Cc3ccc(c(c3)C(=O)N4CCN(CC4)C(=O)C5CC5)F','CNCc1ccc(cc1)c2c3c4c(cc(cc4[nH]2)F)C(=O)NCC3','c1cc2cn(nc2c(c1)C(=O)N)c3ccc(cc3)[C@@H]4CCCNC4','Cn1c(ncn1)[C@H]2c3c4c(cc(cc4N[C@@H]2c5ccc(cc5)F)F)c(=O)[nH]n3']
def func(s):
    #读入查询分子，计算指纹
    nicotine = ch.MolFromSmiles(s)
    nicotine_fingerprint = ch.GetMorganFingerprint(nicotine, 2)

    #计算分子库每个分子指纹
    mols_fps = [(m, ch.GetMorganFingerprint(m, 2)) for m in mols]

    # 计算相似度并排序，输出最相似的的前20个分子

    mols_nicotinesim = [(m, DataStructs.TanimotoSimilarity(fp, nicotine_fingerprint))
                        for m, fp in mols_fps]
    sorted_mols_nicotinesim = sorted(mols_nicotinesim, key=lambda x: x[1], reverse=True)
    for i in range(len(mols)):
        # print(s,ch.MolToSmiles(mols[i]),mols_nicotinesim[i][1])
        similar = mols_nicotinesim[i][1]
        print(similar)
for i in range(len(l)):
    func(l[i])