from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import RDConfig
import json

def getcore(ds_loc):
    #count
    suppl=Chem.SmilesMolSupplier(ds_loc)
    ms = [x for x in suppl if x] #要拆解的分子列表，mol
    cores=[] #母核列表 
    decorators=[]
    c_d={}
    for mol in suppl:
        cores.append(MurckoScaffold.GetScaffoldForMol(mol)) #返回mol
    for i in range(len(suppl)):
        if ms[i].HasSubstructMatch(cores[i]):
            rs = Chem.ReplaceCore(ms[i], cores[i]) #rs去除了母核
            res = Chem.GetMolFrags(rs, asMols=True) # res是侧链列表
            res_sm=[Chem.MolToSmiles(x) for x in res]
            decorators.append(res_sm) #返回一维列表,mol
        else:
            decorators.append([]) #错误，返回nul
        c_d[Chem.MolToSmiles(cores[i])]=decorators[i]
    # save
    save_lc=ds_loc[0:ds_loc.rindex('.')]+'_disjoin.json'
    with open(save_lc, "w") as f:
        f.write(json.dumps(c_d, ensure_ascii=False, indent=4, separators=(',', ':')))
    # return you can delete this
    return c_d #<core1:[chain1,chain2]>,<>....

getcore('ORGANIC\model\Valueki_result.smi')


# m = Chem.MolFromSmiles('O=C(NCc1cc(OC)c(O)cc1)Cc1cocc1CC')
# core = MurckoScaffold.GetScaffoldForMol(m)
# m_core = [m, core]
# Draw.MolsToGridImage(m_core, subImgSize=(250, 250))
