import rdkit
from rdkit import Chem
from rdkit.Chem import QED
import numpy
import pandas as pd
import matplotlib.pyplot as plt
import SAscoretest

def maxminnorm(x):
    # 使用均值方差归一化
    mu = numpy.average(x)
    sigma = numpy.std(x)
    x = (x - mu) / sigma
    return x

def countQED(x):
    #输入你的文件的绝对路径，这次的文件是smiles格式。
    suppl=Chem.SmilesMolSupplier(x)
    QEDds=[]
    for mol in suppl:
        QEDds.append( (Chem.QED.qed(mol)))
    print('Length=',len(QEDds),'Mean QED=',numpy.mean(QEDds))
    return numpy.array(QEDds)


def show_plt(inputdata):
    n, bins, patches = plt.hist(x=inputdata, bins='auto', color='#FF1493',alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('QED distribution')
    plt.text(23, 45, r'$\mu=15, b=3$')
    maxfreq = n.max()
    # plt.savefig('test.png')
    plt.show()
    return

if __name__ == '__main__':
    ds_path = 'ORGANIC\model\Valueki_result.smi'
    QEDarray = maxminnorm(countQED(ds_path))
    sascorearray = maxminnorm(SAscoretest.count_SA_score(ds_path))
    sum_score = QEDarray-2*sascorearray
    b = sum_score.tolist()
    print(sum_score)
    print(b.index(max(b)))
    print(b.index(min(b)))

# show_plt(countQED('ORGANIC\experiments\Value_Ki.smi'))
