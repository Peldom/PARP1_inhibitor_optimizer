import rdkit
from rdkit import Chem
from rdkit.Chem import RDConfig
import sys,os,numpy
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer as sa
import matplotlib.pyplot as plt


def count_SA_score(x):
    suppl=Chem.SmilesMolSupplier(x)
    SAscores=[]
    for mol in suppl:
        SAscores.append(sa.calculateScore((mol)))
    print('Length=',len(SAscores),'Mean SA_score=',numpy.mean(SAscores))
    return numpy.array(SAscores)

def show_plt(inputdata):
    n, bins, patches = plt.hist(x=inputdata, bins='auto', color='#1E90FF',alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('SA_Score distribution')
    plt.text(23, 45, r'$\mu=15, b=3$')
    maxfreq = n.max()
    # plt.savefig('test.png')
    plt.show()
    return

if __name__ == '__main__':
    show_plt(count_SA_score('ORGANIC\experiments\Value_Ki.smi'))