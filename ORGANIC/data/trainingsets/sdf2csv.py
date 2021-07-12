from rdkit import Chem
from rdkit.Chem.ChemUtils import SDFToCSV


def sdf2csv(sdf_name, csv_name):
    '''
    Convert sdf to csv
    '''
    f = open(csv_name, 'w')
    suppl = Chem.SDMolSupplier(sdf_name)
    SDFToCSV.Convert(suppl, f)
    f.close()
    print('Done.')
sdf2csv('csv.sdf','sdf.csv')