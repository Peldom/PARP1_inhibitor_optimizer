## Other
import copy
import itertools

from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog
from tqdm import tqdm
import os
import sys

import pickle
import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse



## Rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem import QED
from rdkit.Chem import Draw
import rdkit.Chem.rdchem as rdchem

# Setup SA scorer
from rdkit.Chem import RDConfig

from model import check_valency
from collections import namedtuple

MCNNInput = namedtuple('MCNNInput', ['node_type', 'adj_type', 'graph_size', 'properties', 'smile', 'validity_masks'])


sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


## Pytorch
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset


def calculateSA(m):
    sascore = sascorer.calculateScore(m)
    return sascore

def cyclescore(mol):
    cycle_list = nx.cycle_basis(nx.Graph(Chem.rdmolops.GetAdjacencyMatrix(mol, useBO=True)))
    if len(cycle_list) == 0:
        cycle_length = 0
    else:
        cycle_length = max([len(j) for j in cycle_list])
    if cycle_length <= 6:
        cycle_length = 0
    else:
        cycle_length = cycle_length - 6
    cycle_score = -cycle_length
    return cycle_score

def penalized_logp(mol):
    """
    Reward that consists of log p penalized by SA and # long cycles,
    as described in (Kusner et al. 2017). Scores are normalized based on the
    statistics of 250k_rndm_zinc_drugs_clean.smi dataset
    :param mol: rdkit mol object
    :return: float
    """
    # normalization constants, statistics from 250k_rndm_zinc_drugs_clean.smi
    logP_mean = 2.4570953396190123
    logP_std = 1.434324401111988
    SA_mean = -3.0525811293166134
    SA_std = 0.8335207024513095
    cycle_mean = -0.0485696876403053
    cycle_std = 0.2860212110245455

    log_p = Descriptors.MolLogP(mol)
    SA = -calculateSA(mol)

    # cycle score
    cycle_list = nx.cycle_basis(nx.Graph(
        Chem.rdmolops.GetAdjacencyMatrix(mol)))
    if len(cycle_list) == 0:
        cycle_length = 0
    else:
        cycle_length = max([len(j) for j in cycle_list])
    if cycle_length <= 6:
        cycle_length = 0
    else:
        cycle_length = cycle_length - 6
    cycle_score = -cycle_length

    normalized_log_p = (log_p - logP_mean) / logP_std
    normalized_SA = (SA - SA_mean) / SA_std
    normalized_cycle = (cycle_score - cycle_mean) / cycle_std

    return normalized_log_p + normalized_SA + normalized_cycle

def split_dataset(dataset, val_prop=.1, test_prop=.1, seed=42):
    assert val_prop + test_prop < 1.0
    test_seed = 42  # Do not change.

    full_size = len(dataset)
    val_size = int(val_prop * full_size)
    test_size = int(test_prop * full_size)
    train_size = full_size - val_size - test_size

    # We want the test set to stay fixed even if we change the seed.
    torch.manual_seed(test_seed)
    train_val, test = torch.utils.data.random_split(dataset, [train_size + val_size, test_size])
    torch.manual_seed(seed)
    train, val = torch.utils.data.random_split(train_val, [train_size, val_size])
    return train, val, test


def rare_event(mol):
    rare_event_cnt = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.TRIPLE:
            rare_event_cnt += 1

    for atom in mol.GetAtoms():
        if atom.GetSymbol() in {'Br', 'Cl', 'F'}:
            rare_event_cnt += 1

    return rare_event_cnt

class SmilesDataset(Dataset):
    file_smi = '.smi'
    file_pickle = '.pickle'
    file_bfs = '_bfs'
    END_GRAPH = None  # The end-graph token. (Can be anything that evaluates to False)
    bond_types = (0, rdchem.BondType.SINGLE, rdchem.BondType.DOUBLE, rdchem.BondType.TRIPLE)  # 0 for NO_BOND.
    properties = {'hmolw': Descriptors.HeavyAtomMolWt,  # FIXME: Is not the same as the statistics!
                  'molw': Descriptors.MolWt,
                  'logp': Descriptors.MolLogP,
                  'plogp': penalized_logp,
                  'SA': calculateSA,
                  'cycle': cyclescore,
                  'qed': QED.qed,
                  'rare': rare_event}

    def __init__(self, smi_path, validity_masks=True, conditioning=None):
        if conditioning is None:
            conditioning = []
        self.validity_masks = validity_masks
        self.conditioning = sorted(conditioning)
        # It might be faster if we create validity masks as preprocessing, but it will probably take too much memory.

        # Change .smi to .pickle, and look up to see if the preprocessed dataset exists.
        pickle_path = smi_path[:-len(self.file_smi)] + self.file_pickle
        if not os.path.exists(pickle_path) or os.path.getmtime(pickle_path) < os.path.getmtime(smi_path):
            try:
                os.remove(pickle_path)
                print('SMILES file newer than pickle file - processing SMILES file to create an updated pickle file.')
            except FileNotFoundError:
                print('Found no pickle file - processing SMILES file to create pickle file.')
            self.preprocess(smi_path, pickle_path)

        print('Loading dataset', pickle_path)
        with open(pickle_path, "rb") as f:
            self.__dict__.update(pickle.load(f))

        # We store the lists of atom/bond types in _types. Then self.atom_types[i] = self.atoms[self.atom_types[i]].
        self.atom_types = [SmilesDataset.END_GRAPH] + self.atoms  # self.atoms is read from pickle file.
        self.atoms = {atom: idx for idx, atom in enumerate(self.atom_types)}
        self.bonds = {bond: idx for idx, bond in enumerate(SmilesDataset.bond_types)}  # Overrides class variable.
        self.max_num_atoms = self.max_num_atoms + 2  # add 1 for end-graph node. 1 more for some mask stuff.

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        smile = self.smiles[idx]
        mol = Chem.MolFromSmiles(smile, sanitize=False)
        mol = self.process_mol(mol)

        ### NETWORK INPUTS/TARGETS
        node_type, adj_type = self.encode_mol(mol)
        node_type = F.pad(node_type, [0, 1], value=self.atoms[self.END_GRAPH]) # Add END_GRAPH token to nodes.

        # Apply some padding, so tensors are of the same size.
        padded_node_type = node_type.new_zeros((self.max_num_atoms,))
        padded_node_type[:len(node_type)] = node_type
        padded_adj_type = adj_type.new_zeros((self.max_num_atoms, self.max_num_atoms))
        padded_adj_type[:len(adj_type), :len(adj_type)] = adj_type
        graph_size = torch.tensor(len(node_type), dtype=torch.long)  # Number of atoms + 1.

        ### PROPERTIES
        properties = self.get_properties(mol, self.conditioning)

        ### VALIDITY MASKS
        validity_mask = None
        if self.validity_masks:
            validity_mask = self.get_validity_masks(node_type, adj_type)

        return MCNNInput(padded_node_type, padded_adj_type, graph_size, properties, smile, validity_mask)

    @staticmethod
    def get_properties(mol, conditioning):
        properties = []
        for prop_key in conditioning:
            # Grab mappings from molecule to a property, and collect result into a property vector values.
            properties.append(SmilesDataset.properties[prop_key](mol))
        return torch.tensor(properties, dtype=torch.float)

    def get_validity_masks(self, node_type, adj_type, return_mol=False):
        max_node_unroll = self.max_num_atoms
        max_edge_unroll = self.max_graph_bandwidth
        num_masks = (max_node_unroll  # Number of nodes.
                     + ((max_edge_unroll - 1) * max_edge_unroll) // 2
                     + (max_node_unroll - max_edge_unroll) * max_edge_unroll)
        num_mask_edge = num_masks - max_node_unroll

        valid_edge_mask = torch.ones((num_mask_edge, len(self.bond_types)), dtype=torch.bool)
        edge_step = 0
        rw_mol = Chem.RWMol()
        # We have to enumerate the atoms in the right order.
        for i, atom_idx in enumerate(node_type):
            atom_type = self.atom_types[atom_idx]
            if atom_type:
                rw_mol.AddAtom(Chem.Atom(atom_type))
            else:
                break

            edge_total = min(i, self.max_graph_bandwidth)  # The number of edges to sample.
            start = max(0, i - self.max_graph_bandwidth)
            for j in range(start, start + edge_total):
                # Now to make the valid edge prediction mask.
                for b, edge_type in enumerate(self.bond_types[1:], start=1):
                    # If this is an actual bond. (NO_BOND is always valid)
                    rw_mol.AddBond(i, j, edge_type)
                    valid_edge_mask[edge_step, b] = check_valency(rw_mol)
                    rw_mol.RemoveBond(i, j)

                bond_type = self.bond_types[adj_type[i, j]]
                if bond_type:
                    rw_mol.AddBond(i, j, bond_type)
                edge_step += 1
        # Now we have recreated the molecule.

        assert check_valency(rw_mol)  # Gave issues with Zinc250, which turned out to contain ions/charged atoms.
        # This meant our model actually learned to generate invalid molecules.

        if return_mol: # For plotting.
            return rw_mol
        return valid_edge_mask

    @staticmethod
    def reorder(adj_type):
        # Reverse Cuthill-McKee
        csc = scipy.sparse.csc_matrix(adj_type)
        # order = np.ascontiguousarray(scipy.sparse.csgraph.reverse_cuthill_mckee(csc, symmetric_mode=True)[::-1])
        order = np.ascontiguousarray(scipy.sparse.csgraph.reverse_cuthill_mckee(csc)[::-1])
        return order

    def encode_mol(self, mol, reorder=True):
        # Make atom type matrix
        node_type = torch.as_tensor([self.atoms[atom.GetSymbol()] for atom in mol.GetAtoms()], dtype=torch.long)

        # Adjacency matrix (we assume all bond types are integer, i.e. no aromatic 1.5 bond.)
        adj_type = torch.as_tensor(Chem.GetAdjacencyMatrix(mol, useBO=True), dtype=torch.long)

        if not reorder:
            return node_type, adj_type

        # Performs BFS re-ordering to minimize adjacency matrix bandwidth.
        order = self.reorder(adj_type)

        return node_type[order], adj_type[np.ix_(order, order)]

    def preprocess(self, smi_path, pickle_path):
        smiles = pd.read_csv(smi_path, sep='\n', header=None)[0].values
        all_atoms = set()

        # Some statistics to capture.
        num_atoms = []
        graph_bandwidths = []  # The bandwidth of the adjacency tensors after reordering.
        weight = []
        logP = []
        SA = []
        cycle =  []
        all_smiles = []
        unique_smiles = set()
        # We could use multiprocessing to speed this up. But it takes only around 1 hour.
        for _, smile in enumerate(tqdm(smiles, desc='Processing ' + smi_path)):
            mol = self.process_mol(Chem.MolFromSmiles(smile))
            ### Find unique atoms
            atoms = mol.GetAtoms()

            # Keep track of all the unique atoms seen.
            all_atoms.update(atom.GetSymbol() for atom in atoms)

            # Calculate graph bandwidth. (Used for bond prediction)
            adjmat = Chem.GetAdjacencyMatrix(mol, useBO=True)
            order = self.reorder(adjmat)
            foo, bar = adjmat[np.ix_(order, order)].nonzero()  # rows and cols where non-zeros are located
            graph_bandwidths.append(int(max(foo - bar)) + 1)  # add 1 for correctness

            # Other statistics
            weight.append(Descriptors.HeavyAtomMolWt(mol))
            logP.append(Descriptors.MolLogP(mol))
            SA.append(calculateSA(mol))
            cycle.append(cyclescore(mol))
            num_atoms.append(len(atoms))

            # Convert processed molecule back to SMILES.
            # mols.append(mol)
            smile = Chem.MolToSmiles(mol, isomericSmiles=False)
            all_smiles.append(smile)  # It is saved as a Canonical SMILES.
            unique_smiles.add(smile)

        pt = Chem.GetPeriodicTable()
        processed_data = {'smiles': all_smiles,
                          'atoms': sorted(list(all_atoms), key = lambda atom: pt.GetAtomicNumber(atom)),
                          'max_num_atoms': max(num_atoms),
                          'max_graph_bandwidth': max(graph_bandwidths),
                          'mean_num_atoms': np.mean(num_atoms), 'std_num_atoms': np.std(num_atoms),
                          'mean_bandwidth': np.mean(graph_bandwidths), 'std_bandwidth': np.std(graph_bandwidths),
                          'mean_logP': np.mean(logP), 'std_logP': np.std(logP),
                          'mean_hMolWt': np.mean(weight), 'std_hMolWt': np.std(weight),
                          'mean_SA': np.mean(SA), 'std_SA': np.std(SA),
                          'mean_cycle': np.mean(cycle), 'std_cycle': np.std(cycle)}

        with open(pickle_path, "wb") as f:
            pickle.dump(processed_data, f)

        with open(pickle_path[:-len(self.file_pickle)] + '_unique' + self.file_pickle, "wb") as f:
            pickle.dump(unique_smiles, f)

    @staticmethod
    def process_mol(mol):
        ### convert smiles to the molecular form we want
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol)
        Chem.RemoveStereochemistry(mol)  # (-In case it's present) Our model doesn't support Stereochemistry.
        mol = Chem.RemoveHs(mol, sanitize=False)

        # From https://www.rdkit.org/docs/GettingStartedInPython.html:
        # "Note: as of this writing (Aug 2008), the smiles provided when one requests kekuleSmiles are not canonical.
        # The limitation is not in the SMILES generation, but in the kekulization itself."
        return mol

    @staticmethod
    def bfs_seq(G, start_id):
        '''
        get a bfs node sequence
        :param G:
        :param start_id:
        :return:
        '''
        dictionary = dict(nx.bfs_successors(G, start_id))
        start = [start_id]
        output = [start_id]
        while len(start) > 0:
            next = []
            while len(start) > 0:
                current = start.pop(0)
                neighbor = dictionary.get(current)
                if neighbor is not None:
                    next = next + neighbor
            output = output + next
            start = next
        return output

def smile_to_matrix(smile):
    # Converts smile to image matrix
    m = Chem.MolFromSmiles(smile)
    im = Draw.MolToImage(m)
    M = np.array(im)  # OBS - this is saved as a RGBA!
    return M


if __name__ == "__main__":

    #smi_path = 'dataset/250k_rndm_zinc_drugs_clean_sorted.smi'

    ### Process dataset
    # DataSet(smi_path,bfs=True)
    #ds = SmilesDataset(smi_path, bfs=False)
    #ds[0]

    triple_bonds(Chem.MolFromSmiles('C#C'))