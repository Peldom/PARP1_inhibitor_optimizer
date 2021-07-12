import re
import gzip

# coding=utf-8

"""
Vocabulary helper class
"""

import re
import numpy as np


class Vocabulary:
    """Stores the tokens and their conversion to one-hot vectors."""

    def __init__(self, tokens=None, starting_id=0):
        """
        Instantiates a Vocabulary instance.
        :param tokens: A list of tokens (str).
        :param starting_id: The value for the starting id.
        :return:
        """
        self._tokens = {}
        self._current_id = starting_id

        if tokens:
            for token, idx in tokens.items():
                self._add(token, idx)
                self._current_id = max(self._current_id, idx + 1)

    def __getitem__(self, token_or_id):
        """
        Retrieves the if the token is given or a token if the id is given.
        :param token_or_id: A token or an id.
        :return: An id if a token was given or a token if an id was given.
        """
        return self._tokens[token_or_id]

    def add(self, token):
        """
        Adds a token to the vocabulary.
        :param token: Token to add.
        :return: The id assigned to the token. If the token was already there,
                 the id of that token is returned instead.
        """
        if not isinstance(token, str):
            raise TypeError("Token is not a string")
        if token in self:
            return self[token]
        self._add(token, self._current_id)
        self._current_id += 1
        return self._current_id - 1

    def update(self, tokens):
        """
        Adds many tokens at once.
        :param tokens: A list of tokens.
        :return: The ids of the tokens added.
        """
        return [self.add(token) for token in tokens]

    def __delitem__(self, token_or_id):
        """
        Deletes a (token, id) tuple, given a token or an id.
        :param token_or_id: A token or an id.
        :return:
        """
        other_val = self._tokens[token_or_id]
        del self._tokens[other_val]
        del self._tokens[token_or_id]

    def __contains__(self, token_or_id):
        """
        Checks whether a token is contained in the vocabulary.
        :param a token or an id to check
        :return : True if it is contained, otherwise False.
        """
        return token_or_id in self._tokens

    def __eq__(self, other_vocabulary):
        """
        Compares two vocabularies.
        :param other_vocabulary: Other vocabulary to be checked.
        :return: True if they are the same.
        """
        return self._tokens == other_vocabulary._tokens  # pylint: disable=W0212

    def __len__(self):
        """
        Calculates the length (number of tokens) of the vocabulary.
        :return : The number of tokens.
        """
        return len(self._tokens) // 2

    def encode(self, tokens):
        """
        Encodes a list of tokens, encoding them in 1-hot encoded vectors.
        :param tokens: Tokens to encode.
        :return : An numpy array with the tokens encoded.
        """
        ohe_vect = np.zeros(len(tokens), dtype=np.float32)
        try:
            for i, token in enumerate(tokens):
                ohe_vect[i] = self._tokens[token]
        except KeyError:
            return None
        else:
            return ohe_vect

    def decode(self, ohe_vect):
        """
        Decodes a one-hot encoded vector matrix to a list of tokens.
        :param : A numpy array with some encoded tokens.
        :return : An unencoded version of the input array.
        """
        tokens = []
        for ohv in ohe_vect:
            tokens.append(self[ohv])
        return tokens

    def _add(self, token, idx):
        if idx not in self._tokens:
            self._tokens[token] = idx
            self._tokens[idx] = token
        else:
            raise ValueError("IDX already present in vocabulary")

    def tokens(self):
        """
        Returns the tokens from the vocabulary.
        :return: A list of tokens.
        """
        return [t for t in self._tokens if isinstance(t, str)]


class SMILESTokenizer:
    """Deals with the tokenization and untokenization of SMILES."""

    REGEXPS = {
        "brackets": re.compile(r"(\[[^\]]*\])"),
        "2_ring_nums": re.compile(r"(%\d{2})"),
        "brcl": re.compile(r"(Br|Cl)")
    }
    REGEXP_ORDER = ["brackets", "2_ring_nums", "brcl"]

    def tokenize(self, smiles, with_begin_and_end=True):
        """
        Tokenizes a SMILES string.
        :param smiles: A SMILES string.
        :param with_begin_and_end: Appends a begin token and prepends an end token.
        :return : A list with the tokenized version.
        """
        def split_by(smiles, regexps):
            if not regexps:
                return list(smiles)
            regexp = self.REGEXPS[regexps[0]]
            splitted = regexp.split(smiles)
            tokens = []
            for i, split in enumerate(splitted):
                if i % 2 == 0:
                    tokens += split_by(split, regexps[1:])
                else:
                    tokens.append(split)
            return tokens

        tokens = split_by(smiles, self.REGEXP_ORDER)
        if with_begin_and_end:
            tokens = ["^"] + tokens + ["$"]
        return tokens

    def untokenize(self, tokens):
        """
        Untokenizes a SMILES string.
        :param tokens: List of tokens.
        :return : A SMILES string.
        """
        smi = ""
        for token in tokens:
            if token == "$":
                break
            if token != "^":
                smi += token
        return smi


def create_vocabulary(smiles_list, tokenizer):
    """
    Creates a vocabulary for the SMILES syntax.
    :param smiles_list: A list with SMILES.
    :param tokenizer: Tokenizer to use.
    :return: A vocabulary instance with all the tokens in the smiles_list.
    """
    tokens = set()
    for smi in smiles_list:
        tokens.update(tokenizer.tokenize(smi, with_begin_and_end=False))

    vocabulary = Vocabulary()
    vocabulary.update(["^"] + sorted(tokens))
    return vocabulary
def read_csv_file(file_path, ignore_invalid=True, num=-1):
    """
    Reads a SMILES file.
    :param file_path: Path to a CSV file.
    :param ignore_invalid: Ignores invalid lines (empty lines)
    :param num: Parse up to num rows.
    :return: An iterator with the rows.
    """
    with open_file(file_path, "rt") as csv_file:
        for i, row in enumerate(csv_file):
            if i == num:
                break
            fields = row.rstrip().split("\t")
            if fields:
                yield fields
            elif not ignore_invalid:
                yield None
def read_smi_file(file_path, ignore_invalid=True, num=-1):
    """
    Reads a SMILES file.
    :param file_path: Path to a SMILES file.
    :param ignore_invalid: Ignores invalid lines (empty lines)
    :param num: Parse up to num SMILES.
    :return: A list with all the SMILES.
    """
    return map(lambda fields: fields[0], read_csv_file(file_path, ignore_invalid, num))

def open_file(path, mode="r", with_gzip=False):
    """
    Opens a file depending on whether it has or not gzip.
    :param path: Path where the file is located.
    :param mode: Mode to open the file.
    :param with_gzip: Open as a gzip file anyway.
    """
    open_func = open
    if path.endswith(".gz") or with_gzip:
        open_func = gzip.open
    return open_func(path, mode)
# if __name__ == '__main__':
#     input_smiles_path='Value_Ki.smi'
def create_vocab(input_smiles_path):
    smiles_list = read_smi_file(input_smiles_path)

    tokenizer = SMILESTokenizer()
    vocabulary = create_vocabulary(smiles_list, tokenizer=tokenizer)

    tokens = vocabulary.tokens()
    # print("Vocabulary contains %d tokens: %s", len(tokens), tokens)
    return tokens
    # input_smiles_path = 'Value_Ki.csv'
    # smiles_list = read_csv_file(input_smiles_path)
    #
    # tokenizer = SMILESTokenizer()
    # vocabulary = create_vocabulary(smiles_list, tokenizer=tokenizer)
    #
    # tokens = vocabulary.tokens()
    # print("Vocabulary contains %d tokens: %s", len(tokens), tokens)
# vo=create_vocab('Value_Ki.smi')


def build_vocab(smiles, pad_char='_', start_char='^'):
    i = 1
    char_dict, ord_dict = {start_char: 0}, {0: start_char}
    for smile in smiles:
        # for c in smile:
        if smile not in char_dict:
            char_dict[smile] = i
            ord_dict[i] = smile
            i += 1
    char_dict[pad_char], ord_dict[i] = i, pad_char
    return char_dict, ord_dict
# char,ord=build_vocab(vo)

