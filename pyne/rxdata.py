import re
import numpy as np
import collections
# import pyne.endf as endf

class RxLib(object):
    """RxLib is a parent type which implements an abstract representation of
    nuclear data. Eventually it will be able to represent ENDF, ACE, and other
    filetypes.
    """
    def __init__(self, data):
        self.data = data

    def write(self, filename, file_type_out):
        pass

class DoubleSpinDict(collections.MutableMapping):
    """DoubleSpinDict is a dictionary that takes in half-spin numbers and represents
    them as integers. This avoids floating point problems causing key errors in the 
    dictionary.

    Parameters
    ----------
    spin_dict: a dictionary where the keys are (spi, L, j) tuples
        No need to pre-double spin to get an integer.
    """
    def __init__(self, spin_dict):
        self.dict = spin_dict

    def __len__(self):
        return len(self.dict)

    def __iter__(self):
        return self.dict.iterkeys()

    def __contains__(self, key):
        return self.double_spin(key) in self.dict

    def __getitem__(self, key):
        return self.dict.get(self.double_spin(key))

    def __setitem__(self, key, value):
        self.dict[self.double_spin(key)] = value

    def __delitem__(self, key):
        del self.dict[self.double_spin(key)]

    def double_spin(self, key):
        return (int(round(2.0*key[0])), key[1], key[2])
