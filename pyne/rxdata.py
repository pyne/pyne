import collections

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
from warnings import warn
from pyne.utils import QA_warn

QA_warn(__name__)


class RxLib(object):
    """RxLib is a parent type that implements an abstract representation of
    nuclear data. Eventually it will be able to represent ENDF, ACE, and other
    filetypes.
    """

    def __init__(self, data):
        self.data = data

    def write(self, filename, file_type_out):
        pass


class DoubleSpinDict(collectionsAbc.MutableMapping):
    """Sanitizes input, avoiding errors arising from half-integer values of spin.

    Parameters
    ----------
    spin_dict: a dictionary where the keys are (spi, [L], [j]) tuples.
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
        try:
            if len(key) == 1:
                return int(round(2.0 * key[0]))
            if len(key) == 2:
                return (int(round(2.0 * key[0])), key[1])
            if len(key) == 3:
                return (int(round(2.0 * key[0])), key[1], key[2])
        except TypeError:
            return key
