import re
import numpy as np
import collections
# import pyne.endf as endf

class RxLib(object):
    def __init__(self, data):
        self.data = data
        self.RxMaterials = {}
        
    
    # def canonicalize(mat_state, data_type, rx_name):
        # return mat_state + 1, data_type + 1, rx_name + 1
    # def canonical_to_ENDF(mat_state, data_type, rx_name):
        # return mat_state, data_type, rx_name
    # def canonical_to_ACE(mat_state, data_type, rx_name):
        # return mat_state, data_type, rx_name

    def write(self, filename, file_type_out):
        # just a placeholder at this point
        if file_type_out.lower() == 'endf':
            self.write_to_file(filename)
            # return filename
        else:
            return filename, file_type_out
        pass

class RxMaterial(RxLib):
    def __init__(self, mat_id, header, data, contents):
        self.mat_id = mat_id
        self.header = header
        self.data = data
        self.contents = contents

class double_spin_dict(collections.MutableMapping):


    def __init__(self, double_spin_dict):
        self.dict = double_spin_dict
        
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

