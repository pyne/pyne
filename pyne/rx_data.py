import re
import numpy as np
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
