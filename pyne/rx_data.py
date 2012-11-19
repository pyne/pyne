import re
import numpy as np

class rx_data(object):
    def __init__(self, data):
        self.data = data
    
    # def canonicalize(mat_state, data_type, rx_name):
        # return mat_state + 1, data_type + 1, rx_name + 1
    # def canonical_to_ENDF(mat_state, data_type, rx_name):
        # return mat_state, data_type, rx_name
    # def canonical_to_ACE(mat_state, data_type, rx_name):
        # return mat_state, data_type, rx_name
        
    def get(self, mat_state, data_type, rx_name, file_type):
        # mat_state, data_type, rx_name = canonicalize(mat_state, data_type, rx_name)

        if file_type.lower() == 'endf':
            # mat_state, data_type, rx_name = canonical_to_ENDF(mat_state, data_type, rx_name)
            result = self.read_mfmt(mat_state, data_type, rx_name)

        elif file_type.lower() == 'ace':
            table_name, rx_name = canonical_to_ACE(mat_state, data_type, rx_name)
            result = self.find_table(table_name).reactions[rx_name]
            print 'This is an ACE file and I don\'t know how to deal with it yet.'
            return False

        else:
            print "I don't know what that format is."
            return False

        return result
    
    def write(self, filename, file_type_in, file_type_out):
        # just a placeholder at this point
        # f = open(filename, 'w')
        if file_type_in.lower() == 'endf':
            for mat in self.mats:
                # f.close()
                return filename
        else:
            return filename, file_type
        pass
