import re
import numpy as np

class rx_data(object):
    def __init__(self, data):
        self.data = data
    
    # def canonicalize_reaction_name(mat_state, data_type, rx_name):
        # return mat_state + 1, data_type + 1, rx_name + 1
        
    def get(self, mat_state, data_type, rx_name, file_type):
        # mat_state, data_type, rx_name = canonicalize_reaction_name(mat_state, data_type, rx_name)
        if file_type.lower() == 'endf':
            if (mat_state, data_type, rx_name) in self.mts:
                print "Reading Material %d, MF %d, MT %d" % (mat_state, data_type, rx_name)
                start = (self.mts[(mat_state, data_type, rx_name)][0] - 1) * 6
                stop = (self.mts[(mat_state, data_type, rx_name)][1] - 1)* 6
            else:
                print "Material %d, File %d, MT %d does not exist." % (mat_id, mf, mt)
                return False
        else:
            print "I don't know what that format is."
            return False
        return self.data.flat[start:stop]


class ENDF_Library(rx_data):
    def __init__(self, filename):
        self.fh = open(filename, 'r')
        self.mats = {}
        self.mts = {}
        
        data = self._read_data(filename)
        rx_data.__init__(self, data)
        
        # are there more files to read the headers of?
        self.more_files = True
        
        # tracks theoretical length of file, based on headers
        self.chars_til_now = 0
        
        # tracks how many lines would have been skipped when reading data
        self.offset = 1

        # read ALL the headers!
        while self.more_files:
            self._read_headers()
            
    def _read_data(self, filename):
        print 'Reading data ...'
        data = np.genfromtxt(filename, 
                             delimiter = 11, 
                             usecols = (0, 1, 2, 3, 4, 5), 
                             invalid_raise = False,
                             skip_header = 1,                                    
                             converters = {0: convert,
                                           1: convert, 
                                           2: convert,
                                           3: convert, 
                                           4: convert, 
                                           5: convert})
        return data
    
    def _read_headers(self):
        
        # skip the first line and get the material id
        self.fh.seek(self.chars_til_now+81)
        line = self.fh.readline()
        mat_id = line[67:70]
        
        print 'Reading headers for material ID %d ...' % int(mat_id)
        
        # parse header (all lines with 1451)
        comments = ''
        mf = 1
        stop = self.chars_til_now/81
        
        while re.search('1451 +\d{1,3}', line):
    
            # parse contents section
            if re.match(' +\d{1,2} +\d{1,3} +\d{1,4} +', line):
                
                # while reading data we skip a line at the beginning
                # of every material, so we need an offset
                if int(mat_id) not in self.mats:
                    self.offset -= 1
                    
                self.mats.update({int(mat_id):(self.chars_til_now / 81)})
                
                # accounting for skipped lines between MF's and MT's
                old_mf = mf
                mf, mt = int(line[31:33]), int(line[41:44])
                mt_length = int(line[50:55])
                if old_mf == mf:
                    start = stop + 1
                else:
                    start = stop + 2
                    
                stop = start + mt_length
                self.mts.update({(int(mat_id), mf, mt):(start+self.offset, stop+self.offset)})
                line = self.fh.readline()
            elif re.search('C O N T E N T S', line):
                line = self.fh.readline()
                continue
            # parse comments
            else:
                comments = comments + '\n' + line[0:66]
                line = self.fh.readline()

        # find where end of material is
        self.chars_til_now = (stop + 4)*81
        
        # jump to end of this material         
        self.fh.seek(self.chars_til_now)
        
        # are we at the end of the file?
        if self.fh.readline() == '':
            self.more_files = False
        
        # update materials list
        if mat_id != '':
            self.mats.update({int(mat_id):(self.chars_til_now / 81, comments)})
        
    def get(self, mat_id, mf, mt):
        return rx_data.get(self, mat_id, mf, mt, 'endf')

def convert(string):
    """
    This function converts a number listed on an ENDF tape into a float or int
    depending on whether an exponent is present.
    """
  
    if re.search('[^ 0-9+\-\.]', string):
        return None
    elif string[-2] in '+-':
        return float(string[:-2] + 'e' + string[-2:])
    else:
        return float(string)
        
        
test = ENDF_Library('tests/endftest_small.txt')
print test.get(419, 4, 2)
        
