import os
import nose
from nose.plugins.skip import Skip, SkipTest
import pyne
#try:
from pyne import ensdf_processing
#except:
#  raise SkipTest

def test_delta():
    input_dict = {}
    #input_dict['input_file'] = os.path.join('ensdf_processing', 'delta.dat')
    input_dict['input_file'] = 'ensdf_processing/delta.dat'
    input_dict['output_file'] = 'ensdf_processing/delta.out'
    output_dict = ensdf_processing.delta(input_dict)

#if __name__ == "__main__":
#  nose.runmodule()
if __name__ == "__main__":
    a = test_delta()
