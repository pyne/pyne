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
    input_dict['input_file'] = 'ensdf_processing/delta.dat'
    input_dict['output_file'] = 'ensdf_processing/delta.out'
    output_dict = ensdf_processing.delta(input_dict)

def test_gabs():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gabs.in'
    input_dict['output_file'] = 'ensdf_processing/gabs.rpk'
    input_dict['dataset_file'] = 'ensdf_processing/gabs.dts'
    output_dict = ensdf_processing.gabs(input_dict)

def test_gtol():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gtol.inp'
    input_dict['report_file'] = 'ensdf_processing/gtol.rpt'
    input_dict['new_ensdf_file_with_results'] = 1
    input_dict['output_file'] = 'ensdf_processing/gtol.out'
    input_dict['supress_gamma_comparison'] = 0
    input_dict['supress_intensity_comparison'] = 0
    input_dict['dcc_theory_percent'] = 3
    output_dict = ensdf_processing.gtol(input_dict)

#if __name__ == "__main__":
#  nose.runmodule()
if __name__ == "__main__":
    a = test_delta()
    b = test_gabs()
    c = test_gtol()
