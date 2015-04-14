'''This module accesses various ensdf processing tools'''

import sys
from warnings import warn
from pyne.utils import QAWarning

import numpy as np

if sys.version_info[0] > 2:
    basestring = str

warn(__name__ + " is not yet QA compliant.", QAWarning)

def alphad(inputdict_unchecked):
    """
    This function calculates the alpha hinderance factors and theoretical half 
    lives for even even ground state transitions.

    Input Dictionary Required Key Pair Value:
    @TODO: put in pretty table format
        ensdf_input_file : input file
        output_file : file for output to be written to (doesn't have to exist)

    Full documentation explaining the details of the functionality and physics
    behind ALPHAD can be found at:
        http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/alphad/readme-alphad.pdf
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out

def bricc(inputdict_unchecked):
    """
    This function calculates the conversion electron, electron-positron pair conversion
    coefficients, and the E0 electron factors.

    Input Dictionary Required Key Pair Value:
        input_index_file : input index file
        input_icc_file : input icc file
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out
