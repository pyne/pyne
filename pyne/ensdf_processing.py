'''This module accesses various ensdf processing tools'''

import sys
import os
from warnings import warn
from pyne.utils import QAWarning
import subprocess
import numpy as np

if sys.version_info[0] > 2:
    basestring = str

warn(__name__ + " is not yet QA compliant.", QAWarning)

def path_to_exe(exe_name ):
    exe_path_abs, dp = os.path.split(os.path.abspath(__file__))
    exe_path_abs = os.path.join(exe_path_abs, exe_name)
    sesh_path = sys.argv[0]  
    print(sesh_path)
    print(exe_path_abs)
    exe_path = os.path.relpath(exe_path_abs,sys.argv[0])
    exe_path = os.path.join('./',exe_path)
    exe_path = exe_path.replace("../","",1)
    print(exe_path)
    return exe_path

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

def delta(inputdict_unchecked):
    """
    This function calculates the best values of mixing ratios based of its analysis of
    the angular correlation and conversion coefficient data.

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output_file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    output_file = inputdict_unchecked['output_file']

    exe_path = path_to_exe('delta')
    delta_output = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    delta_output.stdin.write(input_file + '\n' + output_file + '\n')
    delta_output.communicate()[0]
    delta_output.stdin.close()

def gabs(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    dataset_file = inputdict_unchecked['dataset_file']
    output_file = inputdict_unchecked['output_file'] #report file << CHANGE BACK TO REPORT..

    #add option to not get new dataset (currently new dataset is hardprogrammed to yes)

    exe_path = path_to_exe('gabs')
    delta_output = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    delta_output.stdin.write(input_file + '\n' + output_file + '\n' + 'Y' + '\n' + dataset_file )
    delta_output.communicate()[0]
    delta_output.stdin.close()

def gtol(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    report_file = inputdict_unchecked['report_file']
    new_out = inputdict_unchecked['new_ensdf_file_with_results']
    output_file = inputdict_unchecked['output_file'] #output file if report = yes
    supress_g = inputdict_unchecked['supress_gamma_comparison']
    supress_ic = inputdict_unchecked['supress_intensity_comparison']
    dcc_theory = inputdict_unchecked['dcc_theory_percent']

    #add option to not get new dataset (currently new dataset is hardprogrammed to yes)

    exe_path = path_to_exe('gtol')
    delta_output = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    inp = input_file + '\n' + report_file + '\n'
    if (new_out):
        inp = inp + 'Y' + '\n' + output_file + '\n'
    if (supress_g):
        inp = inp + 'Y' + '\n'
    else:
        inp = inp + 'N' + '\n'
    if(supress_ic):
        inp = inp + 'Y' + '\n'
    else:
        inp = inp + 'N' + '\n' + `dcc_theory` + '\n'

    delta_output.stdin.write(inp)
    delta_output.communicate()[0]
    delta_output.stdin.close()

def hsicc(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out

def logft(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out

def pandora(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out

def radlist(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out

def ruler(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    #       call executable
    #       copy output file to specified out
