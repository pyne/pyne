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
    exe_path_abs = os.path.join('./',exe_path_abs)
    #print(exe_path_abs)
    return exe_path_abs

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
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    report_file = inputdict_unchecked['report_file']
    new_out = inputdict_unchecked['rewrite_input_with_hinderance_factor']
    output_file = 'alphad.out'
    if(new_out == 1):    
        output_file = inputdict_unchecked['output_file'] #output file if report = yes
        print("yes to output file..")
    exe_path = path_to_exe('alphad')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    inp = input_file + '\n' + report_file + '\n'
    print inp
    if (new_out == 1):
        inp = inp + 'Y' + '\n' + output_file
    else:
        inp = inp + 'N' + '\n'
    print inp
    proc.stdin.write(inp)
    proc.communicate()[0]
    proc.stdin.close()



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
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_file + '\n')
    proc.communicate()[0]
    proc.stdin.close()

#def gabs(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
#    inputdict = {}
#    input_file = inputdict_unchecked['input_file']
#    dataset_file = inputdict_unchecked['dataset_file']
#    output_file = inputdict_unchecked['output_file'] #report file << CHANGE BACK TO REPORT..

    #add option to not get new dataset (currently new dataset is hardprogrammed to yes)

#    exe_path = path_to_exe('gabs')
#    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
#    proc.stdin.write(input_file + '\n' + output_file + '\n' + 'Y' + '\n' + dataset_file )
#    proc.communicate()[0]
#    proc.stdin.close()

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
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    inp = input_file + '\n' + report_file + '\n'
    if (new_out):
        inp = inp + 'Y' + '\n' + output_file + '\n'
    else:
        inp = inp + 'N' + '\n'
    if (supress_g):
        inp = inp + 'Y' + '\n'
    else:
        inp = inp + 'N' + '\n'
    if(supress_ic):
        inp = inp + 'Y' + '\n'
    else:
        inp = inp + 'N' + '\n' + `dcc_theory` + '\n'
    proc.stdin.write(inp)
    proc.communicate()[0]
    proc.stdin.close()

def bldhst(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    output_table_file = inputdict_unchecked['output_table_file']
    output_index_file = inputdict_unchecked['output_index_file'] #report file << CHANGE BACK TO REPORT..

    exe_path = path_to_exe('bldhst')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_table_file + '\n' + output_index_file )
    proc.communicate()[0]
    proc.stdin.close()

def hsicc(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    data_deck = inputdict_unchecked['data_deck']
    icc_index = inputdict_unchecked['icc_index']
    icc_table = inputdict_unchecked['icc_table']
    complete_report = inputdict_unchecked['complete_report']
    new_card_deck = inputdict_unchecked['new_card_deck']
    comparison_report = inputdict_unchecked['comparison_report']
    multipol_known = inputdict_unchecked['is_multipol_known'] #'Y or CR'

    exe_path = path_to_exe('hsicc')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(data_deck + '\n' + icc_index + '\n' + icc_table + '\n' + complete_report + '\n' + new_card_deck + '\n' + comparison_report + '\n' + multipol_known )
    proc.communicate()[0]
    proc.stdin.close()

def hsmrg(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    data_deck = inputdict_unchecked['data_deck']
    card_deck = inputdict_unchecked['card_deck']
    merged_data_deck = inputdict_unchecked['merged_data_deck']

    exe_path = path_to_exe('hsmrg')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(data_deck + '\n' + card_deck + '\n' + merged_data_deck )
    proc.communicate()[0]
    proc.stdin.close()

def seqhst(inputdict_unchecked):
    #NOTE: changed input file line length to 90 to support longer file paths
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['binary_table_input_file']
    output_file = inputdict_unchecked['sequential_output_file']

    exe_path = path_to_exe('seqhst')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_file )
    proc.communicate()[0]
    proc.stdin.close()

def logft(inputdict_unchecked):
    #NOTE: changed input file line length to 90 to support longer file paths
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    input_data_set = inputdict_unchecked['input_data_set']
    output_report = inputdict_unchecked['output_report']
    data_table = inputdict_unchecked['data_table']
    output_data_set = inputdict_unchecked['output_data_set']

    exe_path = path_to_exe('logft')
    inp = input_data_set + '\n' + output_report + '\n' + data_table + '\n' + output_data_set + '\n'
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    proc.communicate()[0]
    proc.stdin.close()

def pandora(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    print('Executable not yet linked')
    #@todo: get path to executable
    inputdict = {}
    input_data_set = inputdict_unchecked['input_data_set']
    exe_path = path_to_exe('pandora')
    inp = input_data_set + '\n' + '0' + '\n' + '0' + '\n' + '0' + '\n' + '0' + '\n' + '0' + '\n'
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    print proc.communicate()[0]
    proc.stdin.close()

def radd(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """
    #@todo: check dictionary
    inputdict = {}
    atomic_number = inputdict_unchecked['atomic_number']
    neutron_number = inputdict_unchecked['neutron_number']

    exe_path = path_to_exe('radd')
    inp = atomic_number + '\n' + neutron_number + '\n' + 'NO' + '\n'
    print exe_path
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    print proc.communicate()[0]
    proc.stdin.close()

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
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    output_report_file = inputdict_unchecked['output_report_file']
    mode_of_operation = inputdict_unchecked['mode_of_operation']
    assumed_dcc_theory = inputdict_unchecked['assumed_dcc_theory']
    
    exe_path = path_to_exe('ruler')
    ruler_output = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    inp = input_file + '\n' + output_report_file + '\n' + mode_of_operation + '\n' + assumed_dcc_theory
    print(inp)
    ruler_output.stdin.write(inp)
    ruler_output.communicate()[0]
    ruler_output.stdin.close()

