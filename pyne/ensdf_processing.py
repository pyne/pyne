'''This module accesses various ensdf processing tools'''

import sys, os, os.path, shutil, urllib, subprocess, urllib2
import numpy as np

import urllib

from warnings import warn
from pyne.utils import QAWarning

if sys.version_info[0] > 2:
    basestring = str

warn(__name__ + " is not yet QA compliant.", QAWarning)

def path_to_exe(exe_name ):
    exe_path_abs, dp = os.path.split(os.path.abspath(__file__))
    exe_path_abs = os.path.join(exe_path_abs, exe_name)
    exe_path_abs = os.path.join('./',exe_path_abs)
    #print(exe_path_abs)
    return exe_path_abs

def download_exe(exe_path, exe_url):
    response = urllib2.urlopen(exe_url)
    print 'opened url'
    prog = 0
    CHUNK = 32 * 1024
    f = open(exe_path, 'wb')
    print 'opened file'
    while True:
        chunk = response.read(CHUNK)
        prog = prog + (32)
        print 'read chunk'
        print prog
        if not chunk: break
        f.write(chunk)
    f.close()
    return True

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
    rewrite_hinderance = inputdict_unchecked['rewrite_input_with_hinderance_factor']
    output_file = 'alphad.out'
    if(rewrite_hinderance == 1):    
        output_file = inputdict_unchecked['output_file'] #output file if report = yes
    exe_path = path_to_exe('alphad')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    inp = input_file + '\n' + report_file + '\n' + 'Y' + '\n'
    if (rewrite_hinderance == 1):
        inp = inp + 'Y' + '\n' + output_file
    else:
        inp = inp + 'N' + '\n'
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
    
    exe_path = path_to_exe('bricc')
    '''
    try:
        os.remove(exe_path) # purge for testing
        print 'file purged'
    except:
        print 'file not purged'
    if os.path.isfile(exe_path):
        print 'previous executable still present...'
    else:
        print 'executable clean/gone'  

    bricc_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/BrIcc/Linux/BriccV23-Linux.tgz"
    downloaded = download_exe(exe_path, bricc_url)
    '''
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    output_file = inputdict_unchecked['output_file']

    exe_path = path_to_exe('delta')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_file + '\n')
    proc.communicate()[0]
    proc.stdin.close()


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
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_file + '\n' + 'Y' + '\n' + dataset_file)
    proc.communicate()[0]
    proc.stdin.close()

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
    #@todo: get path to executable
    inputdict = {}
    input_data_set = inputdict_unchecked['input_data_set']
    exe_path = path_to_exe('pandora')
    inp = input_data_set + '\n' + \
        '0' + '\n' + \
        '0' + '\n' + \
        '0' + '\n' + \
        '0' + '\n' + \
        '0' + '\n'
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    print proc.communicate()[0]
    proc.stdin.close()
    if os.path.isfile('pandora.err'):  
        if 'output_err' in inputdict_unchecked:
            shutil.copyfile('pandora.err', inputdict_unchecked['output_err'])
    if os.path.isfile('pandora.gam'):
        if 'output_gam' in inputdict_unchecked:
            shutil.copyfile('pandora.gam', inputdict_unchecked['output_gam'])
    if os.path.isfile('pandora.gle'):  
        if 'output_gle' in inputdict_unchecked:
            shutil.copyfile('pandora.gle', inputdict_unchecked['output_gle'])
    if os.path.isfile('pandora.lev'):  
        if 'output_lev' in inputdict_unchecked:
            shutil.copyfile('pandora.lev', inputdict_unchecked['output_lev'])
    if os.path.isfile('pandora.rad'):  
        if 'output_rad' in inputdict_unchecked:
            shutil.copyfile('pandora.rad', inputdict_unchecked['output_rad'])
    if os.path.isfile('pandora.rep'):  
        if 'output_rep' in inputdict_unchecked:
            shutil.copyfile('pandora.rep', inputdict_unchecked['output_rep'])
    if os.path.isfile('pandora.xrf'):  
        if 'output_xrf' in inputdict_unchecked:
            shutil.copyfile('pandora.xrf', inputdict_unchecked['output_xrf'])
    if os.path.isfile('pandora.out'):  
        if 'output_out' in inputdict_unchecked:
            shutil.copyfile('pandora.out', inputdict_unchecked['output_out'])

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
    output_file = inputdict_unchecked['output_file']

    exe_path = path_to_exe('radd')
    inp = atomic_number + '\n' + neutron_number + '\n' + 'NO' + '\n'
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    radd_output = proc.communicate()[0]
    proc.stdin.close()
    f = open(output_file, 'w')
    f.write(radd_output)
    f.close()

def radlist(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)
    """

    print('Executable not yet linked')

    exe_path = path_to_exe('radlist')
    #try:
    #    os.remove(exe_path) # purge for testing
    #radlist_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/radlst/unx/radlist"
    #print exe_path
    #downloaded = download_exe(exe_path, radllist_url)

    inputdict = {}
    output_rad_listing = inputdict_unchecked['output_radiation_listing']
    output_endf_like_file = inputdict_unchecked['output_endf_like_file']
    output_file_for_nudat = inputdict_unchecked['output_file_for_nudat']
    output_mird_listing = inputdict_unchecked['output_mird_listing']
    calculate_continua = inputdict_unchecked['calculate_continua']
    input_file = inputdict_unchecked['input_file']
    output_radlst_file = inputdict_unchecked['output_radlst_file']
    input_radlst_data_table = inputdict_unchecked['input_radlst_data_table']
    input_masses_data_table = inputdict_unchecked['input_masses_data_table']
    output_ensdf_file = inputdict_unchecked['output_ensdf_file']

    inp = output_rad_listing + '\n' + output_endf_like_file + '\n' + output_file_for_nudat +\
          output_mird_listing + '\n' + calculate_continua + '\n' + input_file +\
          output_radlst_file + '\n' + input_radlst_data_table + '\n' + input_masses_data_table +\
          output_ensdf_file
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(inp)
    radd_output = proc.communicate()[0]
    proc.stdin.close()

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
    ruler_output.stdin.write(inp)
    ruler_output.communicate()[0]
    ruler_output.stdin.close()

