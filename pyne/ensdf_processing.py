'''This module accesses various ensdf processing tools'''

import sys, os, os.path, shutil, subprocess, tarfile
import numpy as np

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

from warnings import warn
from pyne.utils import QAWarning

if sys.version_info[0] > 2:
    basestring = str

warn(__name__ + " is not yet QA compliant.", QAWarning)

def path_to_exe(exe_name ):
    exe_path_abs, dp = os.path.split(os.path.abspath(__file__))
    exe_path_abs = os.path.join(exe_path_abs, exe_name)
    exe_path_abs = os.path.join('./',exe_path_abs)
    return exe_path_abs

def verify_download_exe(exe_path, exe_url, compressed = 0, decomp_path = '', dl_size = 0):
    if not os.path.exists(exe_path):
        print 'fetching executable'

        response = urllib2.urlopen(exe_url)
        prog = 0
        CHUNK = 32 * 1024
        f = open(exe_path, 'wb')
        while True:
            chunk = response.read(CHUNK)
            prog = prog + (256)
            if dl_size != 0:
                print 'Download progress: %d/100' % (100.0 * (float(prog) / float(dl_size)))
            if not chunk: break
            f.write(chunk)
        f.close()
        # set proper permissions on newly downloaded file
        os.chmod(exe_path, 0744)
        if compressed:
            tfile = tarfile.open(exe_path, 'r:gz')
            tfile.extractall(decomp_path)

def alphad(inputdict_unchecked):
    """
    This function calculates the alpha hinderance factors and theoretical half 
    lives for even even ground state transitions.

    Input Dictionary Required Key Pair Value:
    @TODO: put in pretty table format
        ensdf_input_file : input file
        output_file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if ALPHAD completes successfully.

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
    return inputdict_unchecked

def bricc(inputdict_unchecked):
    """
    This function calculates the conversion electron, electron-positron pair conversion
    coefficients, and the E0 electron factors.

    Input Dictionary Required Key Pair Value:
        input_line : input line to be used, as if bricc was being used in the interactive
                     command line style.  This line will be passed in verbatim to the
                     bricc executable, and default output file names will be used.  For
                     generating new records, or any operation that produces output files,
                     '<CR>' should be appended on the end of input_line, to force default
                     file names to be used.
    Output Dictionary Values:
        input_line : user defined input
        output_file_directory : the directory all produced bricc output files will be 
                                located.
        bricc_output : data printed to command line.  Useful for interactive use.
    NOTE:
        All the various ouptput files bricc can generate are found in the
        'output_file_directory' path.  '<CR>' must be appended to 'input_line' for this
        to work properly. 
    """
    
    exe_path = path_to_exe('bricc')
    exe_dir = path_to_exe('')
    compressed_exe_path = exe_path + '.tar.gz'

    bricc_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/BrIcc/Linux/BriccV23-Linux.tgz"
    decomp_exe_path = path_to_exe('')
    decomp_options = ['bricc', '.tgz', True]
    verify_download_exe(compressed_exe_path, bricc_url, compressed = True, decomp_path = decomp_exe_path, dl_size = 127232)
    #@todo: check dictionary
    
    # check if BriIccHome environment variable has been set (needed by BRICC executable)
    if not os.environ.get('BrIccHome'):
        os.environ['BrIccHome'] = str(exe_dir)
    
    input_line = inputdict_unchecked['input_line']
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_line)
    casdfads = proc.communicate()[0]
    proc.stdin.close()

    output_dict = inputdict_unchecked
    output_dict['output_file_directory'] = exe_dir
    output_dict['bricc_output'] = casdfads
    return output_dict

    
def delta(inputdict_unchecked):
    """
    This function calculates the best values of mixing ratios based of its analysis of
    the angular correlation and conversion coefficient data.

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output_file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if DELTA completes successfully.
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
    return inputdict_unchecked

def gabs(inputdict_unchecked):
    """
    This function ...

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        dataset_file : dataset file to be used
        output file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if GABS completes successfully.
    """
    exe_path = path_to_exe('gabs') 

    gabs_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/gabs/unx/gabs"
    verify_download_exe(exe_path, gabs_url, dl_size = 8704)
    
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
        input_file : input ensdf file.
        report_file : desired gtol report file path.
        new_ensdf_file_with_results : boolean, if true then a new ensdf file with results
                                      will be created.
        output_file : desired gtol output file path.
        supress_gamma_comparison : boolean, if true the gamma comparison will be suppressed.
        dcc_theory_percent : double, specifies the dcc theory percentage to be used.


    Output Dictionary Values:
        Everything in input dictionary is returned if GTOL completes successfully.
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
    return inputdict_unchecked

def bldhst(inputdict_unchecked):
    """
    This program builds a direct access file of the internal conversion coefficient 
    table. (BLDHST readme)

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file.
        output_table_file : desired output table file path.
        output_index_file : desired output index file path.


    Output Dictionary Values:
        Everything in input dictionary is returned if BLDHST completes successfully.
    """
    #@todo: check dictionary
    inputdict = {}
    input_file = inputdict_unchecked['input_file']
    output_table_file = inputdict_unchecked['output_table_file']
    output_index_file = inputdict_unchecked['output_index_file']

    exe_path = path_to_exe('bldhst')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_table_file + '\n' + output_index_file )
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked

def hsicc(inputdict_unchecked):
    """
    This program calculates internal conversion coefficients. (HSICC readme)

    Input Dictionary Required Key Pair Value:
        data_deck : data deck to be used for hsicc program.
        icc_index : icc index to be used for hsicc program.
        icc_table : icc table to be used for the hsicc program.
        complete_report : desired report file path for hsicc program.
        new_card_deck : desired new card deck file path for hsicc program.
        comparison_report : desired comparison report path for hsicc program.
        is_multipol_known : 1 if multipol is known, 0 otherwise.

    Output Dictionary Values:
        Everything in input dictionary is returned if HSICC completes successfully.
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
    proc.stdin.write(data_deck + '\n' + icc_index + '\n' + icc_table + '\n' + \
        complete_report + '\n' + new_card_deck + '\n' + comparison_report + '\n' + multipol_known )
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked

def hsmrg(inputdict_unchecked):
    """
    This program merges new gamma records created by HSICC with the original input 
    data.  (HSICC readme)

    Input Dictionary Required Key Pair Value:
        data_deck : data deck file path for hsmrg to use.
        card_deck : card deck file path for hsmrg to use.
        merged_data_deck : desired merged data deck file path created by hsmrg.

    Output Dictionary Values:
        Everything in input dictionary is returned if HSMRG completes successfully.
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
    return inputdict_unchecked

def seqhst(inputdict_unchecked):
    """
    This program recreates a sequential file of the internal conversion table from the 
    direct access file.  (HSICC readme)

    Input Dictionary Required Key Pair Value:
        binary_table_input_file : binary table input file path.
        sequential_output_file : desired path of sequential output file.

    Output Dictionary Values:
        Everything in input dictionary is returned if SEQHST completes successfully.
    """
    #NOTE: changed input file line length to 90 to support longer file paths
    inputdict = {}
    input_file = inputdict_unchecked['binary_table_input_file']
    output_file = inputdict_unchecked['sequential_output_file']

    exe_path = path_to_exe('seqhst')
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_file + '\n' + output_file )
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked

def logft(inputdict_unchecked):
    #NOTE: changed input file line length to 90 to support longer file paths
    """
    This program calculates log ft values for beta and electron-capture decay, average beta energies, 
    and capture fractions.  (LOGFT readme)

    Input Dictionary Required Key Pair Value:
        input_data_set : path to input data file.
        output_report : desired path to output report file.
        data_table : path to data table.
        output_data_set : desired path to output data set.

    Output Dictionary Values:
        Everything in input dictionary is returned if LOGFT completes successfully.
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
    return inputdict_unchecked

def pandora(inputdict_unchecked):
    """
    This program performs a least-squares fit to the gamma-energies to obtain level energies 
    and calculates the net feeding to levels. (PANDORA readme)

    Input Dictionary Required Key Pair Value:
        xxxxxxx

    Output Dictionary Values:
        Everything in input dictionary is returned if PANDORA completes successfully.
    """
    #@todo: get path to executable
    inputdict = {}
    input_data_set = inputdict_unchecked['input_data_set']
    exe_path = path_to_exe('pandora')

    # create temp file for fortran file to write temporary things to
    tmpfile_path = path_to_exe('tmpfil.tmp')
    tmpfile = open(tmpfile_path, 'w+')
    open('tmpfil.tmp', 'w+')
    tmpfile.close()
    inp = input_data_set + '\n' + \
         ('1' if inputdict_unchecked['level_report_and_files_sorted'] else '0') + '\n' + \
         ('1' if inputdict_unchecked['gamma_report_and_files_sorted'] else '0') + '\n' + \
         ('1' if inputdict_unchecked['radiation_report_and_files_sorted'] else '0') + '\n' + \
         ('1' if inputdict_unchecked['cross_reference_output'] else '0') + '\n'
    if inputdict_unchecked['cross_reference_output']:
        inp = inp + 'pandora.out' + '\n'
    inp = inp + '1' if  inputdict_unchecked['supress_warning_messages'] else '0 + \n'
    print inp
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
        else:
            print 'could not find output in dictionary'
    else:
        print 'could not find output file'
    return inputdict_unchecked



def radd(inputdict_unchecked):
    """
    This code (RadD.FOR) deduces the radius parameter (r 0 ) for odd-odd and odd-A nuclei 
    using the even-even radii [1] as input parameters. These radii deduced for odd-A and 
    odd-odd nuclides can be used in the calculation of alpha hindrance factors. In this 
    procedure, it is assumed that radius parameter ( r 0 Z , N ) for odd-Z and odd-N 
    nuclides lies midway between the radius parameters of adjacent even-even neighbors 
    calculates reduced transition probabilities. (RADD readme)

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if RADD completes successfully.
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
    return inputdict_unchecked

def radlist(inputdict_unchecked):
    """
    This program calculates atomic & nuclear radiations and checks energy balance.  
    (RADLIST readme)
    
    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if RADLIST completes successfully.
    """
    exe_path = path_to_exe('radlist')
    radlist_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/radlst/unx/radlist"
    print exe_path
    verify_download_exe(exe_path, radlist_url, dl_size = 8704)
    
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
    return inputdict_unchecked
    

def ruler(inputdict_unchecked):
    """
    This program calculates reduced transition probabilities. (RULER readme)

    Input Dictionary Required Key Pair Value:
        input_file : input ensdf file
        output file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if RULER completes successfully.
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
    return inputdict_unchecked
