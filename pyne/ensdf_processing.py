'''This module accesses various ensdf processing tools'''

import sys, os, shutil, subprocess, tarfile
import numpy as np
from warnings import warn
from pyne.utils import QAWarning

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

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
            #if dl_size != 0:
            #    print 'Download progress: %d/100' % (100.0 * (float(prog) / float(dl_size)))
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
        ensdf_input_file : input file
        output_file : file for output to be written to (doesn't have to exist)

    Output Dictionary Values:
        Everything in input dictionary is returned if ALPHAD completes successfully.

    Full documentation explaining the details of the functionality and physics
    behind ALPHAD can be found at:
        http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/alphad/readme-alphad.pdf
    """
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
    
    # check if BriIccHome environment variable has been set (needed by BRICC executable)
    if not os.environ.get('BrIccHome'):
        os.environ['BrIccHome'] = str(exe_dir)
    
    input_line = inputdict_unchecked['input_line']
    proc = subprocess.Popen([exe_path],stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    proc.stdin.write(input_line + '\n' + 'exit' + '\n')
    bricc_output = proc.communicate()[0]
    proc.stdin.close()

    output_dict = inputdict_unchecked
    output_dict['output_file_directory'] = exe_dir
    output_dict['bricc_output'] = bricc_output
    return output_dict
