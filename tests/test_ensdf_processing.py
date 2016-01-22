import filecmp, numpy, os, shutil
from pyne import ensdf_processing

import nose

# path to folder for temporary test files.
tmp_path = 'ensdf_processing/tmp'

def test_alphad():
    create_tmp()
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/alphad/ref_a228.ens'
    input_dict['report_file'] = tmp_path + '/tmp_alphad.rpt'
    input_dict['rewrite_input_with_hinderance_factor'] = 1
    input_dict['output_file'] = tmp_path + '/tmp_alphad.out'
    output_dict = ensdf_processing.alphad(input_dict)
    exceptions = [[2, 'DATE RUN']]
    file_comp(input_dict['report_file'],'ensdf_processing/alphad/ref_a228.ens.alphad.rpt', exceptions)
    cleanup_tmp()

# Tests gabs output for 80Br sample input.  Date and file names are expected to be different
# in the test output and reference file sets.
def test_gabs():
    create_tmp()
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gabs/ref_gabs_80Br.in'
    input_dict['output_file'] = tmp_path + '/tmp_gabs_80Br.rpt'
    input_dict['dataset_file'] = tmp_path + '/tmp_gabs_80Br.new'
    output_dict = ensdf_processing.gabs(input_dict)
    exceptions_output = [[4,0],[1, '  * * * GABS Version 11 '], 
                         [1, '        Current date: '],
                         [1, '        ENSDF input file: '],
                         [1, '        new ENSDF file:']]
    exceptions_dataset = [[4,0]]
    d_report1 = file_comp(input_dict['output_file'],
                          'ensdf_processing/gabs/ref_gabs_80Br.rpt', exceptions_output)
    d_report2 = file_comp(input_dict['dataset_file'],
                          'ensdf_processing/gabs/ref_gabs_80Br.new', exceptions_dataset)
    cleanup_tmp()


def test_logft():
    create_tmp()
    input_dict = {}
    input_dict['input_data_set'] = 'ensdf_processing/logft/ref_logft.inp'
    input_dict['output_report'] = tmp_path + '/tmp_logft.rpt'
    input_dict['data_table'] = 'ensdf_processing/logft/ref_logft.dat'
    input_dict['output_data_set'] = tmp_path + '/tmp_logft.new'
    output_dict = ensdf_processing.logft(input_dict)
    ref_output_data_set = 'ensdf_processing/logft/ref_logft.new'
    d_data = file_comp(input_dict['output_data_set'], ref_output_data_set, [])
    cleanup_tmp()

def test_radd():
    create_tmp()
    input_dict = {}
    input_dict['atomic_number'] = '86'
    input_dict['neutron_number'] = '113'
    input_dict['output_file'] = tmp_path + '/tmp_output.out'
    ensdf_processing.radd(input_dict)
    ref_output = 'ensdf_processing/radd/ref_output.out'
    d_report = file_comp(input_dict['output_file'], ref_output, [])
    cleanup_tmp()

def test_radlist():
    create_tmp()
    input_dict = {}
    input_dict['output_radiation_listing'] = 'Y'
    input_dict['output_ensdf_like_file'] = 'N'
    input_dict['output_file_for_nudat'] = 'N'
    input_dict['output_mird_listing'] = 'N'
    input_dict['calculate_continua'] = 'N'
    input_dict['input_file'] = 'ensdf_processing/radlst/ref_radlst.inp'
    input_dict['output_radlst_file'] = tmp_path + '/tmp_radlst.rpt'
    input_dict['input_radlst_data_table'] = 'ensdf_processing/radlst/ref_mednew.dat'
    input_dict['output_ensdf_file'] = tmp_path + '/tmp_ensdf.rpt'
    output_dict = ensdf_processing.radlist(input_dict)
    ref_output_radlst_file = 'ensdf_processing/radlst/ref_radlst.rpt'
    ref_output_ensdf_file = 'ensdf_processing/radlst/ref_ensdf.rpt'
    # exceptions contain lines in the ouptut that can have a tolerable precision difference
    radlst_exceptions = [[1, '1PROGRAM RADLST 5.5 [ 5-OCT-88].  RUN ON'], [3, 66], [3, 135], [3, 713], [3, 714], [3, 760],[3,  944]]
    ensdf_exceptions = [[3, 341], [3, 351], [3, 357]]
    d_radlst = file_comp(input_dict['output_radlst_file'], ref_output_radlst_file, radlst_exceptions)
    d_ensdf = file_comp(input_dict['output_ensdf_file'], ref_output_ensdf_file, ensdf_exceptions)
    cleanup_tmp()
    os.remove('atomic.dat')

def test_ruler():
    create_tmp()
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/ruler/ref_ruler.inp'
    input_dict['output_report_file'] = tmp_path + '/tmp_ruler.rpt'
    input_dict['mode_of_operation'] = 'R'
    input_dict['assumed_dcc_theory'] = '1.4'
    output_dict = ensdf_processing.ruler(input_dict)
    ref_output = 'ensdf_processing/ruler/ref_ruler.rpt'
    exceptions = [[1, '         INPUT FILE:'], [1, 'RULER Version 3.2d [20-Jan-2009]']]
    d_report = file_comp(input_dict['output_report_file'], ref_output, exceptions)
    cleanup_tmp()

def create_tmp():
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

def cleanup_tmp():
    shutil.rmtree(tmp_path)

def file_comp(file_out, file_ref, exceptions):
    '''
    xceptions format: [type, options]
        type 1: prefix of length n.
            options: 'prefix'.
        type 2: general line ignore.
            options: 'prefix'
        type 3: minor precision issue ignore.
            options: line number of line with potential precision issue.
        type 4: carriage return vs. non standard return type.
            options: line number of return.
    '''
    f_out = open(file_out, 'r')
    f_ref = open(file_ref, 'r')
    diff_lines = numpy.array([])
    line_num = 0
    for line_out in f_out:
        line_ref = f_ref.readline()
        if(line_ref != line_out):
            ignore = False
            for i in range(0, len(exceptions)):
                if exceptions[i][0] == 1:
                    if line_out[0:len(exceptions[i][1])] == exceptions[i][1]:
                        ignore = True
                elif exceptions[i][0] == 2:
                    if exceptions[i][1] in line_out:
                          ignore = True
                elif exceptions[i][0] == 3:
                    # ignores select lines to allow for tolerable differences in output precision
                    if exceptions[i][1] == line_num:
                        ignore = True
                elif exceptions[i][0] == 4:
                    if len(line_ref[:-1]) == len(line_out):
                    # special exception for lines with possible carriage return instead of standard 
                    #line feed return
                        if line_ref[:-2] == line_out[:-1]:
                            if map(bin,bytearray(line_ref[len(line_ref)-1])) \
                            == map(bin,bytearray(line_out[len(line_out)-1])):
                                ignore = True
            if not ignore:
                raise Exception('ENSDF Processing: Incorrect output generated, file: ' + file_ref)
        line_num = line_num + 1
    f_out.close()
    f_ref.close()
    diff_dict = {}
    diff_dict['differences_lines'] = diff_lines
    return diff_dict

#  nose.runmodule()
if __name__ == "__main__":
    alphad = test_alphad()
    gabs = test_gabs()
    logft = test_logft()
    radd = test_radd()
    radlst = test_radlist()
    ruler = test_ruler()
