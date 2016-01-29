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

def test_delta():
    create_tmp()
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/delta/ref_inp.dat'
    input_dict['output_file'] = tmp_path + '/tmp_delta.dat'
    output_dict = ensdf_processing.delta(input_dict)
    # exceptions contain lines in the ouptut that can have a tolerable precision difference
    exceptions = [[3, 82], [3, 89], [3, 119], [3, 202], [3, 209], [3, 213],[3,  217],[3,  229], \
                 [3,  232], [3, 236], [3, 243], [3, 318], [3, 355], [3, 458]]
    d_ouptut = file_comp(input_dict['output_file'],'ensdf_processing/delta/ref_delta.rpt', exceptions)
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

def test_gtol():
    create_tmp()
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gtol/ref_gtol.inp'
    input_dict['report_file'] = tmp_path + '/tmp_gtol.rpt'
    input_dict['new_ensdf_file_with_results'] = 0
    input_dict['output_file'] = tmp_path + '/tmp_gtol.out'
    input_dict['supress_gamma_comparison'] = 1
    input_dict['supress_intensity_comparison'] = 1
    input_dict['dcc_theory_percent'] = 1.4
    output_dict = ensdf_processing.gtol(input_dict)
    ref_output_report = 'ensdf_processing/gtol/ref_gtol.rpt'
    exceptions = [[1, 'DATE:'], [1, 'INPUT-FILE name:'], [1, 'TIME:']]
    d_report = file_comp(input_dict['report_file'], ref_output_report, exceptions)
    cleanup_tmp()

def test_hsicc():
    create_tmp()
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsicc/ref_hsicc_data.tst'
    input_dict['icc_index'] = 'ensdf_processing/hsicc/ref_hsicc_iccndx.dat'
    input_dict['icc_table'] = 'ensdf_processing/hsicc/ref_hsicc_icctbl.dat'
    input_dict['complete_report'] = tmp_path + '/tmp_out_hsicc_hscalc.lst'
    input_dict['new_card_deck'] = tmp_path + '/tmp_out_hsicc_cards.new'
    input_dict['comparison_report'] = tmp_path + '/tmp_out_hsicc_compar.lst'
    input_dict['is_multipol_known'] = 'Y'
    output_dict = ensdf_processing.hsicc(input_dict)
    ref_report = 'ensdf_processing/hsicc/ref_hscalc.lst'
    ref_card_deck = 'ensdf_processing/hsicc/ref_cards.new'
    ref_comparison_report = 'ensdf_processing/hsicc/ref_compar.lst'

    exceptions = [[3, 55], [3, 70], [3, 83], [3, 107], [3, 131], [3, 151]]
    d_report = file_comp(input_dict['complete_report'], ref_report, exceptions)
    d_card_deck = file_comp(input_dict['new_card_deck'], ref_card_deck, [])
    d_comparison_report = file_comp(input_dict['comparison_report'], ref_comparison_report, [])
    cleanup_tmp()

def test_hsmrg():
    create_tmp()
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsmrg/ref_hsmrg_data.tst'
    input_dict['card_deck'] = 'ensdf_processing/hsmrg/ref_hsmrg_cards.new'
    input_dict['merged_data_deck'] = tmp_path + '/tmp_out_cards.mrg'
    output_dict = ensdf_processing.hsmrg(input_dict)
    ref_deck = 'ensdf_processing/hsmrg/ref_cards.mrg'
    d_report = file_comp(input_dict['merged_data_deck'], ref_deck, [])
    cleanup_tmp()

def test_seqhst():
    create_tmp()
    input_dict = {}
    input_dict['binary_table_input_file'] = 'ensdf_processing/seqhst/ref_seqhst_icctbl.dat'
    input_dict['sequential_output_file'] = tmp_path + '/tmp_out_iccseq.dat'
    output_dict = ensdf_processing.seqhst(input_dict)
    ref_sequence = 'ensdf_processing/seqhst/ref_iccseq.dat'
    d_report = file_comp(input_dict['sequential_output_file'], ref_sequence, [])
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
    Exceptions format: [type, options]
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
    delta = test_delta()
    gtol = test_gtol()
    hsicc = test_hsicc()
    hsmrg = test_hsmrg()
    seqhst = test_seqhst()
    logft = test_logft()
    radd = test_radd()
    ruler = test_ruler()
