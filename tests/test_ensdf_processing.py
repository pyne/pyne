import os, filecmp, numpy
import nose
from nose.plugins.skip import Skip, SkipTest
import pyne
from pyne import ensdf_processing
#except:
#  raise SkipTest

def test_alphad():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/alphad.inp'
    input_dict['report_file'] = 'ensdf_processing/alphad.rpt'
    input_dict['rewrite_input_with_hinderance_factor'] = 1
    input_dict['output_file'] = 'ensdf_processing/alphad.out'
    output_dict = ensdf_processing.alphad(input_dict)
    print filecmp.cmp('ensdf_processing/alphad.rpt','ensdf_processing/alphad_correct.rpt')

    d_report = comp_file_with_date_difference('ensdf_processing/alphad.rpt','ensdf_processing/alphad_correct.rpt',0)


def test_delta():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/delta.dat'
    input_dict['output_file'] = 'ensdf_processing/delta.out'
    output_dict = ensdf_processing.delta(input_dict)
    print filecmp.cmp('ensdf_processing/delta.out','ensdf_processing/compare/delta_ref.rpt')
    d_report = comp_file_with_date_difference('ensdf_processing/alphad.rpt','ensdf_processing/alphad_correct.rpt',0)

def test_brick():
    print("not working..")
    print("brick only has executable no source..")

#def test_gabs_80Br():
#    input_dict = {}
#    input_dict['input_file'] = 'ensdf_processing/gabs_80Br.in'
#    input_dict['output_file'] = 'ensdf_processing/gabs_80Br.rpt'
#    input_dict['dataset_file'] = 'ensdf_processing/gabs_80Br.new'
#    output_dict = ensdf_processing.gabs(input_dict)
#    d_report1 = comp_file_with_date_difference(input_dict['output_file'],'ensdf_processing/compare/gabs_80Br_ref.rpt',0)
#    d_report2 = comp_file_with_date_difference(input_dict['dataset_file'],'ensdf_processing/compare/gabs_80Br_ref.new',0)

def test_gtol():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gtol/gtol.inp'
    input_dict['report_file'] = 'ensdf_processing/gtol/gtol.rpt'
    input_dict['new_ensdf_file_with_results'] = 1
    input_dict['output_file'] = 'ensdf_processing/gtol/gtol.out'
    input_dict['supress_gamma_comparison'] = 0
    input_dict['supress_intensity_comparison'] = 0
    input_dict['dcc_theory_percent'] = 3
    output_dict = ensdf_processing.gtol(input_dict)
    ref_output_report = 'ensdf_processing/gtol/gtol_ref.rpt'
    ref_output = 'ensdf_processing/gtol/gtol_ref.out'
    #d_report = comp_file_with_date_difference(input_dict['report_file'],ref_output_report,0)
    #d_output = comp_file_with_date_difference(input_dict['output_file'],ref_output,0)

def test_bldhst():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/bldhst/bldhst_iccseq.dat'
    input_dict['output_table_file'] = 'ensdf_processing/bldhst/bldhst_icctbl.dat'
    input_dict['output_index_file'] = 'ensdf_processing/bldhst/bldhst_iccndx.dat'
    output_dict = ensdf_processing.bldhst(input_dict)
    ref_table = 'ensdf_processing/bldhst/ref_icctbl.dat'
    ref_index = 'ensdf_processing/bldhst/ref_iccndx.dat'
    d_table = comp_file_with_date_difference(input_dict['output_table_file'],ref_table,0)
    d_index = comp_file_with_date_difference(input_dict['output_index_file'],ref_index,0)

def test_hsicc():
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsicc/hsicc_data.tst'
    input_dict['icc_index'] = 'ensdf_processing/hsicc/hsicc_iccndx.dat'
    input_dict['icc_table'] = 'ensdf_processing/hsicc/hsicc_icctbl.dat'
    input_dict['complete_report'] = 'ensdf_processing/hsicc/out_hsicc_hscalc.lst'
    input_dict['new_card_deck'] = 'ensdf_processing/hsicc/out_hsicc_cards.new'
    input_dict['comparison_report'] = 'ensdf_processing/hsicc/out_hsicc_compar.lst'
    input_dict['is_multipol_known'] = 'Y'
    output_dict = ensdf_processing.hsicc(input_dict)
    ref_report = 'ensdf_processing/hsicc/ref_hscalc.lst'
    ref_card_deck = 'ensdf_processing/hsicc/ref_cards.new'
    ref_comparison_report = 'ensdf_processing/hsicc/ref_compar.lst'
    d_report = comp_file_with_date_difference(input_dict['complete_report'],ref_report,0)
    d_card_deck = comp_file_with_date_difference(input_dict['new_card_deck'],ref_card_deck,0)
    d_comparison_report = comp_file_with_date_difference(input_dict['comparison_report'],ref_comparison_report,0)

def test_hsmrg():
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsicc_data.tst'
    input_dict['card_deck'] = 'ensdf_processing/cards.new'
    input_dict['merged_data_deck'] = 'ensdf_processing/cards.mrg'
    output_dict = ensdf_processing.hsmrg(input_dict)
        

def test_seqhst():
    input_dict = {}
    input_dict['binary_table_input_file'] = 'ensdf_processing/hsicc_icctbl.dat'
    input_dict['sequential_output_file'] = 'ensdf_processing/seqhst_iccseq.dat'
    output_dict = ensdf_processing.seqhst(input_dict)

def test_logft_functional():
    input_dict = {}
    input_dict['input_data_set'] = 'ensdf_processing/logft_data.tst'
    input_dict['output_report'] = 'ensdf_processing/logft.rpt'
    input_dict['data_table'] = 'ensdf_processing/logft.dat'
    input_dict['output_data_set'] = 'ensdf_processing/logft.new'
    output_dict = ensdf_processing.logft(input_dict)
    #check output files present?

def test_logft_outputs():
    input_dict = {}
    input_dict['input_data_set'] = 'ensdf_processing/logft_data.tst'
    input_dict['output_report'] = 'ensdf_processing/logft.rpt'
    input_dict['data_table'] = 'ensdf_processing/logft.dat'
    input_dict['output_data_set'] = 'ensdf_processing/logft.new'
    output_dict = ensdf_processing.logft(input_dict)
    ref_output_report = 'ensdf_processing/compare/logft_ref.rpt'
    ref_output_data_set = 'ensdf_processing/compare/logft_ref.new'
    d_report = comp_file_with_date_difference(input_dict['output_report'],ref_output_report,0)
    d_data = comp_file_with_date_difference(input_dict['output_data_set'], ref_output_data_set,0)

    # It is known that line 5 is the date.  Ignore difference at line 5
    if(d_report['discrepancy']):
        print d_report['differences_lines']
    if(d_data['discrepancy']):
        print d_data['differences_lines']

def test_pandora():
    print("not working")

def test_radlist():
    print("not working")

def test_ruler():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/ruler.inp'
    input_dict['output_report_file'] = 'ensdf_processing/ruler.rpt'
    input_dict['mode_of_operation'] = 'R'
    input_dict['assumed_dcc_theory'] = '3'
    output_dict = ensdf_processing.ruler(input_dict)

def comp_file_with_date_difference(file_out, file_ref, num_diff_lines):
    f_out = open(file_out, 'r')
    f_ref = open(file_ref, 'r')
    diff_lines = numpy.array([])
    line_num = 0
    for line_out in f_out:
        line_ref = f_ref.readline()
        if(line_ref != line_out):
            print  line_num
            print '     line_out is: ' + line_out
            print '     line_ref is: ' + line_ref
            diff_lines = numpy.append(diff_lines, line_num)
        line_num = line_num + 1

    diff_dict = {}
    diff_dict['differences_lines'] = diff_lines
    return diff_dict
            


#  nose.runmodule()
if __name__ == "__main__":
    #a1 = test_alphad()
    #a = test_delta()
    #b = test_gabs_80Br()
    #c = test_gtol()
    #d = test_bldhst()
    #nc = test_hsicc()
    n = test_hsmrg()
    #l = test_seqhst()
    #z = test_logft_functional()
    #z = test_logft_outputs()
    # pandora test needed
    # radlist test needed
    #r = test_ruler()
