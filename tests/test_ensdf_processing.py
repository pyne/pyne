import os
import nose
from nose.plugins.skip import Skip, SkipTest
import pyne
#try:
from pyne import ensdf_processing
#except:
#  raise SkipTest

def test_alphad():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/alphad.inp'
    input_dict['report_file'] = 'ensdf_processing/alphad.rpt'
    input_dict['rewrite_input_with_hinderance_factor'] = 'Y'
    input_dict['output_file'] = 'ensdf_processing/alphad.out'
    output_dict = ensdf_processing.alphad(input_dict)

def test_delta():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/delta.dat'
    input_dict['output_file'] = 'ensdf_processing/delta.out'
    output_dict = ensdf_processing.delta(input_dict)

def test_brick():
    print("not working..")
    print("brick only has executable no source..")

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

def test_bldhst():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/bldhst_iccseq.dat'
    input_dict['output_table_file'] = 'ensdf_processing/bldhst_icctbl.dat'
    input_dict['output_index_file'] = 'ensdf_processing/bldhst_iccndx.dat'
    output_dict = ensdf_processing.bldhst(input_dict)

def test_hsicc():
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsicc_data.tst'
    input_dict['icc_index'] = 'ensdf_processing/hsicc_iccndx.dat'
    input_dict['icc_table'] = 'ensdf_processing/hsicc_icctbl.dat'
    input_dict['complete_report'] = 'ensdf_processing/hsicc_hscalc.lst'
    input_dict['new_card_deck'] = 'ensdf_processing/hsicc_cards.new'
    input_dict['comparison_report'] = 'ensdf_processing/hsicc_compar.lst'
    input_dict['is_multipol_known'] = 'Y'
    output_dict = ensdf_processing.hsicc(input_dict)

def test_hsmrg():
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsicc_data.tst'
    input_dict['card_deck'] = 'ensdf_processing/cards.new'
    input_dict['merged_data_deck'] = 'ensdf_processing/cards.mrg'
    output_dict = ensdf_processing.hsmrg(input_dict)
        

def test_seqhst():
    input_dict = {}
    input_dict['binary_table_input_file'] = 'ensdf_processing/seqhst_icctbl.dat'
    input_dict['sequential_output_file'] = 'ensdf_processing/seqhst_iccseq.dat'
    output_dict = ensdf_processing.seqhst(input_dict)

def test_logft():
    print("not working")

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


#if __name__ == "__main__":
#  nose.runmodule()
if __name__ == "__main__":
    a1 = test_alphad()
    a = test_delta()
    b = test_gabs()
    c = test_gtol()
    d = test_bldhst()
    nc = test_hsicc()
    n = test_hsmrg()
    l = test_seqhst()
    # logft test needed
    # pandora test needed
    # radlist test needed
    r = test_ruler()
