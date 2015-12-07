import os, filecmp, numpy
import nose
from nose.plugins.skip import Skip, SkipTest
import pyne
from pyne import ensdf_processing
#except:
#  raise SkipTest

def test_alphad():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/alphad/ref_a228.ens'
    input_dict['report_file'] = 'ensdf_processing/alphad/tmp_alphad.rpt'
    input_dict['rewrite_input_with_hinderance_factor'] = 1
    input_dict['output_file'] = 'ensdf_processing/alphad/tmp_alphad.out'
    output_dict = ensdf_processing.alphad(input_dict)
    exceptions = [[2, 'DATE RUN']]
    file_comp('ensdf_processing/alphad/tmp_alphad.rpt','ensdf_processing/alphad/ref_a228.ens.alphad.rpt', exceptions)
    #print filecmp.cmp('ensdf_processing/alphad/tmp_alphad.rpt','ensdf_processing/alphad/ref_a228.ens.alphad.rpt')
    #print comp_file_with_date_difference('ensdf_processing/alphad/tmp_alphad.rpt','ensdf_processing/alphad/ref_a228.ens.alphad.rpt',0)

def test_delta():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/delta/ref_inp.dat'
    input_dict['output_file'] = 'ensdf_processing/delta/tmp_delta.dat'
    output_dict = ensdf_processing.delta(input_dict)
    print filecmp.cmp(input_dict['output_file'],'ensdf_processing/delta/ref_delta.rpt')
    print file_comp(input_dict['output_file'],'ensdf_processing/delta/ref_delta.rpt', [])

def test_bricc():
    print("not working..")
    print("bricc only has executable no source..")
    input_dict = {}
    output_dict = ensdf_processing.bricc(input_dict)

def test_gabs_80Br():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gabs/ref_gabs_80Br.in'
    input_dict['output_file'] = 'ensdf_processing/gabs/tmp_gabs_80Br.rpt'
    input_dict['dataset_file'] = 'ensdf_processing/gabs/tmp_gabs_80Br.new'
    output_dict = ensdf_processing.gabs(input_dict)
    #d_report1 = file_comp(input_dict['output_file'],'ensdf_processing/gabs/ref_gabs_80Br.rpt',[])
    #d_report2 = file_comp(input_dict['dataset_file'],'ensdf_processing/gabs/ref_gabs_80Br.new',[])

def test_gtol():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/gtol/ref_gtol.inp'
    input_dict['report_file'] = 'ensdf_processing/gtol/tmp_gtol.rpt'
    input_dict['new_ensdf_file_with_results'] = 0
    input_dict['output_file'] = 'ensdf_processing/gtol/tmp_gtol.out'
    input_dict['supress_gamma_comparison'] = 1
    input_dict['supress_intensity_comparison'] = 1
    input_dict['dcc_theory_percent'] = 1.4
    output_dict = ensdf_processing.gtol(input_dict)
    ref_output_report = 'ensdf_processing/gtol/ref_gtol.rpt'
    exceptions = [[1, 'DATE:'], [1, 'INPUT-FILE name:'], [1, 'TIME:']]
    d_report = file_comp(input_dict['report_file'], ref_output_report, exceptions)

def test_bldhst():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/bldhst/ref_bldhst_iccseq.dat'
    input_dict['output_table_file'] = 'ensdf_processing/bldhst/tmp_bldhst_icctbl.dat'
    input_dict['output_index_file'] = 'ensdf_processing/bldhst/tmp_bldhst_iccndx.dat'
    output_dict = ensdf_processing.bldhst(input_dict)
    ref_table = 'ensdf_processing/bldhst/ref_icctbl.dat'
    ref_index = 'ensdf_processing/bldhst/ref_iccndx.dat'
    d_table = file_comp(input_dict['output_table_file'], ref_table, [])
    d_index = file_comp(input_dict['output_index_file'], ref_index, [])

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
    d_report = file_comp(input_dict['complete_report'], ref_report, [])
    print d_report
    d_card_deck = file_comp(input_dict['new_card_deck'], ref_card_deck, [])
    print d_card_deck
    d_comparison_report = file_comp(input_dict['comparison_report'], ref_comparison_report, [])
    print d_comparison_report

def test_hsmrg():
    input_dict = {}
    input_dict['data_deck'] = 'ensdf_processing/hsmrg/hsmrg_data.tst'
    input_dict['card_deck'] = 'ensdf_processing/hsmrg/hsmrg_cards.new'
    input_dict['merged_data_deck'] = 'ensdf_processing/hsmrg/out_cards.mrg'
    output_dict = ensdf_processing.hsmrg(input_dict)
    ref_deck = 'ensdf_processing/hsmrg/ref_cards.mrg'
    d_report = file_comp(input_dict['merged_data_deck'], ref_deck, [])

def test_seqhst():
    input_dict = {}
    input_dict['binary_table_input_file'] = 'ensdf_processing/seqhst/seqhst_icctbl.dat'
    input_dict['sequential_output_file'] = 'ensdf_processing/seqhst/out_iccseq.dat'
    output_dict = ensdf_processing.seqhst(input_dict)
    ref_sequence = 'ensdf_processing/seqhst/ref_iccseq.dat'
    d_report = file_comp(input_dict['sequential_output_file'], ref_sequence, [])

def test_logft():
    input_dict = {}
    input_dict['input_data_set'] = 'ensdf_processing/logft_data.tst'
    input_dict['output_report'] = 'ensdf_processing/logft.rpt'
    input_dict['data_table'] = 'ensdf_processing/logft.dat'
    input_dict['output_data_set'] = 'ensdf_processing/logft.new'
    output_dict = ensdf_processing.logft(input_dict)
    ref_output_data_set = 'ensdf_processing/compare/logft_ref.new'
    d_data = file_comp(input_dict['output_data_set'], ref_output_data_set, [])

def test_pandora():
    input_dict = {}
    input_dict['input_data_set'] = 'ensdf_processing/pandora/pandora.inp'
    input_dict['output_err'] = 'ensdf_processing/pandora/tmp_pandora.err'
    input_dict['output_gam'] = 'ensdf_processing/pandora/tmp_pandora.gam'
    input_dict['output_gle'] = 'ensdf_processing/pandora/tmp_pandora.gle'
    input_dict['output_lev'] = 'ensdf_processing/pandora/tmp_pandora.lev'
    input_dict['output_rad'] = 'ensdf_processing/pandora/tmp_pandora.rad'
    input_dict['output_xrf'] = 'ensdf_processing/pandora/tmp_pandora.xrf'
    input_dict['output_out'] = 'ensdf_processing/pandora/tmp_pandora.out'
    output_dict = ensdf_processing.pandora(input_dict)
    print output_dict

def test_radd():
    input_dict = {}
    input_dict['atomic_number'] = '86'
    input_dict['neutron_number'] = '113'
    input_dict['output_file'] = 'ensdf_processing/radd/tmp_output.out'
    ensdf_processing.radd(input_dict)
    ref_output = 'ensdf_processing/radd/ref_output.out'
    d_report = file_comp(input_dict['output_file'], ref_output, [])

def test_radlist():
    print("implement once download finished")
    #input_dict = {}
    #output_dict = ensdf_processing.radlist(input_dict)

def test_ruler():
    input_dict = {}
    input_dict['input_file'] = 'ensdf_processing/ruler/ref_ruler.inp'
    input_dict['output_report_file'] = 'ensdf_processing/ruler/tmp_ruler.rpt'
    input_dict['mode_of_operation'] = 'R'
    input_dict['assumed_dcc_theory'] = '1.4'
    output_dict = ensdf_processing.ruler(input_dict)
    ref_output = 'ensdf_processing/ruler/ref_ruler.rpt'
    exceptions = [[1, '         INPUT FILE:'], [1, 'RULER Version 3.2d [20-Jan-2009]']]
    d_report = file_comp(input_dict['output_report_file'], ref_output, exceptions)

def file_comp(file_out, file_ref, exceptions):
    # exceptions format: [type, options]
    #   type 1: prefix of length n
    #       options: 'prefix'
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
            if not ignore:
            	#raise Exception('ENSDF Processing: Incorrect output generated, file: ' + file_ref)
            	print 'difference found %i', line_num
            	print '     line_out is: ' + line_out
            	print '     line_ref is: ' + line_ref
            	print len(line_out)
            	print len(line_ref)
            	diff_lines = numpy.append(diff_lines, line_num)
        line_num = line_num + 1
    f_out.close()
    f_ref.close()
    diff_dict = {}
    diff_dict['differences_lines'] = diff_lines
    #if len(diff_lines) == 0:
    	#print file_ref + ' is the same'
    return diff_dict

#  nose.runmodule()
if __name__ == "__main__":
    alphad = test_alphad()
    #b = test_bricc() # FINISH DOWNLOAD
    #g = test_gabs_80Br() # FINISH DOWNLOAD
    c = test_gtol()
    #d1 = test_delta() # DIFFERENCES
    d = test_bldhst()
    #nc = test_hsicc() # DIFFERENCES
    n = test_hsmrg()
    l = test_seqhst()
    z = test_logft()
    c = test_radd()
    #p = test_pandora() # FIX BUG
    #ra = test_radlist() # FINISH DOWNLOAD
    r = test_ruler()
