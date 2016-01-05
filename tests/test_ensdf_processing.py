import os, filecmp, numpy, shutil
import nose
from nose.plugins.skip import Skip, SkipTest
import pyne
from pyne import ensdf_processing

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
    d_report1 = file_comp(input_dict['output_file'],'ensdf_processing/gabs/ref_gabs_80Br.rpt',exceptions_output)
    d_report2 = file_comp(input_dict['dataset_file'],'ensdf_processing/gabs/ref_gabs_80Br.new',exceptions_dataset)
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
               		# special exception for lines with possible carriage return instead of standard line feed return
              			if line_ref[:-2] == line_out[:-1]:
	              			if map(bin,bytearray(line_ref[len(line_ref)-1])) == map(bin,bytearray(line_out[len(line_out)-1])):
        	      				ignore = True
            if not ignore:
            	raise Exception('ENSDF Processing: Incorrect output generated, file: ' + file_ref)
            	#print 'difference found %i', line_num
            	#print '     line_out is: ' + line_out
            	#print '     line_ref is: ' + line_ref
            	#print len(line_out)
            	#print len(line_ref)
            	#diff_lines = numpy.append(diff_lines, line_num)
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

