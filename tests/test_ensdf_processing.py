import filecmp, numpy, os, shutil, io
from pyne import ensdf_processing

import nose
from nose.plugins.skip import SkipTest

# path to folder for temporary test files.
tmp_path = "ensdf_processing/tmp"


def test_alphad():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/alphad/ref_a228.ens"
    input_dict["report_file"] = tmp_path + "/tmp_alphad.rpt"
    input_dict["rewrite_input_with_hinderance_factor"] = 1
    input_dict["output_file"] = tmp_path + "/tmp_alphad.out"
    try:
        output_dict = ensdf_processing.alphad(input_dict)
    except OSError:
        raise SkipTest
    exceptions = [[2, "DATE RUN"]]
    file_comp(
        input_dict["report_file"],
        "ensdf_processing/alphad/ref_a228.ens.alphad.rpt",
        exceptions,
    )
    cleanup_tmp()


def optional_t_bricc_interactive():
    create_tmp()
    input_dict = {}
    input_dict["input_type"] = "interactive"
    input_dict["element"] = "44"
    try:
        output_dict = ensdf_processing.bricc(input_dict)
    except OSError:
        raise SkipTest
    bricc_out_tmp = tmp_path + "/tmp_bricc_out.out"
    bricc_out_ref = "ensdf_processing/bricc/ref_bricc_44.out"
    bricc_outfile = open(tmp_path + "/tmp_bricc_out.out", "w+")
    bricc_outfile.write(output_dict["bricc_output"].decode("utf-8"))
    file_comp(bricc_out_tmp, bricc_out_ref, [])
    cleanup_tmp()


def optional_t_bricc_evaluation():
    create_tmp()
    input_dict = {}
    input_dict["input_type"] = "evaluation"
    input_dict["input_file"] = "ensdf_processing/bricc/ref_a228.ens"
    input_dict["BrIccNH"] = 0
    input_dict["calculation_report"] = tmp_path + "/tmp_bricc_calculation_report"
    input_dict["comparison_report"] = tmp_path + "/tmp_bricc_comparison_report"
    try:
        output_dict = ensdf_processing.bricc(input_dict)
    except OSError:
        raise SkipTest
    bricc_comparison_ref = "ensdf_processing/bricc/ref_a228_comparison_report"
    file_comp(input_dict["comparison_report"], bricc_comparison_ref, [])
    cleanup_tmp()


def test_bldhst():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/bldhst/ref_bldhst_iccseq.dat"
    input_dict["output_table_file"] = tmp_path + "/tmp_bldhst_icctbl.dat"
    input_dict["output_index_file"] = tmp_path + "/tmp_bldhst_iccndx.dat"
    try:
        output_dict = ensdf_processing.bldhst(input_dict)
    except OSError:
        raise SkipTest
    ref_table = "ensdf_processing/bldhst/ref_icctbl.dat"
    ref_index = "ensdf_processing/bldhst/ref_iccndx.dat"
    d_table = file_comp(input_dict["output_table_file"], ref_table, [])
    d_index = file_comp(input_dict["output_index_file"], ref_index, [])
    cleanup_tmp()


def test_delta():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/delta/ref_inp.dat"
    input_dict["output_file"] = tmp_path + "/tmp_delta.dat"
    try:
        output_dict = ensdf_processing.delta(input_dict)
    except OSError:
        raise SkipTest
    # exceptions contain lines in the ouptut that can have a tolerable
    # precision difference
    exceptions = [
        [3, 82],
        [3, 89],
        [3, 119],
        [3, 202],
        [3, 209],
        [3, 213],
        [3, 217],
        [3, 229],
        [3, 232],
        [3, 236],
        [3, 243],
        [3, 318],
        [3, 355],
        [3, 458],
    ]
    d_ouptut = file_comp(
        input_dict["output_file"], "ensdf_processing/delta/ref_delta.rpt", exceptions
    )
    cleanup_tmp()


# Tests gabs output for 80Br sample input.  Date and file names are expected to
# be different in the test output and reference file sets.
def optional_t_gabs():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/gabs/ref_gabs_80Br.in"
    input_dict["output_file"] = tmp_path + "/tmp_gabs_80Br.rpt"
    input_dict["dataset_file"] = tmp_path + "/tmp_gabs_80Br.new"
    output_dict = ensdf_processing.gabs(input_dict)
    exceptions_output = [
        [4, 0],
        [1, "  * * * GABS Version 11 "],
        [1, "        Current date: "],
        [1, "        ENSDF input file: "],
        [1, "        new ENSDF file:"],
    ]
    exceptions_dataset = [[4, 0]]
    d_report1 = file_comp(
        input_dict["output_file"],
        "ensdf_processing/gabs/ref_gabs_80Br.rpt",
        exceptions_output,
    )
    d_report2 = file_comp(
        input_dict["dataset_file"],
        "ensdf_processing/gabs/ref_gabs_80Br.new",
        exceptions_dataset,
    )
    cleanup_tmp()


def test_gtol():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/gtol/ref_gtol.inp"
    input_dict["report_file"] = tmp_path + "/tmp_gtol.rpt"
    input_dict["new_ensdf_file_with_results"] = 0
    input_dict["output_file"] = tmp_path + "/tmp_gtol.out"
    input_dict["supress_gamma_comparison"] = 1
    input_dict["supress_intensity_comparison"] = 1
    input_dict["dcc_theory_percent"] = 1.4
    try:
        output_dict = ensdf_processing.gtol(input_dict)
    except OSError:
        raise SkipTest
    ref_output_report = "ensdf_processing/gtol/ref_gtol.rpt"
    exceptions = [[1, "DATE:"], [1, "INPUT-FILE name:"], [1, "TIME:"]]
    d_report = file_comp(input_dict["report_file"], ref_output_report, exceptions)
    cleanup_tmp()


def test_hsicc():
    create_tmp()
    input_dict = {}
    input_dict["data_deck"] = "ensdf_processing/hsicc/ref_hsicc_data.tst"
    input_dict["icc_index"] = "ensdf_processing/hsicc/ref_hsicc_iccndx.dat"
    input_dict["icc_table"] = "ensdf_processing/hsicc/ref_hsicc_icctbl.dat"
    input_dict["complete_report"] = tmp_path + "/tmp_out_hsicc_hscalc.lst"
    input_dict["new_card_deck"] = tmp_path + "/tmp_out_hsicc_cards.new"
    input_dict["comparison_report"] = tmp_path + "/tmp_out_hsicc_compar.lst"
    input_dict["is_multipol_known"] = "Y"
    try:
        output_dict = ensdf_processing.hsicc(input_dict)
    except OSError:
        raise SkipTest
    ref_report = "ensdf_processing/hsicc/ref_hscalc.lst"
    ref_card_deck = "ensdf_processing/hsicc/ref_cards.new"
    ref_comparison_report = "ensdf_processing/hsicc/ref_compar.lst"

    exceptions = [[3, 55], [3, 70], [3, 83], [3, 107], [3, 131], [3, 151]]
    d_report = file_comp(input_dict["complete_report"], ref_report, exceptions)
    d_card_deck = file_comp(input_dict["new_card_deck"], ref_card_deck, [])
    d_comparison_report = file_comp(
        input_dict["comparison_report"], ref_comparison_report, []
    )
    cleanup_tmp()


def test_hsmrg():
    create_tmp()
    input_dict = {}
    input_dict["data_deck"] = "ensdf_processing/hsmrg/ref_hsmrg_data.tst"
    input_dict["card_deck"] = "ensdf_processing/hsmrg/ref_hsmrg_cards.new"
    input_dict["merged_data_deck"] = tmp_path + "/tmp_out_cards.mrg"
    try:
        output_dict = ensdf_processing.hsmrg(input_dict)
    except OSError:
        raise SkipTest
    ref_deck = "ensdf_processing/hsmrg/ref_cards.mrg"
    d_report = file_comp(input_dict["merged_data_deck"], ref_deck, [])
    cleanup_tmp()


def test_seqhst():
    create_tmp()
    input_dict = {}
    input_dict[
        "binary_table_input_file"
    ] = "ensdf_processing/seqhst/ref_seqhst_icctbl.dat"
    input_dict["sequential_output_file"] = tmp_path + "/tmp_out_iccseq.dat"
    try:
        output_dict = ensdf_processing.seqhst(input_dict)
    except OSError:
        raise SkipTest
    ref_sequence = "ensdf_processing/seqhst/ref_iccseq.dat"
    d_report = file_comp(input_dict["sequential_output_file"], ref_sequence, [])
    cleanup_tmp()


def test_logft():
    create_tmp()
    input_dict = {}
    input_dict["input_data_set"] = "ensdf_processing/logft/ref_logft.inp"
    input_dict["output_report"] = tmp_path + "/tmp_logft.rpt"
    input_dict["data_table"] = "ensdf_processing/logft/ref_logft.dat"
    input_dict["output_data_set"] = tmp_path + "/tmp_logft.new"
    try:
        output_dict = ensdf_processing.logft(input_dict)
    except OSError:
        raise SkipTest
    ref_output_data_set = "ensdf_processing/logft/ref_logft.new"
    d_data = file_comp(input_dict["output_data_set"], ref_output_data_set, [])
    cleanup_tmp()


def test_radd():
    create_tmp()
    input_dict = {}
    input_dict["atomic_number"] = "86"
    input_dict["neutron_number"] = "113"
    input_dict["output_file"] = tmp_path + "/tmp_output.out"
    try:
        ensdf_processing.radd(input_dict)
    except OSError:
        raise SkipTest
    ref_output = "ensdf_processing/radd/ref_output.out"
    d_report = file_comp(input_dict["output_file"], ref_output, [])
    cleanup_tmp()


def optional_t_radlist():
    create_tmp()
    input_dict = {}
    input_dict["output_radiation_listing"] = "Y"
    input_dict["output_ensdf_like_file"] = "N"
    input_dict["output_file_for_nudat"] = "N"
    input_dict["output_mird_listing"] = "N"
    input_dict["calculate_continua"] = "N"
    input_dict["input_file"] = "ensdf_processing/radlst/ref_radlst.inp"
    input_dict["output_radlst_file"] = tmp_path + "/tmp_radlst.rpt"
    input_dict["input_radlst_data_table"] = "ensdf_processing/radlst/ref_mednew.dat"
    input_dict["output_ensdf_file"] = tmp_path + "/tmp_ensdf.rpt"
    output_dict = ensdf_processing.radlist(input_dict)
    ref_output_radlst_file = "ensdf_processing/radlst/ref_radlst.rpt"
    ref_output_ensdf_file = "ensdf_processing/radlst/ref_ensdf.rpt"
    # exceptions contain lines in the ouptut that can have a tolerable precision difference
    radlst_exceptions = [
        [1, "1PROGRAM RADLST 5.5 [ 5-OCT-88].  RUN ON"],
        [3, 66],
        [3, 135],
        [3, 713],
        [3, 714],
        [3, 760],
        [3, 944],
    ]
    ensdf_exceptions = [[3, 341], [3, 351], [3, 357]]
    d_radlst = file_comp(
        input_dict["output_radlst_file"], ref_output_radlst_file, radlst_exceptions
    )
    d_ensdf = file_comp(
        input_dict["output_ensdf_file"], ref_output_ensdf_file, ensdf_exceptions
    )
    cleanup_tmp()
    os.remove("atomic.dat")


def test_ruler():
    create_tmp()
    input_dict = {}
    input_dict["input_file"] = "ensdf_processing/ruler/ref_ruler.inp"
    input_dict["output_report_file"] = tmp_path + "/tmp_ruler.rpt"
    input_dict["mode_of_operation"] = "R"
    input_dict["assumed_dcc_theory"] = "1.4"
    try:
        output_dict = ensdf_processing.ruler(input_dict)
    except OSError:
        raise SkipTest
    ref_output = "ensdf_processing/ruler/ref_ruler.rpt"
    exceptions = [[1, "         INPUT FILE:"], [1, "RULER Version 3.2d [20-Jan-2009]"]]
    d_report = file_comp(input_dict["output_report_file"], ref_output, exceptions)
    cleanup_tmp()


def create_tmp():
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)


def cleanup_tmp():
    shutil.rmtree(tmp_path)


def file_comp(file_out, file_ref, exceptions):
    """
    Exceptions format: [type, options]
        type 1: prefix of length n.
            options: 'prefix'.
        type 2: general line ignore.
            options: 'prefix'
        type 3: minor precision issue ignore.
            options: line number of line with potential precision issue.
        type 4: carriage return vs. non standard return type.
            options: line number of return.
    """
    f_out = io.open(file_out, "r", encoding="latin-1")
    f_ref = io.open(file_ref, "r", encoding="latin-1")
    diff_lines = numpy.array([])
    line_num = 0
    for line_out in f_out:
        line_ref = f_ref.readline()
        if line_ref != line_out:
            ignore = False
            for i in range(0, len(exceptions)):
                if exceptions[i][0] == 1:
                    if line_out[0 : len(exceptions[i][1])] == exceptions[i][1]:
                        ignore = True
                elif exceptions[i][0] == 2:
                    if exceptions[i][1] in line_out:
                        ignore = True
                elif exceptions[i][0] == 3:
                    # ignores select lines to allow for tolerable differences
                    # in output precision
                    if exceptions[i][1] == line_num:
                        ignore = True
                elif exceptions[i][0] == 4:
                    if len(line_ref[:-1]) == len(line_out):
                        # special exception for lines with possible carriage return
                        # instead of standard line feed return
                        if line_ref[:-2] == line_out[:-1]:
                            ignore = True
            if not ignore:
                raise Exception(
                    "ENSDF Processing: Incorrect output generated, file: " + file_ref
                )
        line_num = line_num + 1
    f_out.close()
    f_ref.close()
    diff_dict = {}
    diff_dict["differences_lines"] = diff_lines
    return diff_dict


#  nose.runmodule()
if __name__ == "__main__":
    alphad = test_alphad()
    bldhst = test_bldhst()
    delta = test_delta()
    gtol = test_gtol()
    hsicc = test_hsicc()
    hsmrg = test_hsmrg()
    seqhst = test_seqhst()
    logft = test_logft()
    radd = test_radd()
    ruler = test_ruler()

    ### The following code runs the additional setup code and runs tests
    ### for each.  They are not included in the main test suite because
    ### downloading the additional execuables can take upwards of 120 seconds.
    # ensdf_processing.setup_additional_downloads()
    # gabs = optional_t_gabs()
    # bricc_1 = optional_t_bricc_interactive()
    # bricc_2 = optional_t_bricc_evaluation()
    # radlist = optional_t_radlist()
