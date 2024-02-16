"""This module accesses various ensdf processing tools"""

import sys, os, shutil, subprocess, tarfile
from warnings import warn
from pyne.utils import QA_warn

try:
    import urllib.request as urllib
except ImportError:
    import urllib2 as urllib

if sys.version_info[0] > 2:
    basestring = str

QA_warn(__name__)


class SetupIncompleteError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return (
            "This ensdf analysis module requires additional setup.  Run \
                    pyne.ensdf_processing.setup_additional_downloads to download the \
                    additional executables and complete this modules instalation. \
                    (following executable not yet installed/downloaded: "
            + repr(self.value)
        )


def path_to_exe(exe_name):
    exe_path_abs, dp = os.path.split(os.path.abspath(__file__))
    exe_path_abs = os.path.join(exe_path_abs, exe_name)
    exe_path_abs = os.path.join("./", exe_path_abs)
    return exe_path_abs


def download_exe(exe_path, exe_url, compressed=0, decomp_path="", dl_size=0):
    if not os.path.exists(exe_path):
        msg = "Downloading {0!r} to {1!r}".format(exe_url, exe_path)
        print(msg)
        req = urllib.Request(exe_url, headers={"User-Agent": "Mozilla/5.0"})
        f = urllib.urlopen(req, timeout=30.0)
        try:
            html = f.read()
        finally:
            f.close()
        with open(exe_path, "wb") as f:
            f.write(html)

        # set proper permissions on newly downloaded file
        os.chmod(exe_path, 744)
        if compressed:
            os.chmod(exe_path, 777)
            tfile = tarfile.open(exe_path, "r:gz")
            tfile.extractall(decomp_path)


def setup_additional_downloads():
    print("Downloading BRICC executable")
    setup_bricc()
    print("Downloading GABS executable")
    setup_gabs()
    print("Downloading RADLIST executable")
    setup_radlist()


def setup_bricc():
    exe_path = path_to_exe("bricc")
    exe_dir = path_to_exe("")[:-1]
    compressed_exe_path = exe_path + ".tar.gz"

    bricc_url = (
        "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/"
        + "BrIcc/Linux/BriccV23-Linux.tgz"
    )
    decomp_exe_path = path_to_exe("")
    decomp_options = ["bricc", ".tgz", True]
    download_exe(
        compressed_exe_path,
        bricc_url,
        compressed=True,
        decomp_path=decomp_exe_path,
        dl_size=127232,
    )


def setup_gabs():
    exe_path = path_to_exe("gabs")
    gabs_url = "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/" + "gabs/unx/gabs"
    download_exe(exe_path, gabs_url, dl_size=8704)


def setup_radlist():
    exe_path = path_to_exe("radlist")
    radlist_url = (
        "http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/" + "radlst/unx/radlist"
    )
    download_exe(exe_path, radlist_url, dl_size=8704)


def verify_download_exe(exe_path):
    if not os.path.exists(exe_path):
        raise SetupIncompleteError(exe_path)


def alphad(inputdict_unchecked):
    """
    This function calculates the alpha hinderance factors and theoretical half
    lives for even even ground state transitions. (alphad readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            ensdf_input_file : string, input file
            output_file : string, file for output to be written to (doesn't
                          have to exist)

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if ALPHAD completes
        successfully.

    Full documentation explaining the details of the functionality and physics
    behind ALPHAD can be found at:
    http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/alphad/readme-alphad.pdf
    """
    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    report_file = inputdict_unchecked["report_file"]
    rewrite_hinderance = inputdict_unchecked["rewrite_input_with_hinderance_factor"]
    output_file = "alphad.out"
    if rewrite_hinderance == 1:
        output_file = inputdict_unchecked["output_file"]  # output file if report = yes
    exe_path = path_to_exe("alphad")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + report_file + "\n" + "Y" + "\n"
    if rewrite_hinderance == 1:
        inp = inp + "Y" + "\n" + output_file
    else:
        inp = inp + "N" + "\n"
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def bricc(inputdict_unchecked):
    """
    This function calculates the conversion electron, electron-positron pair conversion
    coefficients, and the E0 electron factors.

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            xx
            xx

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if BRICC completes
        successfully.
        Additional dictionary entries including:
            output_file_directory : string, the directory all produced bricc
                                    output files will be located.
            bricc_output : string, Only for interactive use: data printed to
                           command line.

    Notes
    -----
        All the various ouptput files bricc can generate are found in the
        output_file_directory

    """
    exe_path = path_to_exe("bricc")
    exe_dir = path_to_exe("")[:-1]
    compressed_exe_path = exe_path + ".tar.gz"

    verify_download_exe(compressed_exe_path)
    # check if BriIccHome environment variable has been set (needed by BRICC executable)
    if not os.environ.get("BrIccHome"):
        os.environ["BrIccHome"] = str(exe_dir)

    input_type = inputdict_unchecked["input_type"]
    output_dict = inputdict_unchecked
    inp = (
        (
            inputdict_unchecked["calculation_report"]
            if inputdict_unchecked.has_key("calculation_report")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["G_SG_records"]
            if inputdict_unchecked.has_key("G_SG_records")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["comparison_report"]
            if inputdict_unchecked.has_key("comparison_report")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["conversion_coefficients"]
            if inputdict_unchecked.has_key("conversion_coefficients")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["calculate_conversion_coefficients"]
            if inputdict_unchecked.has_key("calculate_conversion_coefficients")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["lowest_cc_value"]
            if inputdict_unchecked.has_key("lowest_cc_value")
            else ""
        )
        + "\n"
        + (
            inputdict_unchecked["assumed_mr_value"]
            if inputdict_unchecked.has_key("assumed_mr_value")
            else ""
        )
        + "\n\r\n"
    )

    if input_type == "interactive":
        input_element = inputdict_unchecked["element"]
        inp = input_element + "\n" + "exit" + "\n"
        proc = subprocess.Popen(
            [exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE
        )
        proc.stdin.write(inp.encode("utf-8"))
        output_dict["bricc_output"] = proc.communicate()[0]
        proc.stdin.close()
    elif input_type == "evaluation":
        input_file = inputdict_unchecked["input_file"]
        briccnh = inputdict_unchecked["BrIccNH"]
        if briccnh:
            proc = subprocess.Popen(
                exe_path,
                [input_file, "BrIccNH"],
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
            )
            proc.stdin.write(inp.encode("utf-8"))
            proc.communicate()[0]
            proc.stdin.close()
        else:
            proc = subprocess.Popen(
                [exe_path, input_file], stdout=subprocess.PIPE, stdin=subprocess.PIPE
            )
            proc.stdin.write(inp.encode("utf-8"))
            proc.communicate()[0]
            proc.stdin.close()

    output_dict["output_file_directory"] = exe_dir
    return output_dict


def bldhst(inputdict_unchecked):
    """
    This program builds a direct access file of the internal conversion
    coefficient table. (BLDHST readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, input ensdf file.
            output_table_file : string, desired output table file path.
            output_index_file : string, desired output index file path.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if DELTA completes successfully.
    """
    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    output_table_file = inputdict_unchecked["output_table_file"]
    output_index_file = inputdict_unchecked["output_index_file"]

    exe_path = path_to_exe("bldhst")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + output_table_file + "\n" + output_index_file
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def delta(inputdict_unchecked):
    """
    This function calculates the best values of mixing ratios based of its
    analysis of the angular correlation and conversion coefficient data.

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, path to input ensdf file.
            output_file : string, path to file for output write (doesn't have to
                          exist).

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if DELTA completes
        successfully.
    """
    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    output_file = inputdict_unchecked["output_file"]

    exe_path = path_to_exe("delta")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + output_file + "\n"
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def gabs(inputdict_unchecked):
    """
    This program calculates Gamma-ray absolute intensity and normalization
    (GABS readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, input ensdf file
            dataset_file : string, dataset file to be used
            output file : string, file for output to be written to
                          (doesn't have to exist)

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if GABS completes
        successfully.
    """
    exe_path = path_to_exe("gabs")
    verify_download_exe(exe_path)

    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    dataset_file = inputdict_unchecked["dataset_file"]
    output_file = inputdict_unchecked["output_file"]

    exe_path = path_to_exe("gabs")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + output_file + "\n" + "Y" + "\n" + dataset_file
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()


def gtol(inputdict_unchecked):
    """
    GTOL uses gamma-ray energies to derive a set of least-squares adjusted level
    energies.

    The net feeding at each level is calculated from the input gamma intensities
    and conversion coefficients. (GTOL readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, input ensdf file.
            report_file : string, desired gtol report file path.
            new_ensdf_file_with_results : boolean, if true then a new ensdf file
                                          with results will be created.
            output_file : string, desired gtol output file path.
            supress_gamma_comparison : boolean, if true the gamma comparison will
                                       be suppressed.
            dcc_theory_percent : double, specifies the dcc theory percentage to
                                 be used.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if GTOL completes successfully.
    """
    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    report_file = inputdict_unchecked["report_file"]
    new_out = inputdict_unchecked["new_ensdf_file_with_results"]
    output_file = inputdict_unchecked["output_file"]  # output file if report = yes
    supress_g = inputdict_unchecked["supress_gamma_comparison"]
    supress_ic = inputdict_unchecked["supress_intensity_comparison"]
    dcc_theory = inputdict_unchecked["dcc_theory_percent"]

    exe_path = path_to_exe("gtol")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + report_file + "\n"
    if new_out:
        inp = inp + "Y" + "\n" + output_file + "\n"
    else:
        inp = inp + "N" + "\n"
    if supress_g:
        inp = inp + "Y" + "\n"
    else:
        inp = inp + "N" + "\n"
    if supress_ic:
        inp = inp + "Y" + "\n"
    else:
        inp = inp + "N" + "\n" + dcc_theory + "\n"
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def hsicc(inputdict_unchecked):
    """
    This program calculates internal conversion coefficients. (HSICC readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            data_deck : string, data deck to be used for hsicc program.
            icc_index : string, icc index to be used for hsicc program.
            icc_table : string, icc table to be used for the hsicc program.
            complete_report : string, desired report file path for hsicc
                              program.
            new_card_deck : string, desired new card deck file path for hsicc
                            program.
            comparison_report : string, desired comparison report path for
                                hsicc program.
            is_multipol_known : int, 1 if multipol is known, 0 otherwise.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if HSICC completes
        successfully.
    """
    inputdict = {}
    data_deck = inputdict_unchecked["data_deck"]
    icc_index = inputdict_unchecked["icc_index"]
    icc_table = inputdict_unchecked["icc_table"]
    complete_report = inputdict_unchecked["complete_report"]
    new_card_deck = inputdict_unchecked["new_card_deck"]
    comparison_report = inputdict_unchecked["comparison_report"]
    multipol_known = inputdict_unchecked["is_multipol_known"]  #'Y or CR'

    exe_path = path_to_exe("hsicc")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = (
        data_deck
        + "\n"
        + icc_index
        + "\n"
        + icc_table
        + "\n"
        + complete_report
        + "\n"
        + new_card_deck
        + "\n"
        + comparison_report
        + "\n"
        + multipol_known
    )
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def hsmrg(inputdict_unchecked):
    """
    This program merges new gamma records created by HSICC with the original
    input data.  (HSICC readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            data_deck : string, data deck file path for hsmrg to use.
            card_deck : string, card deck file path for hsmrg to use.
            merged_data_deck : string, desired merged data deck file path
                               created by hsmrg.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if HSMRG completes
        successfully.
    """
    inputdict = {}
    data_deck = inputdict_unchecked["data_deck"]
    card_deck = inputdict_unchecked["card_deck"]
    merged_data_deck = inputdict_unchecked["merged_data_deck"]

    exe_path = path_to_exe("hsmrg")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = data_deck + "\n" + card_deck + "\n" + merged_data_deck
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def seqhst(inputdict_unchecked):
    """
    This program recreates a sequential file of the internal conversion table
    from the direct access file.  (HSICC readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            binary_table_input_file : string, binary table input file path.
            sequential_output_file : string, desired path of sequential output
                                     file.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if SEQHST completes
        successfully.
    """
    # NOTE: changed input file line length to 90 to support longer file paths
    #      in fortran source.
    inputdict = {}
    input_file = inputdict_unchecked["binary_table_input_file"]
    output_file = inputdict_unchecked["sequential_output_file"]

    exe_path = path_to_exe("seqhst")
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    inp = input_file + "\n" + output_file
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def logft(inputdict_unchecked):
    # NOTE: changed input file line length to 90 to support longer file paths
    #      in fortran source.
    """
    This program calculates log ft values for beta and electron-capture decay,
    average beta energies, and capture fractions.  (LOGFT readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_data_set : string, path to input data file.
            output_report : string, desired path to output report file.
            data_table : string, path to data table.
            output_data_set : string, desired path to output data set.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if LOGFT completes
        successfully.
    """
    inputdict = {}
    input_data_set = inputdict_unchecked["input_data_set"]
    output_report = inputdict_unchecked["output_report"]
    data_table = inputdict_unchecked["data_table"]
    output_data_set = inputdict_unchecked["output_data_set"]

    exe_path = path_to_exe("logft")
    inp = (
        input_data_set
        + "\n"
        + output_report
        + "\n"
        + data_table
        + "\n"
        + output_data_set
        + "\n"
    )
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    proc.stdin.write(inp.encode("utf-8"))
    proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def radd(inputdict_unchecked):
    """
    This code (RadD.FOR) deduces the radius parameter (r 0 ) for odd-odd and
    odd-A nuclei using the even-even radii [1] as input parameters.

    These radii deduced for odd-A and odd-odd nuclides can be used in the
    calculation of alpha hindrance factors. In this procedure, it is assumed
    that radius parameter ( r 0 Z , N ) for odd-Z and odd-N nuclides lies midway
    between the radius parameters of adjacent even-even neighbors calculates
    reduced transition probabilities. (RADD readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, input ensdf file
            output file : string, file for output to be written to (doesn't
                          have to exist)

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if RADD completes
        successfully.
    """
    inputdict = {}
    atomic_number = inputdict_unchecked["atomic_number"]
    neutron_number = inputdict_unchecked["neutron_number"]
    output_file = inputdict_unchecked["output_file"]

    # Create symlinks to the two binaries the radd executables uses.
    ak04_path = path_to_exe("98AK04.in")
    ele_path = path_to_exe("ELE.in")
    ak04_set = False
    ele_set = False
    if not os.path.exists("98AK04.in"):
        os.symlink(ak04_path, "98AK04.in")
        ak04_set = True
    if not os.path.exists("ELE.in"):
        os.symlink(ele_path, "ELE.in")
        ele_set = True

    exe_path = path_to_exe("radd")
    inp = atomic_number + "\n" + neutron_number + "\n" + "NO" + "\n"
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    proc.stdin.write(inp.encode("utf-8"))
    radd_output = proc.communicate()[0]
    proc.stdin.close()
    with open(output_file, "w") as f:
        f.write(radd_output.decode("utf-8"))

    if ak04_set:
        os.remove("98AK04.in")
    if ele_set:
        os.remove("ELE.in")
    return inputdict_unchecked


def radlist(inputdict_unchecked):
    """
    This program calculates atomic & nuclear radiations and checks energy
    balance. (RADLIST readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            output_radiation_listing : string, 'Y' if output radiation listing
                                       is desired, else 'N'.
            output_ensdf_like_file : string, 'Y' if output ensdf like file is
                                     desired, else 'N'.
            output_file_for_nudat : string, 'Y' if output file for nudat is desired,
                                    else 'N'.
            output_mird_listing : string, 'Y' if output mird listing is desired,
                                  else 'N'.
            calculate_continua : string, 'Y' if calculate continua is desired,
                                 else 'N'.
            input_file : string, input ensdf file.
            output_radlst_file : string, path to desired output radlst file.
            input_radlst_data_table : string, path to input radlst data table
                                      (mednew.dat location).
            input_masses_data_table : string, (optional) path to input masses
                                      data table.
            output_ensdf_file : string, path to desired output ensdf file.

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if RADLIST completes
        successfully.
    """
    exe_path = path_to_exe("radlist")
    print(exe_path)
    verify_download_exe(exe_path)

    inputdict = {}
    output_rad_listing = inputdict_unchecked["output_radiation_listing"]
    output_endf_like_file = inputdict_unchecked["output_ensdf_like_file"]
    output_file_for_nudat = inputdict_unchecked["output_file_for_nudat"]
    output_mird_listing = inputdict_unchecked["output_mird_listing"]
    calculate_continua = inputdict_unchecked["calculate_continua"]
    input_file = inputdict_unchecked["input_file"]
    output_radlst_file = inputdict_unchecked["output_radlst_file"]
    input_radlst_data_table = inputdict_unchecked["input_radlst_data_table"]
    if "input_masses_data_table" in inputdict_unchecked:
        input_masses_data_table = inputdict_unchecked["input_masses_data_table"]
    else:
        input_masses_data_table = ""
    output_ensdf_file = inputdict_unchecked["output_ensdf_file"]
    inp = (
        output_rad_listing
        + "\n"
        + output_endf_like_file
        + "\n"
        + output_file_for_nudat
        + "\n"
        + output_mird_listing
        + "\n"
        + calculate_continua
        + "\n"
        + input_file
        + "\n"
        + output_radlst_file
        + "\n"
        + input_radlst_data_table
        + "\n"
        + input_masses_data_table
        + "\n"
        + output_ensdf_file
    )
    proc = subprocess.Popen([exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    proc.stdin.write(inp.encode("utf-8"))
    radd_output = proc.communicate()[0]
    proc.stdin.close()
    return inputdict_unchecked


def ruler(inputdict_unchecked):
    """
    This program calculates reduced transition probabilities. (RULER readme)

    Parameters
    ----------
    inputdict_unchecked : dictionary
        dictionary that must have the following key-pair values:
            input_file : string, input ensdf file
            output file : string, file for output to be written to (doesn't
                          have to exist)

    Returns
    -------
    rtn : dictionary
        Everything in input dictionary is returned if RULER completes
        successfully.
    """
    inputdict = {}
    input_file = inputdict_unchecked["input_file"]
    output_report_file = inputdict_unchecked["output_report_file"]
    mode_of_operation = inputdict_unchecked["mode_of_operation"]
    assumed_dcc_theory = inputdict_unchecked["assumed_dcc_theory"]

    exe_path = path_to_exe("ruler")
    ruler_output = subprocess.Popen(
        [exe_path], stdout=subprocess.PIPE, stdin=subprocess.PIPE
    )
    inp = (
        input_file
        + "\n"
        + output_report_file
        + "\n"
        + mode_of_operation
        + "\n"
        + assumed_dcc_theory
    )
    ruler_output.stdin.write(inp.encode("utf-8"))
    ruler_output.communicate()[0]
    ruler_output.stdin.close()
    return inputdict_unchecked
