import os
import warnings
from nose.tools import assert_equal, assert_almost_equal, assert_true
from nose.plugins.skip import SkipTest
import numpy as np
from numpy.testing import assert_array_equal
import multiprocessing
import filecmp
import sys
from shutil import copyfile
import tables as tb

from pyne.mcnp import Meshtal
from pyne.material import Material
from pyne.utils import QAWarning, file_block_almost_same
from pyne.alara import response_to_hdf5, response_hdf5_to_mesh, _make_response_dtype
from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB

if not HAVE_PYMOAB:
    raise SkipTest

if sys.version_info[0] > 2:
    izip = zip
else:
    from itertools import izip

warnings.simplefilter("ignore", QAWarning)

thisdir = os.path.dirname(__file__)

responses = [
    "decay_heat",
    "specific_activity",
    "alpha_heat",
    "beta_heat",
    "gamma_heat",
    "wdr",
    "photon_source",
]


def _generate_exp_h5(filename, response, exp_h5_filename):
    """
    This function is used to generate expected h5 file for different responses.
    Supported responses is defined at the begining of this file.
    Filename could can be an output.txt contains specific response for a file
    contains multiple responses.
    """
    # generate expected h5 file
    f = open(filename, "r")
    f.seek(0)
    dt = _make_response_dtype(response)
    filters = tb.Filters(complevel=1, complib="zlib")
    h5f = tb.open_file(exp_h5_filename, "w", filters=filters)
    tab = h5f.create_table("/", "data", dt)
    rows = np.empty(12, dtype=dt)
    rows[0] = (0, "h-3", "shutdown", 9.5258e-18)
    rows[1] = (0, "h-3", "1000 s", 9.5258e-18)
    rows[2] = (0, "h-3", "12 h", 9.5251e-18)
    rows[3] = (0, "h-3", "3 d", 9.5214e-18)
    rows[4] = (0, "n-16", "shutdown", 3.1588e-9)
    rows[5] = (0, "n-16", "1000 s", 0.0000e0)
    rows[6] = (0, "n-16", "12 h", 0.0000e0)
    rows[7] = (0, "n-16", "3 d", 0.0000e0)
    rows[8] = (0, "TOTAL", "shutdown", 3.1588e-9)
    rows[9] = (0, "TOTAL", "1000 s", 9.5258e-18)
    rows[10] = (0, "TOTAL", "12 h", 9.5251e-18)
    rows[11] = (0, "TOTAL", "3 d", 9.5214e-18)
    tab.append(rows[:])
    # close the file
    h5f.close()
    f.close()


def test_response_to_hdf5():
    """
    This function test alara.response_to_hdf5, with response of:
        - decay_heat
        - specific_activity
        - alpha_heat
        - beta_heat
        - gamma_heat
        - wdr
        - photon_source
    """
    # skip test if h5diff not exist
    is_h5diff = os.system("which h5diff")
    if is_h5diff != 0:
        raise SkipTest

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(
            thisdir,
            "files_test_activation_responses",
            "".join([response, "_output.txt"]),
        )
        h5_filename = os.path.join(
            thisdir, "files_test_activation_responses", "".join([response, ".h5"])
        )
        response_to_hdf5(filename, response)

        # generate expected h5 file
        exp_h5_filename = os.path.join(
            thisdir,
            "files_test_activation_responses",
            "".join(["exp_", response, ".h5"]),
        )
        _generate_exp_h5(filename, response, exp_h5_filename)

        # compare two h5 files
        command = "".join(["h5diff ", h5_filename, " ", exp_h5_filename])
        diff_flag = os.system(command)
        # return value 0 if no difference, 1 if differences found, 2 if error
        assert_equal(diff_flag, 0)

        # remove generated files
        os.remove(h5_filename)
        os.remove(exp_h5_filename)


def test_responses_to_hdf5_multiple():
    """
    This function test alara.response_to_hdf5, read output.txt of multiple responses:
        - decay_heat
        - specific_activity
        - alpha_heat
        - beta_heat
        - gamma_heat
        - wdr
        - photon_source
    """
    # skip test if h5diff not exist
    is_h5diff = os.system("which h5diff")
    if is_h5diff != 0:
        raise SkipTest

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(
            thisdir, "files_test_activation_responses", "multiple_output.txt"
        )
        h5_filename = os.path.join(
            thisdir, "files_test_activation_responses", "".join([response, ".h5"])
        )
        response_to_hdf5(filename, response)

        # generate expected h5 file
        exp_h5_filename = os.path.join(
            thisdir,
            "files_test_activation_responses",
            "".join(["exp_", response, ".h5"]),
        )
        _generate_exp_h5(filename, response, exp_h5_filename)

        # compare two h5 files
        command = "".join(["h5diff ", h5_filename, " ", exp_h5_filename])
        diff_flag = os.system(command)
        # return value 0 if no difference, 1 if differences found, 2 if error
        assert_equal(diff_flag, 0)

        # remove generated files
        os.remove(h5_filename)
        os.remove(exp_h5_filename)


def test_response_hdf5_to_mesh():
    """Tests the function photon source_h5_to_mesh."""

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(
            thisdir,
            "files_test_activation_responses",
            "".join([response, "_output.txt"]),
        )
        h5_filename = os.path.join(
            thisdir, "files_test_activation_responses", "".join([response, ".h5"])
        )
        response_to_hdf5(filename, response)
        assert_true(os.path.exists(h5_filename))

        mesh = Mesh(structured=True, structured_coords=[[0, 1], [0, 1], [0, 1]])

        tags = {("h-3", "shutdown"): "tag1", ("TOTAL", "12 h"): "tag2"}
        response_hdf5_to_mesh(mesh, h5_filename, tags, response)

        # create lists of lists of expected results
        tag1_answers = [9.5258e-18]
        tag2_answers = [9.5251e-18]

        ves = list(mesh.structured_iterate_hex("xyz"))
        for i, ve in enumerate(ves):
            assert_equal(mesh.tag1[ve], tag1_answers[i])
            assert_equal(mesh.tag2[ve], tag2_answers[i])

        if os.path.isfile(h5_filename):
            os.remove(h5_filename)


def _activation_responses_test_step1(activation_responses_run_dir):
    os.chdir(thisdir)
    # copy ../scripts/activation_responses.py to activation_responses_run_dir/activation_responses.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(activation_responses_run_dir, "activation_responses.py")
    copyfile(os.path.join(folderpath, "scripts", "activation_responses.py"), dst)

    # run activation_responses step1
    os.chdir(activation_responses_run_dir)
    os.system("python activation_responses.py step1")

    # output files of activation_responses step1
    alara_inp = os.path.join(activation_responses_run_dir, "alara_inp")
    alara_matlib = os.path.join(activation_responses_run_dir, "alara_matlib")
    alara_fluxin = os.path.join(activation_responses_run_dir, "alara_fluxin")
    blank_mesh = os.path.join(activation_responses_run_dir, "blank_mesh.h5m")
    step1_file = os.path.join(
        activation_responses_run_dir, "activation_responses_step1.h5m"
    )

    exp_alara_inp = os.path.join(activation_responses_run_dir, "exp_alara_inp")
    exp_alara_matlib = os.path.join(activation_responses_run_dir, "exp_alara_matlib")
    exp_alara_fluxin = os.path.join(activation_responses_run_dir, "exp_alara_fluxin")

    # compare the output file of step1
    f1 = filecmp.cmp(alara_inp, exp_alara_inp)
    f2 = file_block_almost_same(alara_matlib, exp_alara_matlib)
    f3 = filecmp.cmp(alara_fluxin, exp_alara_fluxin)

    # remove test output files
    os.remove(alara_inp)
    os.remove(alara_fluxin)
    os.remove(alara_matlib)
    os.remove(blank_mesh)
    os.remove(step1_file)
    os.remove(dst)

    assert_equal(f1, True)
    assert_equal(f2, True)
    assert_equal(f3, True)


def _activation_responses_test_step2(activation_responses_run_dir):
    # skip test if h5diff not exist
    is_h5diff = os.system("which h5diff")
    if is_h5diff != 0:
        raise SkipTest

    os.chdir(thisdir)
    # copy ../scripts/activation_responses.py to activation_responses_run_dir/activation_responses.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(activation_responses_run_dir, "activation_responses.py")
    copyfile(os.path.join(folderpath, "scripts", "activation_responses.py"), dst)

    # output files of activation_responses step1
    alara_inp = os.path.join(activation_responses_run_dir, "alara_inp")
    copyfile(os.path.join(activation_responses_run_dir, "exp_alara_inp"), alara_inp)
    blank_mesh = os.path.join(activation_responses_run_dir, "blank_mesh.h5m")
    copyfile(
        os.path.join(activation_responses_run_dir, "exp_blank_mesh.h5m"), blank_mesh
    )

    # run activation_responses step2
    os.chdir(activation_responses_run_dir)
    os.system("python activation_responses.py step2")

    response = "decay_heat"
    # output files of activation_responses step2
    h5_filename = os.path.join(activation_responses_run_dir, "".join([response, ".h5"]))
    exp_h5_filename = os.path.join(
        activation_responses_run_dir, "".join(["exp_", response, ".h5"])
    )
    h5m_filename = os.path.join(
        activation_responses_run_dir, "".join([response, "_1.h5m"])
    )
    exp_h5m_filename = os.path.join(
        activation_responses_run_dir, "".join(["exp_", response, "_1.h5m"])
    )

    # compare the results
    # compare two h5 files
    command = "".join(["h5diff ", h5_filename, " ", exp_h5_filename])
    diff_flag4 = os.system(command)
    # compare two h5m files
    command = "".join(
        [
            "h5diff ",
            h5m_filename,
            " ",
            exp_h5m_filename,
            "".join([" tstt/tags/", response, " /tstt/tags/", response]),
        ]
    )

    diff_flag5 = os.system(command)

    # remove test generated files
    os.remove(blank_mesh)
    os.remove(alara_inp)
    os.remove(h5_filename)
    os.remove(h5m_filename)
    os.remove(dst)

    # return value 0 if no difference, 1 if differences found, 2 if error
    assert_equal(diff_flag4, 0)
    assert_equal(diff_flag5, 0)


def test_activation_responses_script():
    # skip test without dagmc
    try:
        from pyne import dagmc
    except ImportError:
        raise SkipTest

    activation_responses_run_dir = os.path.join(
        thisdir, "files_test_activation_responses", "activation_responses_examples"
    )
    _activation_responses_test_step1(activation_responses_run_dir)
    _activation_responses_test_step2(activation_responses_run_dir)
