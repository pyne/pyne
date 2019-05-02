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
from pyne.utils import QAWarning
from pyne.alara import response_to_hdf5, response_hdf5_to_mesh
from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB
if not HAVE_PYMOAB:
    raise SkipTest

if sys.version_info[0] > 2:
    izip = zip
else:
    from itertools import izip

warnings.simplefilter("ignore", QAWarning)

thisdir = os.path.dirname(__file__)

responses = ['decay_heat', 'specific_activity', 'alpha_heat', 'beta_heat',
             'gamma_heat', 'wdr', 'photon_source']
#responses = ['specific_activity']

def _generate_exp_h5(filename, response, exp_h5_filename):
    """
    This function is used to generate expected h5 file for different responses.
    Supported responses is defined at the begining of this file.
    Filename could can be an output.txt contains specific response for a file
    contains multiple responses.
    """
    # generate expected h5 file
    f = open(filename, 'r')
    f.seek(0)
    dt = np.dtype([
        ('idx', np.int64),
        ('nuc', 'S6'),
        ('time', 'S20'),
        (response, np.float64)])
    filters = tb.Filters(complevel=1, complib='zlib')
    h5f = tb.open_file(exp_h5_filename, 'w', filters=filters)
    tab = h5f.create_table('/', 'data', dt)
    rows = np.empty(12, dtype=dt)
    rows[0] = (0, 'h-3', 'shutdown', 9.5258e-18)
    rows[1] = (0, 'h-3', '1000 s', 9.5258e-18)
    rows[2] = (0, 'h-3', '12 h', 9.5251e-18)
    rows[3] = (0, 'h-3', '3 d', 9.5214e-18)
    rows[4] = (0, 'n-16', 'shutdown', 3.1588e-9)
    rows[5] = (0, 'n-16', '1000 s', 0.0000e0)
    rows[6] = (0, 'n-16', '12 h', 0.0000e0)
    rows[7] = (0, 'n-16', '3 d', 0.0000e0)
    rows[8] = (0, 'TOTAL', 'shutdown', 3.1588e-9)
    rows[9] = (0, 'TOTAL', '1000 s', 9.5258e-18)
    rows[10] = (0, 'TOTAL', '12 h', 9.5251e-18)
    rows[11] = (0, 'TOTAL', '3 d', 9.5214e-18)
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
    is_h5diff = os.system('which h5diff')
    if is_h5diff != 0:
        raise SkipTest

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '_output.txt']))
        h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '.h5']))
        response_to_hdf5(filename, response)

        # generate expected h5 file
        exp_h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join(['exp_', response, '.h5']))
        _generate_exp_h5(filename, response, exp_h5_filename)

        # compare two h5 files
        command = ''.join(['h5diff ', h5_filename, ' ', exp_h5_filename])
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
    is_h5diff = os.system('which h5diff')
    if is_h5diff != 0:
        raise SkipTest

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(thisdir, "files_test_activation_responses", 'multiple_output.txt')
        h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '.h5']))
        response_to_hdf5(filename, response)

        # generate expected h5 file
        exp_h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join(['exp_', response, '.h5']))
        _generate_exp_h5(filename, response, exp_h5_filename)

        # compare two h5 files
        command = ''.join(['h5diff ', h5_filename, ' ', exp_h5_filename])
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
        filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '_output.txt']))
        h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '.h5']))
        response_to_hdf5(filename, response)
        assert_true(os.path.exists(h5_filename))
    
        mesh = Mesh(structured=True,
                    structured_coords=[[0, 1], [0, 1], [0, 1]])
    
        tags = {('h-3', 'shutdown'): 'tag1', ('TOTAL', '12 h'): 'tag2'}
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


def _r2s_test_step1(r2s_run_dir):
    os.chdir(thisdir)
    # copy ../scripts/r2s.py to r2s_run_dir/r2s.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(r2s_run_dir, "r2s.py")
    copyfile(os.path.join(folderpath, "scripts", "r2s.py"), dst)

    # run r2s step1
    os.chdir(r2s_run_dir)
    os.system('python r2s.py step1')

    # output files of r2s step1
    alara_inp = os.path.join(r2s_run_dir, "alara_inp")
    alara_matlib = os.path.join(r2s_run_dir, "alara_matlib")
    alara_fluxin = os.path.join(r2s_run_dir, "alara_fluxin")
    blank_mesh = os.path.join(r2s_run_dir, "blank_mesh.h5m")
    step1_file = os.path.join(r2s_run_dir, "r2s_step1.h5m")

    exp_alara_inp = os.path.join(r2s_run_dir, "exp_alara_inp")
    exp_alara_matlib = os.path.join(r2s_run_dir, "exp_alara_matlib")
    exp_alara_fluxin = os.path.join(r2s_run_dir, "exp_alara_fluxin")

    # compare the output file of step1
    f1 = filecmp.cmp(alara_inp, exp_alara_inp)
    f2 = filecmp.cmp(alara_matlib, exp_alara_matlib)
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


def _r2s_test_step2(r2s_run_dir):
    os.chdir(thisdir)
    # copy ../scripts/r2s.py to r2s_run_dir/r2s.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(r2s_run_dir, "r2s.py")
    copyfile(os.path.join(folderpath, "scripts", "r2s.py"), dst)

    # output files of r2s step1
    alara_inp = os.path.join(r2s_run_dir, "alara_inp")
    copyfile(os.path.join(r2s_run_dir, "exp_alara_inp"), alara_inp)
    blank_mesh = os.path.join(r2s_run_dir, "blank_mesh.h5m")
    copyfile(os.path.join(r2s_run_dir, "exp_blank_mesh.h5m"), blank_mesh)

    # run r2s step2
    os.chdir(r2s_run_dir)
    os.system('python r2s.py step2')

    # output files of r2s step2
    e_bounds = os.path.join(r2s_run_dir, "e_bounds")
    p_src = os.path.join(r2s_run_dir, "phtn_src.h5")
    t_p_src = os.path.join(r2s_run_dir, "total_photon_source_intensities.txt")
    src_c1 = os.path.join(r2s_run_dir, "source_1.h5m")

    exp_e_bounds = os.path.join(r2s_run_dir, "exp_e_bounds")
    exp_t_p_src = os.path.join(
        r2s_run_dir, "exp_total_photon_source_intensities.txt")

    # compare the results
    f4 = filecmp.cmp(e_bounds, exp_e_bounds)
    f5 = filecmp.cmp(t_p_src, exp_t_p_src)

    # remove test generated files
    os.remove(blank_mesh)
    os.remove(alara_inp)
    os.remove(e_bounds)
    os.remove(p_src)
    os.remove(t_p_src)
    os.remove(src_c1)
    os.remove(dst)

    assert_equal(f4, True)
    assert_equal(f5, True)


def test_r2s_script():

    # skip test without dagmc
    try:
        from pyne import dagmc
    except ImportError:
        raise SkipTest

    r2s_run_dir = os.path.join(
        thisdir, "files_test_r2s", "r2s_examples", "r2s_run")
    _r2s_test_step1(r2s_run_dir)
    _r2s_test_step2(r2s_run_dir)

