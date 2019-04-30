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
    a = os.system('which h5diff')
    if a != 0:
        raise SkipTest

    for response in responses:
        # read  output.txt and write h5 file
        filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '_output.txt']))
        h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '_output.txt.h5']))
        response_to_hdf5(filename, response)

        # generate expected h5 file
        exp_h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join(['exp_', response, '_output.txt.h5']))
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
        h5_filename = os.path.join(thisdir, "files_test_activation_responses", ''.join([response, '_output.txt.h5']))
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

