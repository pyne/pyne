import os
import warnings
from nose.tools import assert_equal, assert_almost_equal
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
from pyne.alara import response_to_hdf5
from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB
if not HAVE_PYMOAB:
    raise SkipTest

if sys.version_info[0] > 2:
    izip = zip
else:
    from itertools import izip

warnings.simplefilter("ignore", QAWarning)

thisdir = os.path.dirname(__file__)


def test_response_to_hdf5_decay_heat():
    """
    This function test alara.response_to_hdf5, with response of decay_heat.
    """
    response = 'decay_heat'
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
