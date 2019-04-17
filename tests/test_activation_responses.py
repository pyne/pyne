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

from pyne.mcnp import Meshtal
from pyne.material import Material
from pyne.r2s import irradiation_setup, photon_sampling_setup, total_photon_source_intensity
from pyne.utils import QAWarning
from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB
if not HAVE_PYMOAB:
    raise SkipTest

if sys.version_info[0] > 2:
    izip = zip
else:
    from itertools import izip

warnings.simplefilter("ignore", QAWarning)

thisdir = os.path.dirname(__file__)


