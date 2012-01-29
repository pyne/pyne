import os

import numpy as np
import tables as tb

from nose.tools import assert_equal, assert_not_equal, assert_almost_equal, assert_true, \
                       assert_raises
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.xs.cache import xs_cache
from pyne.xs.channels import sigma_f
from pyne.pyne_config import pyne_conf



#
# Test channels, make sure they run rather than test their values.
# This is OK since the underlying functions are very well tested.
#

def test_sigma_f1():
    E_g = np.array([10.0, 7.5, 5.0, 2.5, 0.1])
    E_n = xs_cache['E_n']
    phi_n = np.ones(len(E_n) - 1)

    sig_f = sigma_f('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f('U238', E_g, E_n, phi_n)
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

    sig_f = sigma_f('U235')
    observed = (0.0 <= sig_f).all()
    assert_true(observed)

