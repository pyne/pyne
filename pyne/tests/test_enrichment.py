"""Enrichment Tests"""
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import numpy as np

import pyne
from pyne.material import Material
from pyne import enrichment as enr

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


#
# Tests the cascade helper class.
#

def test_cascade_constructor():
    casc = enr.Cascade()
    assert_almost_equal(casc.alpha, 0.0)
    assert_almost_equal(casc.Mstar, 0.0)
    assert_equal(casc.j, 0)
    assert_equal(casc.k, 0)
    assert_almost_equal(casc.N, 0.0)
    assert_almost_equal(casc.M, 0.0)
    assert_almost_equal(casc.x_feed_j, 0.0)
    assert_almost_equal(casc.x_prod_j, 0.0)
    assert_almost_equal(casc.x_tail_j, 0.0)
    assert_equal(casc.mat_feed, Material())
    assert_equal(casc.mat_prod, Material())
    assert_equal(casc.mat_tail, Material())
    assert_almost_equal(casc.l_t_per_feed, 0.0)
    assert_almost_equal(casc.swu_per_feed, 0.0)
    assert_almost_equal(casc.swu_per_prod, 0.0)

def test_alpha():
    casc = enr.Cascade()
    casc.alpha = 1.05
    assert_equal(casc.alpha, 1.05)

def test_Mstar():
    casc = enr.Cascade()
    casc.Mstar = 236.5
    assert_equal(casc.Mstar, 236.5)

def test_j():
    casc = enr.Cascade()
    casc.j = 922350
    assert_equal(casc.j, 922350)

def test_k():
    casc = enr.Cascade()
    casc.k = 922380
    assert_equal(casc.k, 922380)

def test_N():
    casc = enr.Cascade()
    casc.N = 30.0
    assert_equal(casc.N, 30.0)

def test_M():
    casc = enr.Cascade()
    casc.M = 10.0
    assert_equal(casc.M, 10.0)

def test_x_feed_j():
    casc = enr.Cascade(x_feed_j=0.0072)
    assert_equal(casc.x_feed_j, 0.0072)

def test_x_prod_j():
    casc = enr.Cascade()
    casc.x_prod_j = 0.05
    assert_equal(casc.x_prod_j, 0.05)

def test_x_tail_j():
    casc = enr.Cascade()
    casc.x_tail_j = 0.0025
    assert_equal(casc.x_tail_j, 0.0025)

def test_default_uranium_cascade():
    casc = enr.default_uranium_cascade()
    assert_equal(casc.alpha, 1.05)
    assert_equal(casc.Mstar, 236.5)
    assert_equal(casc.j, 922350)
    assert_equal(casc.k, 922380)
    assert_equal(casc.N, 30.0)
    assert_equal(casc.M, 10.0)
    assert_equal(casc.x_feed_j, 0.0072)
    assert_equal(casc.x_prod_j, 0.05)
    assert_equal(casc.x_tail_j, 0.0025)
    assert_equal(casc.mat_feed, Material({922340: 5.5e-05, 922350: 0.0072, 
                                922380: 0.992745}, 1.0, 1.0))

def test_prod_per_feed():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    exp = (xf - xt) / (xp - xt)
    obs = enr.prod_per_feed(xf, xp, xt)
    assert_almost_equal(obs, exp)

def test_tail_per_feed():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    exp = (xf - xp) / (xt - xp)
    obs = enr.tail_per_feed(xf, xp, xt)
    assert_almost_equal(obs, exp)

def test_tail_per_prod():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    exp = (xf - xp) / (xt - xf)
    obs = enr.tail_per_prod(xf, xp, xt)
    assert_almost_equal(obs, exp)

def test_alphastar_i():
    a, ms, mi = 1.05, 236.5, 235.0
    exp = a**(ms - mi)
    obs = enr.alphastar_i(a, ms, mi)
    assert_almost_equal(obs, exp)


#
# Intergration tests which test multicomponent() and ltot_per_feed()
#

def test_sample_feed():
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.06
    feed = Material({
            922320: 1.1 * (10.0**-9),
            922340: 0.00021,
            922350: 0.0092,
            922360: 0.0042,
            922380: 0.9863899989,
            })
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, tolerance=1E-11, max_iter=100)

    assert_almost_equal(casc.mat_prod.comp[922350], 0.06,   5) 
    assert_almost_equal(casc.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0,                 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.11652173913043479, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.88347826086956527, 1.0)

    assert_almost_equal(casc.N / 26.8646352802, 1.0, 5)
    assert_almost_equal(casc.M / 16.6379009423, 1.0, 5)

    assert_almost_equal(casc.Mstar / 236.577085, 1.0, 5)

    assert_almost_equal(casc.l_t_per_feed / 357.388791749,  1.0, 5)
    assert_almost_equal(casc.swu_per_feed / 0.932280175218, 1.0, 5)
    assert_almost_equal(casc.swu_per_prod / 8.0009119515,   1.0, 5)

def test_NU():
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.05
    feed = Material({
            922340: 0.000055,
            922350: 0.00720,
            922380: 0.992745,
            })
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, tolerance=1E-11)

    assert_almost_equal(casc.mat_prod.comp[922350], 0.05,   5) 
    assert_almost_equal(casc.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0,             1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.0989473684211, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.901052631579,  1.0)

    assert_almost_equal(casc.N / 27.1835088212, 1.0, 4)
    assert_almost_equal(casc.M / 13.3875092512, 1.0, 4)

    assert_almost_equal(casc.Mstar / 236.562179, 1.0, 5)

    assert_almost_equal(casc.l_t_per_feed / 288.627270162,  1.0, 5)
    assert_almost_equal(casc.swu_per_feed / 0.761263453429, 1.0, 5)
    assert_almost_equal(casc.swu_per_prod / 7.69362000806,  1.0, 5)

def test_vision():
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.055
    feed = Material({
            922340: 0.000183963025893197,
            922350: 0.00818576605617839,
            922360: 0.00610641667100979,
            922380: 0.985523854246919,
            })
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, tolerance=1E-11)

    assert_almost_equal(casc.mat_prod.comp[922350], 0.055,  5) 
    assert_almost_equal(casc.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0,                 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.10830030583196934, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.89169969416803063, 1.0)

    assert_almost_equal(casc.N / 27.38162850698868, 1.0, 2)
    assert_almost_equal(casc.M / 15.09646512546496, 1.0, 2)

    assert_almost_equal(casc.Mstar / 236.581784, 1.0, 4)

    assert_almost_equal(casc.l_t_per_feed / 326.895568684,  1.0, 4)
    assert_almost_equal(casc.swu_per_feed / 0.85102089049,  1.0, 4)
    assert_almost_equal(casc.swu_per_prod / 7.85797310499,  1.0, 4)

def test_tungsten():
    # This test comes from 'Multicomponent Isotope Separation in Matched
    # Abundance Ratio Cascades Composed of Stages with Large Separation Factors' 
    # by casc. von Halle, 1987.
    orig_casc = enr.Cascade()
    orig_casc.alpha = 1.16306
    orig_casc.Mstar = 181.3
    orig_casc.j = 741800
    orig_casc.k = 741860
    orig_casc.N = 30.0
    orig_casc.M = 10.0
    orig_casc.x_prod_j = 0.5109
    orig_casc.x_tail_j = 0.00014

    feed = Material({
            741800: 0.0014, 
            741820: 0.26416, 
            741830: 0.14409, 
            741840: 0.30618, 
            741860: 0.28417,
            })
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, tolerance=1E-5)

    assert_almost_equal(casc.mat_prod.comp[741800], 0.5109,  5) 
    assert_almost_equal(casc.mat_tail.comp[741800], 0.00014, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0,                   1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.0024669120526274574, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.99753308794737272,   1.0)

    assert_almost_equal(casc.N / 43.557515688533513, 1.0, 2)
    assert_almost_equal(casc.M / 11.49556481009056,  1.0, 2)

    assert_almost_equal(casc.Mstar / 181.164592, 1.0, 4)

    assert_almost_equal(casc.l_t_per_feed / 96.8179316719, 1.0, 3)
    assert_almost_equal(casc.swu_per_feed / 2.22221945305, 1.0, 3)
    assert_almost_equal(casc.swu_per_prod / 900.810164953, 1.0, 3)

if __name__ == "__main__":
    noscasc.main()

