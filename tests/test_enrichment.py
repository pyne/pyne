"""Enrichment Tests"""
import nose
import warnings

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_almost_equal,
    assert_true,
    assert_false,
    with_setup,
)

import os
import warnings
import numpy as np
import math

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

import pyne
from pyne.material import Material
from pyne import enrichment as enr

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")

SOLVERS = ["symbolic", "numeric"]


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
    casc.j = 922350000
    assert_equal(casc.j, 922350000)


def test_k():
    casc = enr.Cascade()
    casc.k = 922380000
    assert_equal(casc.k, 922380000)


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
    assert_equal(casc.j, 922350000)
    assert_equal(casc.k, 922380000)
    assert_equal(casc.N, 30.0)
    assert_equal(casc.M, 10.0)
    assert_equal(casc.x_feed_j, 0.0072)
    assert_equal(casc.x_prod_j, 0.05)
    assert_equal(casc.x_tail_j, 0.0025)
    assert_equal(
        casc.mat_feed,
        Material(
            {922340000: 5.5e-05, 922350000: 0.0072, 922380000: 0.992745}, 1.0, 1.0
        ),
    )


def test_prod():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    feed, prod, tails = 15.1596, 1.5, 13.6596
    swu = 11765.0

    exp = prod
    obs = enr.product(xf, xp, xt, feed=feed)
    assert_almost_equal(obs, exp, places=4)
    obs = enr.product(xf, xp, xt, tails=tails)
    assert_almost_equal(obs, exp, places=4)


def test_feed():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    feed, prod, tails = 15.1596, 1.5, 13.6596
    swu = 11765.0

    exp = feed
    obs = enr.feed(xf, xp, xt, product=prod)
    assert_almost_equal(obs, exp, places=4)
    obs = enr.feed(xf, xp, xt, tails=tails)
    assert_almost_equal(obs, exp, places=4)


def test_tails():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    feed, prod, tails = 15.1596, 1.5, 13.6596
    swu = 11765.0

    exp = tails
    obs = enr.tails(xf, xp, xt, feed=feed)
    assert_almost_equal(obs, exp, places=4)
    obs = enr.tails(xf, xp, xt, product=prod)
    assert_almost_equal(obs, exp, places=4)


def test_value():
    x = 0.0072
    exp = (2 * x - 1) * math.log(x / (1 - x))
    obs = enr.value_func(x)
    assert_almost_equal(exp, obs)


def test_swu():
    xf, xp, xt = 0.0072, 0.05, 0.0025
    feed, prod, tails = 15.1596, 1.5, 13.6596
    swu = 11765.0 / 1e3

    exp = swu
    obs = enr.swu(xf, xp, xt, feed=feed)
    assert_almost_equal(obs, exp, places=3)
    obs = enr.swu(xf, xp, xt, product=prod)
    assert_almost_equal(obs, exp, places=3)
    obs = enr.swu(xf, xp, xt, tails=tails)
    assert_almost_equal(obs, exp, places=3)


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
    exp = a ** (ms - mi)
    obs = enr.alphastar_i(a, ms, mi)
    assert_almost_equal(obs, exp)


#
# Intergration tests which test multicomponent() and ltot_per_feed()
#


def check_sample_feed(solver):
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.06
    feed = Material(
        {
            922320000: 1.1 * (10.0**-9),
            922340000: 0.00021,
            922350000: 0.0092,
            922360000: 0.0042,
            922380000: 0.9863899989,
        }
    )
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, solver=solver, tolerance=1e-11, max_iter=100)

    # print "casc.mat_prod = " + repr(casc.mat_prod)
    # print "casc.mat_prod = " + repr(casc.mat_tail)
    assert_almost_equal(casc.mat_prod.comp[922350000], 0.06, 5)
    assert_almost_equal(casc.mat_tail.comp[922350000], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0, 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.11652173913043479, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.88347826086956527, 1.0)

    assert_almost_equal(casc.N / 26.864660071132583, 1.0, 5)
    assert_almost_equal(casc.M / 16.637884564470365, 1.0, 5)

    assert_almost_equal(casc.Mstar / 236.57708506549994, 1.0, 5)

    assert_almost_equal(casc.l_t_per_feed / 357.3888391866117, 1.0, 5)
    assert_almost_equal(casc.swu_per_feed / 0.9322804173594426, 1.0, 5)
    assert_almost_equal(casc.swu_per_prod / 8.000914029577306, 1.0, 5)


def test_sample_feed():
    for solver in SOLVERS:
        yield check_sample_feed, solver


def check_NU(solver):
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.05
    feed = Material(
        {
            922340000: 0.000055,
            922350000: 0.00720,
            922380000: 0.992745,
        }
    )
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, solver=solver, tolerance=1e-11)

    # print "casc.mat_prod = " + repr(casc.mat_prod)
    assert_almost_equal(casc.mat_prod.comp[922350000], 0.05, 5)
    assert_almost_equal(casc.mat_tail.comp[922350000], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0, 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.0989473684211, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.901052631579, 1.0)

    assert_almost_equal(casc.N / 27.183583424704818, 1.0, 4)
    assert_almost_equal(casc.M / 13.387464890476533, 1.0, 4)

    assert_almost_equal(casc.Mstar / 236.5621860655, 1.0, 5)

    assert_almost_equal(casc.l_t_per_feed / 288.62731727645644, 1.0, 5)
    assert_almost_equal(casc.swu_per_feed / 0.761263453429, 1.0, 5)
    assert_almost_equal(casc.swu_per_prod / 7.69362000806, 1.0, 5)


def test_NU():
    for solver in SOLVERS:
        yield check_NU, solver


def check_vision(solver):
    orig_casc = enr.default_uranium_cascade()
    orig_casc.x_prod_j = 0.055
    feed = Material(
        {
            922340000: 0.000183963025893197,
            922350000: 0.00818576605617839,
            922360000: 0.00610641667100979,
            922380000: 0.985523854246919,
        }
    )
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, solver=solver, tolerance=1e-11)

    assert_almost_equal(casc.mat_prod.comp[922350000], 0.055, 5)
    assert_almost_equal(casc.mat_tail.comp[922350000], 0.0025, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0, 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.10830030583196934, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.89169969416803063, 1.0)

    assert_almost_equal(casc.N / 27.38162850698868, 1.0, 2)
    assert_almost_equal(casc.M / 15.09646512546496, 1.0, 2)

    assert_almost_equal(casc.Mstar / 236.58177606549995, 1.0, 4)

    assert_almost_equal(casc.l_t_per_feed / 326.8956175003255, 1.0, 4)
    assert_almost_equal(casc.swu_per_feed / 0.85102089049, 1.0, 4)
    assert_almost_equal(casc.swu_per_prod / 7.85797310499, 1.0, 4)


def test_vision():
    for solver in SOLVERS:
        yield check_vision, solver


def check_tungsten(solver):
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

    feed = Material(
        {
            741800000: 0.0014,
            741820000: 0.26416,
            741830000: 0.14409,
            741840000: 0.30618,
            741860000: 0.28417,
        }
    )
    orig_casc.mat_feed = feed
    casc = enr.multicomponent(orig_casc, solver=solver, tolerance=1e-7)

    assert_almost_equal(casc.mat_prod.comp[741800000], 0.5109, 5)
    assert_almost_equal(casc.mat_tail.comp[741800000], 0.00014, 5)

    assert_almost_equal(casc.mat_feed.mass / 1.0, 1.0)
    assert_almost_equal(casc.mat_prod.mass / 0.0024669120526274574, 1.0)
    assert_almost_equal(casc.mat_tail.mass / 0.99753308794737272, 1.0)

    assert_almost_equal(casc.N / 43.557515688533513, 1.0, 2)
    assert_almost_equal(casc.M / 11.49556481009056, 1.0, 2)

    assert_almost_equal(casc.Mstar / 181.16425540249995, 1.0, 4)

    assert_almost_equal(casc.l_t_per_feed / 96.81774564292206, 1.0, 3)
    assert_almost_equal(casc.swu_per_feed / 2.22221945305, 1.0, 3)
    assert_almost_equal(casc.swu_per_prod / 900.810164953, 1.0, 3)


def test_tungsten():
    for solver in SOLVERS:
        yield check_tungsten, solver


if __name__ == "__main__":
    nose.runmodule()
