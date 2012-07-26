"""Enrichment Tests"""
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, with_setup

import os
import warnings
import numpy as np

from pyne.material import Material
from pyne import enrichment as enr

#
# Fixtures
#

def setup_enrichment():
    bright_conf.track_nucs = set([922350, 922360, 922380])
    bright_conf.verbosity = 0


def setup_enrichment1():
    bright_conf.track_nucs = set([922320, 922340, 922350, 922360, 922380])
    bright_conf.verbosity = 0


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
                                922380: 0.992745}, 1.0, 'Natural Uranium', 1.0))

"""\
#
# Tests that the Enrichment component constructors work.
#

@with_setup(None, teardown_enrichment)
def test_Enrichment_1():
    e = Enrichment()
    assert_equal(e.name, '')
    assert_equal(e.track_params, set(["M",  "MassFeed", "MassProduct", "MassTails", 
                                      "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))

@with_setup(None, teardown_enrichment)
def test_Enrichment_2():
    e = Enrichment(name="e")
    assert_equal(e.name, 'e')
    assert_equal(e.track_params, set(["M",  "MassFeed", "MassProduct", "MassTails", 
                                      "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))

@with_setup(None, teardown_enrichment)
def test_Enrichment_3():
    casc = bright.enrichment.uranium_enrichment_defaults()
    casc.x_prod_j = 0.1
    e = Enrichment(enrich_params=ep)
    assert_equal(e.name, '')
    assert_equal(e.track_params, set(["M",  "MassFeed", "MassProduct", "MassTails", 
                                      "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
    assert_equal(e.alpha, 1.05)
    assert_equal(e.Mstar, 236.5)
    assert_equal(e.j, 922350)
    assert_equal(e.k, 922380)
    assert_equal(e.N, 30.0)
    assert_equal(e.M, 10.0)
    assert_equal(e.x_prod_j, 0.1)
    assert_equal(e.x_tail_j, 0.0025)

@with_setup(None, teardown_enrichment)
def test_Enrichment_4():
    casc = bright.enrichment.uranium_enrichment_defaults()
    casc.j = 922360
    e = Enrichment(enrich_params=ep, name='e')
    assert_equal(e.name, 'e')
    assert_equal(e.track_params, set(["M",  "MassFeed", "MassProduct", "MassTails", 
                                      "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
    assert_equal(e.alpha, 1.05)
    assert_equal(e.Mstar, 236.5)
    assert_equal(e.j, 922360)
    assert_equal(e.k, 922380)
    assert_equal(e.N, 30.0)
    assert_equal(e.M, 10.0)
    assert_equal(e.x_prod_j, 0.05)
    assert_equal(e.x_tail_j, 0.0025)


#
# Tests that enrichment the fuel cycle component attributes work.
#


@with_setup(None, teardown_enrichment)
def test_track_params():
    e = Enrichment()
    assert_equal(e.track_params, set(["M",  "MassFeed", "MassProduct", "MassTails", 
                                      "Mstar", "N", "SWUperFeed", "SWUperProduct", "TotalPerFeed"]))
    e.track_params = set(["Om nom nom"])
    assert_equal(e.track_params, set(["Om nom nom"]))
                        
@with_setup(None, teardown_enrichment)
def test_alpha():
    e = Enrichment()
    e.alpha = 1.06
    assert_equal(e.alpha, 1.06)

@with_setup(None, teardown_enrichment)
def test_Mstar():
    e = Enrichment()
    e.Mstar = 235.5
    assert_equal(e.Mstar, 235.5)

@with_setup(None, teardown_enrichment)
def test_j():
    e = Enrichment()
    e.j = 922360
    assert_equal(e.j, 922360)

@with_setup(None, teardown_enrichment)
def test_k():
    e = Enrichment()
    e.k = 922390
    assert_equal(e.k, 922390)

@with_setup(None, teardown_enrichment)
def test_N():
    e = Enrichment()
    e.N = 35.0
    assert_equal(e.N, 35.0)

@with_setup(None, teardown_enrichment)
def test_M():
    e = Enrichment()
    e.M = 15.0
    assert_equal(e.M, 15.0)

@with_setup(None, teardown_enrichment)
def test_x_prod_j():
    e = Enrichment()
    e.x_prod_j = 0.1
    assert_equal(e.x_prod_j, 0.1)

@with_setup(None, teardown_enrichment)
def test_x_tail_j():
    e = Enrichment()
    e.x_tail_j = 0.005
    assert_equal(e.x_tail_j, 0.005)

@with_setup(None, teardown_enrichment)
def test_Mstar():
    e = Enrichment()
    e.Mstar = 235.5
    assert_equal(e.Mstar, 235.5)

@with_setup(None, teardown_enrichment)
def test_mat_feed():
    e = Enrichment()
    mat = Material({922380: 0.5})        
    e.mat_feed = mat
    assert_equal(e.mat_feed.mass, 0.5)
    assert_equal(e.mat_feed.comp[922380], 1.0)

@with_setup(None, teardown_enrichment)
def test_mat_prod():
    e = Enrichment()
    mat = Material({922380: 0.5})        
    e.mat_prod = mat
    assert_equal(e.mat_prod.mass, 0.5)
    assert_equal(e.mat_prod.comp[922380], 1.0)

@with_setup(None, teardown_enrichment)
def test_mat_tail():
    e = Enrichment()
    mat = Material({922380: 0.5})        
    e.mat_tail = mat
    assert_equal(e.mat_tail.mass, 0.5)
    assert_equal(e.mat_tail.comp[922380], 1.0)    

@with_setup(None, teardown_enrichment)
def test_N():
    e = Enrichment()
    e.N = 35.0
    assert_equal(e.N, 35.0)

@with_setup(None, teardown_enrichment)
def test_M():
    e = Enrichment()
    e.M = 15.0
    assert_equal(e.M, 15.0)

@with_setup(None, teardown_enrichment)
def test_TotalPerFeed():
    e = Enrichment()
    e.TotalPerFeed = 12.0
    assert_equal(e.TotalPerFeed, 12.0)

@with_setup(None, teardown_enrichment)
def test_SWUperFeed():
    e = Enrichment()
    e.SWUperFeed = 101.0
    assert_equal(e.SWUperFeed, 101.0)

@with_setup(None, teardown_enrichment)
def test_SWUperProduct():
    e = Enrichment()
    e.SWUperProduct = 1010.0
    assert_equal(e.SWUperProduct, 1010.0)


#
# Tests that the fuel cycle component methods work.
#

@with_setup(setup_enrichment, teardown_enrichment)
def test_initialize():
    e = Enrichment()
    assert_equal(e.x_prod_j, 0.05)
    casc = Cascade()
    casc.x_prod_j = 0.1
    e.initialize(ep)
    assert_almost_equal(e.alpha, 0.0)
    assert_almost_equal(e.Mstar, 0.0)
    assert_equal(e.j, 0)
    assert_equal(e.k, 0)
    assert_almost_equal(e.N, 0.0)
    assert_almost_equal(e.M, 0.0)
    assert_equal(e.x_prod_j, 0.1)
    assert_almost_equal(e.x_tail_j, 0.0)


@with_setup(setup_enrichment, teardown_enrichment)
def test_calc_1():
    e = Enrichment()
    e.mat_feed = Material({922350: 0.01, 922380: 0.985, 922360: 0.005})
    e.calc()
    assert_almost_equal(e.mat_prod.comp[922350],  0.05,   6) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 6)

@with_setup(setup_enrichment, teardown_enrichment)
def test_calc_2():
    e = Enrichment()
    e.calc({922350: 0.01, 922380: 0.985, 922360: 0.005})
    assert_almost_equal(e.mat_prod.comp[922350],  0.05,   6) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 6)

@with_setup(setup_enrichment, teardown_enrichment)
def test_calc_3():
    e = Enrichment()
    mat = Material({922350: 0.01, 922380: 0.985, 922360: 0.005})
    e.calc(mat)
    assert_almost_equal(e.mat_prod.comp[922350],  0.05,   6) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 6)

@with_setup(setup_enrichment, teardown_enrichment)
def test_calc_params():
    e = Enrichment()
    mat = Material({922350: 0.01, 922380: 0.985, 922360: 0.005})
    e.calc(mat)
    e.calc_params()

    assert_equal(e.params_prior_calc["MassFeed"],  e.mat_feed.mass)
    assert_equal(e.params_after_calc["MassFeed"], 0.0)

    assert_equal(e.params_prior_calc["MassProduct"],  0.0)
    assert_equal(e.params_after_calc["MassProduct"], e.mat_prod.mass)

    assert_equal(e.params_prior_calc["MassTails"],  0.0)
    assert_equal(e.params_after_calc["MassTails"], e.mat_tail.mass)

    assert_equal(e.params_prior_calc["N"],  e.N)
    assert_equal(e.params_after_calc["N"], e.N)

    assert_equal(e.params_prior_calc["M"],  e.M)
    assert_equal(e.params_after_calc["M"], e.M)

    assert_equal(e.params_prior_calc["Mstar"],  e.Mstar)
    assert_equal(e.params_after_calc["Mstar"], e.Mstar)

    assert_equal(e.params_prior_calc["TotalPerFeed"],  e.TotalPerFeed)
    assert_equal(e.params_after_calc["TotalPerFeed"], e.TotalPerFeed)

    assert_equal(e.params_prior_calc["SWUperFeed"],  e.SWUperFeed)
    assert_equal(e.params_after_calc["SWUperFeed"], 0.0)

    assert_equal(e.params_prior_calc["SWUperProduct"],  0.0)
    assert_equal(e.params_after_calc["SWUperProduct"], e.SWUperProduct)



@with_setup(setup_enrichment1, teardown_enrichment)
def test_sample_feed():
    casc = bright.enrichment.uranium_enrichment_defaults()
    casc.x_prod_j = 0.06
    e = Enrichment(enrich_params=ep, name='e')
    mat = Material({
            922320: 1.1 * (10.0**-9),
            922340: 0.00021,
            922350: 0.0092,
            922360: 0.0042,
            922380: 0.9863899989,
            })
    e.calc(mat)

    assert_almost_equal(e.mat_prod.comp[922350],  0.06,   5) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(e.mat_feed.mass / 1.0,                 1.0)
    assert_almost_equal(e.mat_prod.mass / 0.11652173913043479, 1.0)
    assert_almost_equal(e.mat_tail.mass / 0.88347826086956527, 1.0)

    assert_almost_equal(e.N / 26.8646352802, 1.0, 4)
    assert_almost_equal(e.M / 16.6379009423, 1.0, 4)

    assert_almost_equal(e.Mstar / 236.577085, 1.0, 5)

    assert_almost_equal(e.TotalPerFeed  / 357.388791749,  1.0, 5)
    assert_almost_equal(e.SWUperFeed    / 0.932280175218, 1.0, 5)
    assert_almost_equal(e.SWUperProduct / 8.0009119515,  1.0, 5)


@with_setup(setup_enrichment1, teardown_enrichment)
def test_NU():
    casc = bright.enrichment.uranium_enrichment_defaults()
    casc.x_prod_j = 0.05
    e = Enrichment(enrich_params=ep, name='e')
    mat = Material({
            922340: 0.000055,
            922350: 0.00720,
            922380: 0.992745,
            })
    e.calc(mat)

    assert_almost_equal(e.mat_prod.comp[922350],  0.05,   5) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(e.mat_feed.mass   / 1.0,             1.0)
    assert_almost_equal(e.mat_prod.mass  / 0.0989473684211, 1.0)
    assert_almost_equal(e.mat_tail.mass / 0.901052631579,  1.0)

    assert_almost_equal(e.N / 27.1835088212, 1.0, 4)
    assert_almost_equal(e.M / 13.3875092512, 1.0, 4)

    assert_almost_equal(e.Mstar / 236.562179, 1.0, 5)

    assert_almost_equal(e.TotalPerFeed  / 288.627270162,  1.0, 5)
    assert_almost_equal(e.SWUperFeed    / 0.761263453429, 1.0, 5)
    assert_almost_equal(e.SWUperProduct / 7.69362000806,  1.0, 5)


@with_setup(setup_enrichment1, teardown_enrichment)
def test_VISION():
    casc = bright.enrichment.uranium_enrichment_defaults()
    casc.x_prod_j = 0.055
    e = Enrichment(enrich_params=ep, name='e')
    mat = Material({
            922340: 0.000183963025893197,
            922350: 0.00818576605617839,
            922360: 0.00610641667100979,
            922380: 0.985523854246919,
            })
    e.calc(mat)

    assert_almost_equal(e.mat_prod.comp[922350],  0.055,   5) 
    assert_almost_equal(e.mat_tail.comp[922350], 0.0025, 5)

    assert_almost_equal(e.mat_feed.mass / 1.0,                 1.0)
    assert_almost_equal(e.mat_prod.mass / 0.10830030583196934, 1.0)
    assert_almost_equal(e.mat_tail.mass / 0.89169969416803063, 1.0)

    assert_almost_equal(e.N / 27.38162850698868, 1.0, 2)
    assert_almost_equal(e.M / 15.09646512546496, 1.0, 2)

    assert_almost_equal(e.Mstar / 236.581784, 1.0, 4)

    print e.TotalPerFeed
    print e.SWUperFeed
    print e.SWUperProduct
    assert_almost_equal(e.TotalPerFeed  / 326.895568684,  1.0, 4)
    assert_almost_equal(e.SWUperFeed    / 0.85102089049,  1.0, 4)
    assert_almost_equal(e.SWUperProduct / 7.85797310499,  1.0, 4)


@with_setup(setup_enrichment1, teardown_enrichment)
def test_Tungsten():
    # This test comes from 'Multicomponent Isotope Separation in Matched
    # Abundance Ratio Cascades Composed of Stages with Large Separation Factors' 
    # by E. von Halle, 1987.

    casc = Cascade()
    casc.alpha = 1.16306
    casc.Mstar = 181.3
    casc.j = 741800
    casc.k = 741860
    casc.N = 30.0
    casc.M = 10.0
    casc.x_prod_j = 0.5109
    casc.x_tail_j = 0.00014

    e = Enrichment(enrich_params=ep, name='e')
    mat = Material({
            741800: 0.0014, 
            741820: 0.26416, 
            741830: 0.14409, 
            741840: 0.30618, 
            741860: 0.28417,
            })
    e.calc(mat)

    assert_almost_equal(e.mat_prod.comp[741800],  0.5109,  5) 
    assert_almost_equal(e.mat_tail.comp[741800], 0.00014, 5)

    assert_almost_equal(e.mat_feed.mass / 1.0,                   1.0)
    assert_almost_equal(e.mat_prod.mass / 0.0024669120526274574, 1.0)
    assert_almost_equal(e.mat_tail.mass / 0.99753308794737272,   1.0)

    assert_almost_equal(e.N / 43.557515688533513, 1.0, 2)
    assert_almost_equal(e.M / 11.49556481009056,  1.0, 2)

    assert_almost_equal(e.Mstar / 181.164592, 1.0, 4)

    assert_almost_equal(e.TotalPerFeed  / 96.8179316719, 1.0, 3)
    assert_almost_equal(e.SWUperFeed    / 2.22221945305, 1.0, 3)
    assert_almost_equal(e.SWUperProduct / 900.810164953, 1.0, 3)

if __name__ == "__main__":
    nose.main()
"""
