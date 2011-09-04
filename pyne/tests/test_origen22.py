import os
from StringIO import StringIO

import numpy as np
from nose.tools import assert_equal, assert_true

from pyne import origen22
from pyne.material import Material


def test_sec_to_time_unit():
    assert_equal(origen22.sec_to_time_unit(1.0),  (1.0, 1)) 
    assert_equal(origen22.sec_to_time_unit(10.0),  (10.0, 1)) 

    assert_equal(origen22.sec_to_time_unit(60.0),  (1.0, 2)) 
    assert_equal(origen22.sec_to_time_unit(120.0),  (2.0, 2)) 

    assert_equal(origen22.sec_to_time_unit(np.inf),  (0.0, 6)) 
    assert_equal(origen22.sec_to_time_unit(315569260.0),  (10.0, 5)) 
    assert_equal(origen22.sec_to_time_unit(31556926.0 * 1E7),  (10.0, 8)) 
    assert_equal(origen22.sec_to_time_unit(31556926.0 * 1E10),  (10.0, 9)) 



def test_write_tape4():
    mat = Material({"U235": 0.95, 80160: 0.05})
    tape4 = StringIO()
    origen22.write_tape4(mat, tape4)

    tape4.seek(0)
    observed = tape4.read()
    expected = ("1 80160 5.0000000000E-02   0 0   0 0   0 0\n"
                "2 922350 9.5000000000E-01   0 0   0 0   0 0\n"
                "0 0 0 0\n")

    assert_equal(observed, expected)



def test_out_table_string1():
    obs = origen22._out_table_string(None, None)
    exp = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
    assert_equal(obs, exp)


def test_out_table_string2():
    obs = origen22._out_table_string((False, False, True), None)
    exp = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
    assert_equal(obs, exp)


def test_out_table_string3():
    obs = origen22._out_table_string((False, False, True), range(1, 25))
    exp = "7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7"
    assert_equal(obs, exp)


def test_out_table_string4():
    obs = origen22._out_table_string((False, False, True), [10, 5])
    exp = "8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
    assert_equal(obs, exp)


def test_out_table_string5():
    obs = origen22._out_table_string((True, False, True), [10, 5])
    exp = "8 8 8 8 3 8 8 8 8 3 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
    assert_equal(obs, exp)


def test_write_tape5_irradiation():
    tape5 = StringIO()
    origen22.write_tape5_irradiation("IRP", 100, 0.550, xsfpy_nlb=[204, 205, 206], 
                                     outfile=tape5,
                                     out_table_nes=(False, False, True), 
                                     out_table_laf=(True,  False,  True),  
                                     out_table_num=[5, 10])

    tape5.seek(0)
    observed = tape5.read()

    expected = ("  -1\n"
                "  -1\n"  
                "  -1\n"
                "  CUT     5 1.000E-10 -1\n"
                "  RDA     Make sure thet the library identifier numbers match those in the TAPE9.INP file\n"
                "  LIB     0 1 2 3 204 205 206 9 3 0 4 0\n"
                "  OPTL    8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  OPTA    8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  OPTF    8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  INP     1 -1  0  -1  4  4\n"
                "  RDA     All irradiation (IRF and IRP) cards must be between burnup (BUP) cards.\n"
                "  BUP\n"
                "  IRP     1.0000000000E+02  5.5000000000E-01   1   2   4  2\n"
                "  BUP\n"
                "  OUT     2  1 1 0\n"
                "  END\n")

    assert_equal(observed, expected)



def test_parse_tape6():
    r = origen22.parse_tape6('tape6.test')
    assert(0 < len(r))

    assert_equal(r['time_sec'], 8.64E+06)
    assert_equal(r['flux'], 1.71E+17)
    assert_equal(r['specific_power_MW'], 5.50E-01)
    assert_equal(r['burnup_MWD'], 5.50E+01)
    assert_equal(r['k_inf'], 0.08498)
    assert_equal(r['neutron_production_rate'], 2.97E-04)
    assert_equal(r['neutron_destruction_rate'], 3.50E-03)
    assert_equal(r['total_burnup'], 5.50E+01)
    assert_equal(r['average_flux'], 1.71E+17)
    assert_equal(r['average_specific_power'], 5.50E-01)

    tab_keys = set(['table_{0}'.format(n) for n in range(1, 11) + range(13, 25)])
    assert(tab_keys <=  set(r.keys()))

    for tk in tab_keys:
        assert('nuclide' in  r[tk].keys())
        assert('title' in r[tk]['nuclide'].keys())
        assert('units' in r[tk]['nuclide'].keys())
        assert('data'  in r[tk]['nuclide'].keys())


sample_tape9 = """\
   1      SAMPLE DECAY LIB: ACTIVATION PRODUCTS
   1   10010  6     0.0       0.0       0.0       0.0       0.0       0.0
   1                0.0       0.0       0.0       9.998E+01 1.000E+00 1.000E+00
   1   10020  6     0.0       0.0       0.0       0.0       0.0       0.0
   1                0.0       0.0       0.0       1.500E-02 1.000E+00 1.000E+00
   1   10030  1     3.897E+08 0.0       0.0       0.0       0.0       0.0
   1                0.0       0.0       5.680E-03 0.0       2.000E-07 3.000E-03
   1   10040  1     1.000E-03 0.0       0.0       0.0       0.0       0.0
   1                0.0       0.0       0.0       0.0       3.000E-08 1.000E+00
   1  200410  7     8.100E+01 0.0       1.000E+00 0.0       0.0       0.0
   1                0.0       0.0       2.700E-03 0.0       1.000E-10 3.000E-06
   1  781900  9     6.000E+02 0.0       0.0       0.0       1.000E+00 0.0
   1                0.0       0.0       3.250E+00 1.300E-02 2.000E-14 3.000E-08
  -1
   2      SAMPLE DECAY LIB: ACTINIDES
   2  932410  2     1.600E+01 0.0       0.0       0.0       0.0       0.0
   2                0.0       0.0       4.710E-01 0.0       3.000E-08 1.000E+00
   2  942360  1     8.997E+07 0.0       0.0       0.0       1.000E+00 0.0
   2                8.000E-10 0.0       5.871E+00 0.0       6.000E-13 3.000E-05
   2  942370  4     4.560E+01 0.0       1.000E+00 0.0       3.300E-05 0.0
   2                0.0       0.0       6.220E-02 0.0       2.000E-14 3.000E-08
  -1
   3      SAMPLE DECAY LIB: FISSION PRODUCTS
   3  611460  5     5.500E+00 0.0       6.300E-01 0.0       0.0       0.0
   3                0.0       0.0       8.508E-01 0.0       1.000E-10 3.000E-06
   3  621460  8     7.000E+01 0.0       0.0       0.0       1.000E+00 0.0
   3                0.0       0.0       2.540E+00 0.0       2.000E-14 3.000E-08
   3  691720  3     6.360E+01 0.0       0.0       0.0       0.0       0.0
   3                0.0       0.0       1.880E+00 0.0       1.000E-10 3.000E-06
  -1
 381        SAMPLE ACTIVATION PRODUCT XS LIB
 381   10010 1.550E-04 0.0       0.0       0.0       0.0       0.0         -1.0
 381   10020 3.100E-07 6.757E-04 0.0       0.0       0.0       0.0         -1.0
 381   10030 3.510E-09 0.0       0.0       0.0       0.0       0.0         -1.0
 381   20030 3.107E+00 0.0       0.0       3.255E+00 0.0       0.0         -1.0
 381   30060 1.638E-05 0.0       4.466E-01 7.184E-04 0.0       0.0         -1.0
 381   30070 9.100E-06 0.0       4.481E-03 0.0       0.0       0.0         -1.0
 381   40090 5.200E-06 3.738E-02 2.168E-02 0.0       0.0       0.0         -1.0
 381   40100 5.850E-07 0.0       0.0       0.0       0.0       0.0         -1.0
 381   50100 2.925E-04 0.0       2.114E+00 1.566E-03 0.0       0.0         -1.0
 381   50110 3.071E-05 6.503E-07 8.573E-06 8.093E-08 0.0       0.0         -1.0
 381   60120 3.018E-06 0.0       2.189E-04 4.895E-08 0.0       0.0         -1.0
 381   60130 1.690E-06 0.0       4.516E-04 0.0       0.0       0.0         -1.0
 381   60140 5.850E-10 0.0       0.0       0.0       0.0       0.0         -1.0
 381   70140 4.926E-05 3.198E-07 1.500E-02 1.533E-02 0.0       0.0         -1.0
 381   70150 1.315E-05 8.713E-06 2.426E-05 5.091E-06 0.0       0.0         -1.0
 381   80160 1.041E-07 0.0       1.364E-03 5.538E-06 0.0       0.0         -1.0
 381   80170 7.341E-05 9.378E-06 3.639E-02 3.310E-06 0.0       0.0         -1.0
 381   80180 1.053E-06 0.0       2.605E-04 0.0       0.0       0.0         -1.0
  -1
 382       SAMPLE ACTINIDE XS LIB
 382  862200 1.170E-04 0.0       0.0       0.0       0.0       0.0         -1.0
 382  862220 4.212E-04 0.0       0.0       0.0       0.0       0.0         -1.0
 382  922350 4.202E-01 1.326E-03 1.556E-06 1.637E+00 0.0       0.0         -1.0
 382  922360 4.292E-01 1.475E-03 2.529E-05 1.224E-01 0.0       0.0         -1.0
 382  922370 3.661E-01 3.986E-03 5.062E-05 6.297E-01 0.0       0.0         -1.0
 382  922380 2.125E-01 2.731E-03 2.132E-05 4.976E-02 0.0       0.0         -1.0
  -1
 383       SAMPLE FISSION PRODUCT YIELD
 383   10030 3.510E-09 0.0       0.0       0.0       0.0       0.0          1.0
 383     2.00E-02 2.00E-02 2.00E-02 2.30E-02 1.75E-02 1.75E-02 1.75E-02 1.75E-02
 383   30060 1.638E-05 0.0       4.466E-01 7.184E-04 0.0       0.0          1.0
 383     5.00E-05 5.00E-05 5.00E-05 5.00E-05 5.00E-05 5.00E-05 5.00E-05 5.00E-05
 383   30070 9.100E-06 0.0       4.481E-03 0.0       0.0       0.0          1.0
 383     1.00E-06 1.00E-06 1.00E-06 1.00E-06 1.00E-06 1.00E-06 1.00E-06 1.00E-06
 383   40090 5.200E-06 3.738E-02 2.168E-02 0.0       0.0       0.0          1.0
 383     1.50E-06 1.50E-06 1.50E-06 1.50E-06 1.50E-06 1.50E-06 1.50E-06 1.50E-06
 383   40100 5.850E-07 0.0       0.0       0.0       0.0       0.0          1.0
 383     9.00E-06 9.00E-06 9.00E-06 9.00E-06 9.00E-06 9.00E-06 9.00E-06 9.00E-06
 383   60140 5.850E-10 0.0       0.0       0.0       0.0       0.0          1.0
 383     1.30E-06 1.30E-06 1.30E-06 1.30E-06 1.30E-06 1.30E-06 1.30E-06 1.30E-06
 383  290660 7.897E-02 0.0       0.0       0.0       0.0       0.0          1.0
 383     5.97E-12 3.17E-08 9.26E-09 4.12E-10 3.55E-10 1.71E-09 1.68E-09 1.68E-09
 383  300660 1.040E-03 0.0       1.193E-08 0.0       0.0       0.0          1.0
 383     0.0      0.0      2.55E-11 0.0      0.0      0.0      0.0      0.0
 383  300670 2.600E-02 0.0       3.512E-09 0.0       0.0       0.0          1.0
 383     0.0      2.50E-09 3.84E-10 4.60E-12 2.14E-11 0.0      0.0      0.0
 383  300680 4.001E-03 0.0       1.192E-08 0.0       2.883E-04 0.0         -1.0
 383  310690 2.028E-02 0.0       0.0       0.0       0.0       0.0         -1.0
 383  300700 4.855E-05 0.0       0.0       0.0       5.091E-06 0.0         -1.0
  -1
"""

def test_parse_tape9():
    tape9_file = StringIO(sample_tape9)
    tape9 = origen22.parse_tape9(tape9_file)

    assert_equal(set(tape9.keys()), set([1, 2, 3, 381, 382, 383]))

    # Activation product decay
    deck1 = tape9[1]
    assert_equal(deck1['_type'], 'decay')
    assert_equal(deck1['title'], 'SAMPLE DECAY LIB: ACTIVATION PRODUCTS')
    assert_equal(deck1['half_life'][10020], np.inf)
    assert_equal(deck1['half_life'][10040], 1.000E-03)
    assert_equal(deck1['half_life'][200410], 8.100E+01 * 31556926.0 * 1E3)
    assert_equal(deck1['half_life'][781900], 6.000E+02 * 31556926.0 * 1E9)
    assert_equal(deck1['frac_beta_minus_x'][10010], 0.0)
    assert_equal(deck1['frac_beta_plus_or_electron_capture'][200410], 1.0)
    assert_equal(deck1['frac_beta_plus_or_electron_capture_x'][10010], 0.0)
    assert_equal(deck1['frac_alpha'][781900], 1.0)
    assert_equal(deck1['frac_isomeric_transition'][10020], 0.0)
    assert_equal(deck1['frac_spont_fiss'][10020], 0.0)
    assert_equal(deck1['frac_beta_n'][10020], 0.0)
    assert_equal(deck1['recoverable_energy'][10030], 5.680E-03)
    assert_equal(deck1['frac_natural_abund'][10020], 1.500E-02 * 0.01)
    assert_equal(deck1['inhilation_concentration'][781900], 2.000E-14)
    assert_equal(deck1['ingestion_concentration'][781900], 3.000E-08)


    # Actinide Decay
    deck2 = tape9[2]
    assert_equal(deck2['_type'], 'decay')
    assert_equal(deck2['title'], 'SAMPLE DECAY LIB: ACTINIDES')
    assert_equal(deck2['half_life'][932410], 1.600E+01 * 60.0)
    assert_equal(deck2['half_life'][942370], 4.560E+01 * 86400.0)


    # Fission Product Decay
    deck3 = tape9[3]
    assert_equal(deck3['_type'], 'decay')
    assert_equal(deck3['title'], 'SAMPLE DECAY LIB: FISSION PRODUCTS')
    assert_equal(deck3['half_life'][611460], 5.500E+00 * 31556926.0)
    assert_equal(deck3['half_life'][621460], 7.000E+01 * 31556926.0 * 1E6)
    assert_equal(deck3['half_life'][691720], 6.360E+01 * 3600.0)


    # Activation product cross sections
    deck381 = tape9[381]
    assert_equal(deck381['_type'], 'xsfpy')
    assert_equal(deck381['_subtype'], 'activation_products')
    assert_equal(deck381['title'], 'SAMPLE ACTIVATION PRODUCT XS LIB')
    assert_true(all(['_fiss_yield' not in key for key in deck381.keys()]))

    assert_true('sigma_alpha' in deck381)
    assert_true('sigma_3n' not in deck381)

    assert_true('sigma_p' in deck381)
    assert_true('sigma_f' not in deck381)

    assert_equal(deck381['sigma_gamma'][80170], 7.341E-05)
    assert_equal(deck381['sigma_2n'][80170], 9.378E-06)
    assert_equal(deck381['sigma_alpha'][80170], 3.639E-02)
    assert_equal(deck381['sigma_p'][80170], 3.310E-06)
    assert_equal(deck381['sigma_gamma_x'][80170], 0.0)
    assert_equal(deck381['sigma_2n_x'][80170], 0.0)
    assert_equal(deck381['fiss_yields_present'][80170], False)


    # Actinide cross sections
    deck382 = tape9[382]
    assert_equal(deck382['_type'], 'xsfpy')
    assert_equal(deck382['_subtype'], 'actinides')
    assert_equal(deck382['title'], 'SAMPLE ACTINIDE XS LIB')
    assert_true(all(['_fiss_yield' not in key for key in deck382.keys()]))

    assert_true('sigma_alpha' not in deck382)
    assert_true('sigma_3n' in deck382)

    assert_true('sigma_p' not in deck382)
    assert_true('sigma_f' in deck382)

    assert_equal(deck382['sigma_gamma'][922380], 2.125E-01)
    assert_equal(deck382['sigma_2n'][922380], 2.731E-03)
    assert_equal(deck382['sigma_3n'][922380], 2.132E-05)
    assert_equal(deck382['sigma_f'][922380], 4.976E-02)
    assert_equal(deck382['sigma_gamma_x'][922380], 0.0)
    assert_equal(deck382['sigma_2n_x'][922380], 0.0)
    assert_equal(deck382['fiss_yields_present'][922380], False)


    # Fission product cross sections
    deck383 = tape9[383]
    assert_equal(deck383['_type'], 'xsfpy')
    assert_equal(deck383['_subtype'], 'fission_products')
    assert_equal(deck383['title'], 'SAMPLE FISSION PRODUCT YIELD')
    assert_true(any(['_fiss_yield' in key for key in deck383.keys()]))

    assert_true('sigma_alpha' in deck383)
    assert_true('sigma_3n' not in deck383)

    assert_true('sigma_p' in deck383)
    assert_true('sigma_f' not in deck383)

    assert_equal(deck383['sigma_gamma'][300670], 2.600E-02)
    assert_equal(deck383['sigma_2n'][300670], 0.0)
    assert_equal(deck383['sigma_alpha'][300670], 3.512E-09)
    assert_equal(deck383['sigma_p'][300670], 0.0)
    assert_equal(deck383['sigma_gamma_x'][300670], 0.0)
    assert_equal(deck383['sigma_2n_x'][300670], 0.0)
    assert_equal(deck383['fiss_yields_present'][300670], True)
    assert_equal(deck383['TH232_fiss_yield'][300670], 0.0)
    assert_equal(deck383['U233_fiss_yield'][300670], 2.50E-09)
    assert_equal(deck383['U235_fiss_yield'][300670], 3.84E-10)
    assert_equal(deck383['U238_fiss_yield'][300670], 4.60E-12)
    assert_equal(deck383['PU239_fiss_yield'][300670], 2.14E-11)
    assert_equal(deck383['PU241_fiss_yield'][300670], 0.0)
    assert_equal(deck383['CM245_fiss_yield'][300670], 0.0)
    assert_equal(deck383['CF249_fiss_yield'][300670], 0.0)



def test_merge_tape9():
    tape9_file = StringIO(sample_tape9)
    tape9_file = origen22.parse_tape9(tape9_file)

    tape9_dict = {1: {'_type': 'decay', 'half_life': {10010: 42.0}},
                  2: {'_type': 'decay', '_bad_key': None},
                  3: {'_type': 'decay', 'title': "Sweet Decay"},
                  382: {'_type': 'xsfpy', '_subtype': 'actinides', 'sigma_f': {922350: 16.0}},
                 }

    # merge so that dict takes precedence
    tape9 = origen22.merge_tape9([tape9_dict, tape9_file])

    # run tests
    assert_equal(tape9[1]['half_life'][10010], 42.0)
    assert_true('_bad_key' in tape9[2])
    assert_equal(tape9[3]['title'], "Sweet Decay")
    assert_equal(tape9[382]['sigma_f'][922350], 16.0)

    assert_true('_cards' not in tape9[1])
    assert_true('_cards' not in tape9[2])
    assert_true('_cards' not in tape9[3])
    assert_true('_cards' in tape9[381])
    assert_true('_cards' not in tape9[382])
    assert_true('_cards' in tape9[383])



def test_write_tape9():
    tape9_file = StringIO()

    tape9_dict = {1: {'_type': 'decay', 'half_life': {10010: 42.0}, 'title': 'decay1'},
                  2: {'_type': 'decay', '_bad_key': None, 'title': 'decay2'},
                  3: {'_type': 'decay', 'title': "Sweet Decay"},
                  381: {'_type': 'xsfpy', '_subtype': 'activation_products', 'sigma_gamma': {10010: 12.0}, 'title': 'xs1'},
                  382: {'_type': 'xsfpy', '_subtype': 'actinides', 'sigma_f': {922350: 16.0}, 'title': 'xs2'},
                  383: {'_type': 'xsfpy', '_subtype': 'fission_products', 'sigma_gamma': {10010: 20.0}, 'title': 'xsfpy3'},
                 }

    # Test that basic functionality works
    origen22.write_tape9(tape9_dict, tape9_file)
    tape9_file.seek(0)
    t9str = tape9_file.read()

    # Try to round-trip
    full_tape9_file = StringIO(sample_tape9)
    full_tape9 = origen22.parse_tape9(tape9_file)

    backout_tape9 = StringIO()
    origen22.write_tape9(full_tape9, backout_tape9)
    backout_tape9.seek(0)

    backin_tape9 = origen22.parse_tape9(backout_tape9)
    assert_equal(full_tape9, backin_tape9)
