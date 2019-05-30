"""openmc tests"""
from __future__ import unicode_literals, division
from io import StringIO
import warnings

import os
import nose
import numpy as np
from math import isclose
from nose.tools import assert_equal, assert_true
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne import nucname
from pyne import openmc
from pyne.mesh import HAVE_PYMOAB, Mesh

sample_xs = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
</cross_sections>
""")

sample_xs_with_dir = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <directory>/</directory>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
""")

sample_xs_with_mcnp_id = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="Co-58m.70c" awr="57.4381" location="28897" metastable="1" name="27458.70c" path="endf70b" temperature="2.5301e-08" zaid="27458"/>
  <ace_table alias="Co-58.70c" awr="57.4381" location="28699" name="27058.70c" path="endf70b" temperature="2.5301e-08" zaid="27058"/>
</cross_sections>
""")

cwd = os.getcwd()
def test_ace_table_init():
    atab = openmc.AceTable(zaid='92235', path='U235.ace',
                           cross_sections_path='/tmp/cross_sections.xml')
    assert_equal('92235', atab.zaid)
    assert_equal('U235.ace', atab.path)
    assert_equal('/tmp/U235.ace', atab.abspath)
    assert_equal(nucname.id('U235'), atab.nucid)


def test_ace_table_xml():
    atab = openmc.AceTable(zaid='92235', path='U235.ace',
                           cross_sections_path='/tmp/cross_sections.xml')
    exp = '<ace_table path="U235.ace" zaid="92235"/>'
    obs = atab.xml()
    assert_equal(exp, obs)


def test_cross_sections_read():
    sample_xs.seek(0)
    xs = openmc.CrossSections(sample_xs)
    assert_equal('ascii', xs.filetype)
    assert_true(xs.path is None)

    exp = [openmc.AceTable(alias='H-1.71c', awr='0.999167', location='1',
                           name='1001.71c', path='293.6K/H_001_293.6K.ace',
                           temperature='2.53e-08', zaid='1001'),
           openmc.AceTable(alias='Am-242m.73c', awr='239.9801', location='1',
                           metastable='1', name='95242.73c',
                           path='900K/Am_242_900K.ace', temperature='7.756e-08',
                           zaid='95242'),
           openmc.AceTable(awr='89.1324', location='1', name='ZrZrH.71t',
                           path='tsl/zrzrh.acer', temperature='2.551e-08',
                           zaid='0')]
    assert_equal(exp, xs.ace_tables)


def test_cross_sections_abspath_with_dir():
    xs = openmc.CrossSections(sample_xs_with_dir)
    assert_equal('ascii', xs.filetype)
    assert_equal(xs.path, "/")

    exp_abspaths = ["/293.6K/H_001_293.6K.ace",
                    "/900K/Am_242_900K.ace", "/tsl/zrzrh.acer"]
    obs_abspaths = [table.abspath for table in xs.ace_tables]
    assert_equal(exp_abspaths, obs_abspaths)


def test_cross_sections_mcnp_id():
    xstables = openmc.CrossSections(sample_xs_with_mcnp_id).ace_tables
    mcnp_obs = [table.nucid for table in xstables if table.alias == "Co-58m.70c"][0]
    assert_equal(mcnp_obs, 270580001)
    nucid_obs = [table.nucid for table in xstables if table.alias == "Co-58.70c"][0]
    assert_equal(nucid_obs, 270580000)


def test_cross_sections_roundtrip():
    sample_xs.seek(0)
    xs = openmc.CrossSections(sample_xs)
    sample_xs.seek(0)
    exp = sample_xs.read()
    obs = xs.xml()
    assert_equal(exp, obs)

def test_get_mesh_name():
    mesh_str = """{'mesh 14': /tallies/meshes/mesh 14 (Group) ''\n
                  children := ['dimension' (Array), 'lower_left' (Array), 
                  'type' (Array), 'upper_right' (Array), 'width' (Array)]}"""
    exp_mesh_name = "mesh 14"
    assert_equal(openmc._get_mesh_name(mesh_str), exp_mesh_name)

def test_get_tally_results_from_openmc_sp():
    # energy bin: [0.0, 1.0, 20.0], 2bins
    # 6 volume elemenes
    filename = os.path.join(cwd, "files_test_openmc", "statepoint.ebin2.ves6.h5")
    # dimension of the tally_results: (12, 1, 2)
    tally_results = openmc.get_tally_results_from_openmc_sp(filename,
            tally_num=1)
    exp_tally_results = np.array(
           [[[0.977887, 0.956263 ]],
            [[0.115698, 0.0133861]], 
            [[0.181645, 0.0329949]],
            [[0.938374, 0.880546 ]],
            [[0.219677, 0.0482578]],
            [[0.18134 , 0.0328842]],
            [[4.20556 , 17.6867  ]],
            [[0       , 0        ]],
            [[0       , 0        ]],
            [[3.40724 , 11.6093  ]],
            [[0       , 0        ]],
            [[0       , 0        ]]], dtype=float)
    # compare
    assert_array_equal(tally_results.shape, exp_tally_results.shape)
    for i in range(tally_results.shape[0]):
        for j in range(tally_results.shape[1]):
            for k in  range(tally_results.shape[2]):
                # data in h5 file is 6 digits
                assert_true(isclose(tally_results[i][j][k],
                    exp_tally_results[i][j][k], rel_tol=1e-5))

def test_calc_structured_coords():
    lower_left = np.array([0.0, 0.0, 0.0])
    upper_right = np.array([1.0, 2.0, 3.0])
    dimension = np.array([2, 4, 5])
    exp_structured_coords = np.array([[0.0, 0.5, 1.0],
                                      [0.0, 0.5, 1.0, 1.5, 2.0],
                                      [0.0, 0.6, 1.2, 1.8, 2.4, 3.0]])
    structured_coords = openmc.calc_structured_coords(lower_left,
            upper_right, dimension)
    assert_equal(len(structured_coords), len(exp_structured_coords))
    for i in range(len(exp_structured_coords)):
        assert_array_almost_equal(structured_coords[i],
                exp_structured_coords[i])

def test_meshtally_from_openmc_statepoint():
    if not HAVE_PYMOAB:
        raise SkipTest
    # mesh read from openmc state point file
    # Parameters of the tally and mesh
    # mesh = openmc.Mesh(mesh_id=14, name="n_flux")
    # mesh.dimension= [3, 2, 1]
    # mesh.lower_left = (-40.0, -12.5, -2.5)
    # mesh.upper_right = (40.0, 12.5, 2.5)
    # energy_bins = np.array([0.0, 1.0, 20.0]) * 1e6
    filename = os.path.join(cwd, "files_test_openmc", "statepoint.ebin2.ves6.h5")
    tally_num = 1
    
    mesh = openmc.meshtally_from_openmc_statepoint(filename, tally_num)
    # check mesh attributes
    assert_equal(len(mesh), 6)
    assert(mesh.structured)
    # structured_coords
    assert_array_almost_equal(mesh.structured_coords[0],
            [(-40.0 + x * 80 / 3) for x in range(0, 4)])
    assert_array_almost_equal(mesh.structured_coords[1],
            [(-12.5 + x * 25.0 / 2) for x in range (0, 3)])
    assert_array_almost_equal(mesh.structured_coords[2],
            [(-2.5 + x * 5 / 1) for x in range(0, 2)])
    ve_vol = (80.0/3) * (25.0/2) * (5/1)
    # flux
    exp_n_flux = np.divide(np.array([[0.977887, 4.20556],
                           [0.115698, 0.0],
                           [0.181645, 0.0],
                           [0.938374, 3.40724],
                           [0.219677, 0.0],
                           [0.18134,  0.0]]), ve_vol)
    assert_equal(len(mesh.n_flux[:]), len(exp_n_flux))
    for i in range(len(exp_n_flux)):
        assert_array_almost_equal(mesh.n_flux[i], exp_n_flux[i], decimal=5)
    # error
    exp_n_flux_err = np.divide(np.array([[0.956263, 17.6867],
                               [0.0133861, 0.0],
                               [0.0329949, 0.0],
                               [0.880546, 11.6093],
                               [0.0482578, 0.0],
                               [0.0328842, 0.0]]), ve_vol)
    assert_equal(len(mesh.n_flux_err[:]), len(exp_n_flux_err))
    for i in range(len(exp_n_flux_err)):
        for j in range(exp_n_flux_err.shape[1]):
            assert(isclose(mesh.n_flux_err[i][j],
                exp_n_flux_err[i][j], rel_tol=1e-5))

if __name__ == "__main__":
    nose.runmodule()

