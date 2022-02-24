"""openmc_utils tests"""
from __future__ import unicode_literals, division
from io import StringIO
import warnings
import os
import sys
import nose
from nose.plugins.skip import SkipTest
import numpy as np
from nose.tools import assert_equal, assert_true
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning, is_close

warnings.simplefilter("ignore", QAWarning)

from pyne import nucname
from pyne import openmc_utils
from pyne.mesh import HAVE_PYMOAB, Mesh, MeshTally

sample_xs = StringIO(
    """<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
</cross_sections>
"""
)

sample_xs_with_dir = StringIO(
    """<?xml version="1.0" ?>
<cross_sections>
  <directory>/</directory>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
"""
)

sample_xs_with_mcnp_id = StringIO(
    """<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="Co-58m.70c" awr="57.4381" location="28897" metastable="1" name="27458.70c" path="endf70b" temperature="2.5301e-08" zaid="27458"/>
  <ace_table alias="Co-58.70c" awr="57.4381" location="28699" name="27058.70c" path="endf70b" temperature="2.5301e-08" zaid="27058"/>
</cross_sections>
"""
)

cwd = os.getcwd()


def test_ace_table_init():
    atab = openmc_utils.AceTable(
        zaid="92235", path="U235.ace", cross_sections_path="/tmp/cross_sections.xml"
    )
    assert_equal("92235", atab.zaid)
    assert_equal("U235.ace", atab.path)
    assert_equal("/tmp/U235.ace", atab.abspath)
    assert_equal(nucname.id("U235"), atab.nucid)


def test_ace_table_xml():
    atab = openmc_utils.AceTable(
        zaid="92235", path="U235.ace", cross_sections_path="/tmp/cross_sections.xml"
    )
    exp = '<ace_table path="U235.ace" zaid="92235"/>'
    obs = atab.xml()
    assert_equal(exp, obs)


def test_cross_sections_read():
    sample_xs.seek(0)
    xs = openmc_utils.CrossSections(sample_xs)
    assert_equal("ascii", xs.filetype)
    assert_true(xs.path is None)

    exp = [
        openmc_utils.AceTable(
            alias="H-1.71c",
            awr="0.999167",
            location="1",
            name="1001.71c",
            path="293.6K/H_001_293.6K.ace",
            temperature="2.53e-08",
            zaid="1001",
        ),
        openmc_utils.AceTable(
            alias="Am-242m.73c",
            awr="239.9801",
            location="1",
            metastable="1",
            name="95242.73c",
            path="900K/Am_242_900K.ace",
            temperature="7.756e-08",
            zaid="95242",
        ),
        openmc_utils.AceTable(
            awr="89.1324",
            location="1",
            name="ZrZrH.71t",
            path="tsl/zrzrh.acer",
            temperature="2.551e-08",
            zaid="0",
        ),
    ]
    assert_equal(exp, xs.ace_tables)


def test_cross_sections_abspath_with_dir():
    xs = openmc_utils.CrossSections(sample_xs_with_dir)
    assert_equal("ascii", xs.filetype)
    assert_equal(xs.path, "/")

    exp_abspaths = [
        "/293.6K/H_001_293.6K.ace",
        "/900K/Am_242_900K.ace",
        "/tsl/zrzrh.acer",
    ]
    obs_abspaths = [table.abspath for table in xs.ace_tables]
    assert_equal(exp_abspaths, obs_abspaths)


def test_cross_sections_mcnp_id():
    xstables = openmc_utils.CrossSections(sample_xs_with_mcnp_id).ace_tables
    mcnp_obs = [table.nucid for table in xstables if table.alias == "Co-58m.70c"][0]
    assert_equal(mcnp_obs, 270580001)
    nucid_obs = [table.nucid for table in xstables if table.alias == "Co-58.70c"][0]
    assert_equal(nucid_obs, 270580000)


def test_cross_sections_roundtrip():
    sample_xs.seek(0)
    xs = openmc_utils.CrossSections(sample_xs)
    sample_xs.seek(0)
    exp = sample_xs.read()
    obs = xs.xml()
    assert_equal(exp, obs)


def test_calc_structured_coords():
    lower_left = np.array([0.0, 0.0, 0.0])
    upper_right = np.array([1.0, 2.0, 3.0])
    dimension = np.array([2, 4, 5])
    exp_structured_coords = np.array(
        [[0.0, 0.5, 1.0], [0.0, 0.5, 1.0, 1.5, 2.0], [0.0, 0.6, 1.2, 1.8, 2.4, 3.0]]
    )
    structured_coords = openmc_utils.calc_structured_coords(
        lower_left, upper_right, dimension
    )
    assert_equal(len(structured_coords), len(exp_structured_coords))
    for i in range(len(exp_structured_coords)):
        assert_array_almost_equal(structured_coords[i], exp_structured_coords[i])


def test_get_e_bounds_from_openmc_sp():
    if not HAVE_PYMOAB or sys.version_info[0] == 2:
        raise SkipTest
    try:
        import openmc
    except:
        raise SkipTest
    # energy bin: [0.0, 1.0, 20.0], 2bins
    # 6 volume elemenes
    filename = os.path.join(cwd, "files_test_openmc", "statepoint.10.ebin2.ves6.h5")
    # OpenMC energy unit is eV
    exp_e_bounds = np.array([0.0, 1.0, 20.0]) * 1e6
    e_bounds = openmc_utils.get_e_bounds_from_openmc_sp(filename, tally_id=1)
    assert_array_equal(e_bounds, exp_e_bounds)


def test_flux_changes_order():
    # test mesh dimensions: (3, 2, 1), num_e_groups: 2
    dims = np.array([3, 2, 1])
    result = np.array(
        [
            ["x0y0z0e0", "x0y0z0e1"],
            ["x1y0z0e0", "x1y0z0e1"],
            ["x2y0z0e0", "x2y0z0e1"],
            ["x0y1z0e0", "x0y1z0e1"],
            ["x1y1z0e0", "x1y1z0e1"],
            ["x2y1z0e0", "x2y1z0e1"],
        ]
    )
    exp_result = np.array(
        [
            ["x0y0z0e0", "x0y0z0e1"],
            ["x0y1z0e0", "x0y1z0e1"],
            ["x1y0z0e0", "x1y0z0e1"],
            ["x1y1z0e0", "x1y1z0e1"],
            ["x2y0z0e0", "x2y0z0e1"],
            ["x2y1z0e0", "x2y1z0e1"],
        ]
    )
    result_ = openmc_utils.result_changes_order(result, dims)
    assert_array_equal(result_, exp_result)


def test_get_result_error_from_openmc_sp():
    if not HAVE_PYMOAB or sys.version_info[0] == 2:
        raise SkipTest
    try:
        import openmc
    except:
        raise SkipTest
    filename = os.path.join(
        os.getcwd(), "files_test_openmc", "statepoint.10.ebin2.ves6.h5"
    )
    tally_num = 1
    dims = np.array([3, 2, 1])

    m = MeshTally()
    structured_coords = openmc_utils.get_structured_coords_from_openmc_sp(
        filename, tally_id=tally_num
    )

    super(MeshTally, m).__init__(
        structured_coords=structured_coords, structured=True, mats=()
    )
    num_ves = len(m)
    ve_vol = m.structured_hex_volume(0, 0, 0)

    (
        result,
        rel_err,
        res_tot,
        rel_err_tot,
    ) = openmc_utils.get_result_error_from_openmc_sp(filename, m)
    # read expected data from statepoint.10.ebin2.ves6.h5
    sp = openmc.StatePoint(filename)
    tally = sp.get_tally(id=1)
    flux = tally.get_slice(scores=["flux"])
    num_e_groups = len(flux.mean.flatten()) // num_ves
    exp_result = np.divide(flux.mean.flatten(), ve_vol)
    exp_result = np.reshape(exp_result, newshape=(num_e_groups, num_ves))
    exp_result = exp_result.transpose()
    exp_res_tot = np.sum(exp_result, axis=1)

    exp_rel_err = flux.std_dev / flux.mean
    exp_rel_err = np.reshape(exp_rel_err, newshape=(num_e_groups, num_ves))
    exp_rel_err = exp_rel_err.transpose()

    exp_rel_err_tot = np.zeros_like(exp_res_tot)
    std_dev = np.reshape(flux.std_dev.flatten(), newshape=(num_e_groups, num_ves))
    std_dev = std_dev.transpose()
    var_tot = np.sum(np.square(std_dev), axis=1)
    exp_rel_err_tot = np.sqrt(var_tot) / (exp_res_tot * ve_vol)
    # changes the order of exp results
    exp_result = openmc_utils.result_changes_order(exp_result, dims)
    exp_rel_err = openmc_utils.result_changes_order(exp_rel_err, dims)
    exp_res_tot = openmc_utils.result_changes_order(exp_res_tot, dims)
    exp_rel_err_tot = openmc_utils.result_changes_order(exp_rel_err_tot, dims)
    # compare the data and expected answer
    assert_array_almost_equal(result, exp_result)
    assert_array_almost_equal(rel_err, exp_rel_err)
    assert_array_almost_equal(res_tot, exp_res_tot)
    assert_array_almost_equal(rel_err_tot, exp_rel_err_tot)


def test_create_meshtally():
    if not HAVE_PYMOAB or sys.version_info[0] == 2:
        raise SkipTest
    try:
        import openmc
    except:
        raise SkipTest
    # mesh read from openmc state point file
    # Parameters of the tally and mesh
    # mesh = openmc_utils.Mesh(mesh_id=1, name="n_flux")
    # mesh.dimension= [3, 2, 1]
    # mesh.lower_left = (-40.0, -12.5, -2.5)
    # mesh.upper_right = (40.0, 12.5, 2.5)
    # energy_bins = np.array([0.0, 1.0, 20.0]) * 1e6
    filename = os.path.join(
        os.getcwd(), "files_test_openmc", "statepoint.10.ebin2.ves6.h5"
    )
    tally_num = 1
    dims = np.array([3, 2, 1])
    tag_names = ("n_flux", "n_flux_err", "n_flux_total", "n_flux_total_err")
    mesh = openmc_utils.create_meshtally(
        filename, tally_num, particle="neutron", tag_names=tag_names
    )
    num_ves = len(mesh)
    # check mesh attributes
    assert_equal(num_ves, 6)
    assert mesh.structured
    # structured_coords
    assert_array_almost_equal(
        mesh.structured_coords[0], [(-40.0 + x * 80.0 / 3) for x in range(0, 4)]
    )
    assert_array_almost_equal(
        mesh.structured_coords[1], [(-12.5 + x * 25.0 / 2) for x in range(0, 3)]
    )
    assert_array_almost_equal(
        mesh.structured_coords[2], [(-2.5 + x * 5.0 / 1) for x in range(0, 2)]
    )
    ve_vol = (80.0 / 3) * (25.0 / 2) * (5.0 / 1)

    # read expected data from statepoint.10.ebin2.ves6.h5
    sp = openmc.StatePoint(filename)
    tally = sp.get_tally(id=1)
    flux = tally.get_slice(scores=["flux"])
    num_e_groups = len(flux.mean.flatten()) // num_ves
    exp_result = np.divide(flux.mean.flatten(), ve_vol)
    exp_result = np.reshape(exp_result, newshape=(num_e_groups, num_ves))
    exp_result = exp_result.transpose()
    exp_rel_err = flux.std_dev / flux.mean
    exp_rel_err = np.reshape(exp_rel_err, newshape=(num_e_groups, num_ves))
    exp_rel_err = exp_rel_err.transpose()

    # changes the order of exp results
    exp_result = openmc_utils.result_changes_order(exp_result, dims)
    exp_rel_err = openmc_utils.result_changes_order(exp_rel_err, dims)
    # compare
    assert_array_almost_equal(mesh.n_flux[:], exp_result)
    assert_array_almost_equal(mesh.n_flux_err[:], exp_rel_err)


if __name__ == "__main__":
    nose.runmodule()
