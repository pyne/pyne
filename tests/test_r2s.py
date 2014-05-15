import os
import warnings
from nose.tools import assert_equal, assert_almost_equal
from numpy.testing import assert_array_equal

try:
    from itaps import iMesh
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne.r2s import irradiation_setup, photon_sampling_setup
from pyne.material import Material
from pyne.mesh import Mesh, IMeshTag
from pyne.mcnp import Meshtal

thisdir = os.path.dirname(__file__)

def test_irradiation_setup_structured():

    meshtal = os.path.join(thisdir, "files_test_r2s", "meshtal_2x2x1")
    tally_num = 4
    cell_mats = {2: Material({2004: 1.0}, density=1.0),
                 3: Material({3007: 0.4, 3006: 0.6}, density=2.0)}
    alara_params = "Bogus line for testing" 
    geom = os.path.join(thisdir, "unitbox.h5m")
    num_rays = 9
    grid = True
    flux_tag = "n_flux"
    fluxin = os.path.join(os.getcwd(), "alara_fluxin")
    reverse = True
    alara_inp = os.path.join(os.getcwd(), "alara_inp")
    alara_matlib= os.path.join(os.getcwd(), "alara_matlib")
    output_mesh= os.path.join(os.getcwd(), "r2s_step1.h5m")
 
    irradiation_setup(meshtal, cell_mats, alara_params, tally_num, geom,
                      num_rays, grid, flux_tag, fluxin, reverse, alara_inp,
                      alara_matlib, output_mesh)

    #  expected output files
    exp_alara_inp = os.path.join(thisdir, "files_test_r2s", "exp_alara_inp")
    exp_alara_matlib = os.path.join(thisdir, "files_test_r2s", 
                                             "exp_alara_matlib")
    exp_fluxin = os.path.join(thisdir, "files_test_r2s", "exp_fluxin")

    #  test alara input file
    with open(alara_inp, 'r') as f1, open(exp_alara_inp, 'r') as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    #  test alara matlibe file
    with open(alara_matlib, 'r') as f1, open(exp_alara_matlib) as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    #  test alara fluxin file
    with open(fluxin, 'r') as f1, open(exp_fluxin) as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    # test r2s step 1 output mesh
    fluxes = [[6.93088E-07, 1.04838E-06], [6.36368E-07, 9.78475E-07], 
              [5.16309E-07, 9.86586E-07], [6.36887E-07, 9.29879E-07]]
    errs = [[9.67452E-02, 7.55950E-02], [9.88806E-02, 7.61482E-02],
            [1.04090E-01, 7.69284E-02], [9.75826E-02, 7.54181E-02]]
    tot_fluxes = [1.74147E-06, 1.61484E-06, 1.50290E-06, 1.56677E-06] 
    tot_errs = [6.01522E-02, 6.13336E-02, 6.19920E-02, 5.98742E-02]

    m = Mesh(structured=True, mesh=output_mesh, mats=output_mesh)
    for i, mat, _ in m:
        assert_almost_equal(mat.density, 1.962963E+00)
        assert_equal(len(mat.comp), 3)
        assert_almost_equal(mat.comp[20040000], 1.886792E-02)
        assert_almost_equal(mat.comp[30060000], 5.886792E-01)
        assert_almost_equal(mat.comp[30070000], 3.924528E-01)
        assert_array_equal(m.n_flux[i], fluxes[i])
        assert_array_equal(m.n_flux_err[i], errs[i])
        assert_almost_equal(m.n_flux_total[i], tot_fluxes[i])
        assert_almost_equal(m.n_flux_err_total[i], tot_errs[i])

    os.remove(alara_inp)
    os.remove(alara_matlib)
    os.remove(fluxin)
    os.remove(output_mesh)


def test_photon_sampling_setup_structured():

    phtn_src = os.path.join(thisdir, "files_test_r2s", "phtn_src")
    coords = [[0, 1, 2], [0, 1, 2], [0, 1]]
    m = Mesh(structured=True, structured_coords=coords)
    tags = {(10010000, "1 h"): "tag1", ("TOTAL", "shutdown"): "tag2"}
    photon_sampling_setup(m, phtn_src, tags)

    exp_tag1 = [[1.1, 2.2], [3.3, 4.4], [5.5, 6.6], [7.7, 8.8]]
    exp_tag2 = [[11.1, 12.2], [13.3, 14.4], [15.5, 16.6], [17.7, 18.8]]

    m.tag1 = IMeshTag(2, float)
    m.tag2 = IMeshTag(2, float)

    for i, mat, ve in m:
        assert_array_equal(m.tag1[i], exp_tag1[i])
        assert_array_equal(m.tag2[i], exp_tag2[i])


def test_irradiation_setup_unstructured():

    flux_tag = "n_flux"
    meshtal_file = os.path.join(thisdir, "files_test_r2s", "meshtal_2x2x1")
    meshtal = Meshtal(meshtal_file, {4: (flux_tag, flux_tag + "_err",
                                         flux_tag + "_total",
                                         flux_tag + "_err_total")})
    #  Explicitly make this mesh unstructured, it will now iterate in yxz
    #  order which is MOAB structured mesh creation order.
    meshtal = Mesh(structured=False, mesh=meshtal.tally[4].mesh)
    meshtal_mesh_file = os.path.join(thisdir, "meshtal.h5m")
    meshtal.mesh.save(meshtal_mesh_file)

    cell_mats = {2: Material({2004: 1.0}, density=1.0),
                 3: Material({3007: 0.4, 3006: 0.6}, density=2.0)}
    alara_params = "Bogus line for testing" 
    geom = os.path.join(thisdir, "unitbox.h5m")
    flux_tag = "n_flux"
    fluxin = os.path.join(os.getcwd(), "alara_fluxin")
    reverse = True
    alara_inp = os.path.join(os.getcwd(), "alara_inp")
    alara_matlib= os.path.join(os.getcwd(), "alara_matlib")
    output_mesh= os.path.join(os.getcwd(), "r2s_step1.h5m")
 
    irradiation_setup(flux_mesh=meshtal_mesh_file, cell_mats=cell_mats, 
                      alara_params=alara_params, geom=geom, flux_tag=flux_tag, 
                      fluxin=fluxin, reverse=reverse, alara_inp=alara_inp,
                      alara_matlib=alara_matlib, output_mesh=output_mesh)

    #  expected output files
    exp_alara_inp = os.path.join(thisdir, "files_test_r2s", "exp_alara_inp")
    exp_alara_matlib = os.path.join(thisdir, "files_test_r2s", 
                                             "exp_alara_matlib_un")
    exp_fluxin = os.path.join(thisdir, "files_test_r2s", "exp_fluxin_un")

    #  test alara input file
    with open(alara_inp, 'r') as f1, open(exp_alara_inp, 'r') as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    #  test alara matlibe file
    with open(alara_matlib, 'r') as f1, open(exp_alara_matlib) as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    #  test alara fluxin file
    with open(fluxin, 'r') as f1, open(exp_fluxin) as f2:
        for a, b in zip(f1.readlines(), f2.readlines()):
            assert_equal(a, b)

    # test r2s step 1 output mesh
    fluxes = [[6.93088E-07, 1.04838E-06], [6.36368E-07, 9.78475E-07], 
              [5.16309E-07, 9.86586E-07], [6.36887E-07, 9.29879E-07]]
    errs = [[9.67452E-02, 7.55950E-02], [9.88806E-02, 7.61482E-02],
            [1.04090E-01, 7.69284E-02], [9.75826E-02, 7.54181E-02]]
    tot_fluxes = [1.74147E-06, 1.61484E-06, 1.50290E-06, 1.56677E-06] 
    tot_errs = [6.01522E-02, 6.13336E-02, 6.19920E-02, 5.98742E-02]

    m = Mesh(structured=True, mesh=output_mesh, mats=output_mesh)
    for i, mat, _ in m:
        assert_almost_equal(mat.density, 2.0)
        assert_equal(len(mat.comp), 2)
        assert_almost_equal(mat.comp[30060000], 0.6)
        assert_almost_equal(mat.comp[30070000], 0.4)
        assert_array_equal(m.n_flux[i], fluxes[i])
        assert_array_equal(m.n_flux_err[i], errs[i])
        assert_almost_equal(m.n_flux_total[i], tot_fluxes[i])
        assert_almost_equal(m.n_flux_err_total[i], tot_errs[i])

    os.remove(meshtal_mesh_file)
    os.remove(alara_inp)
    os.remove(alara_matlib)
    os.remove(fluxin)
    os.remove(output_mesh)

def test_photon_sampling_setup_unstructured():

    phtn_src = os.path.join(thisdir, "files_test_r2s", "phtn_src")
    coords = [[0, 1, 2], [0, 1, 2], [0, 1]]
    m = Mesh(structured=True, structured_coords=coords)
    m.structured = False
    tags = {(10010000, "1 h"): "tag1", ("TOTAL", "shutdown"): "tag2"}
    photon_sampling_setup(m, phtn_src, tags)

    exp_tag1 = [[1.1, 2.2], [3.3, 4.4], [5.5, 6.6], [7.7, 8.8]]
    exp_tag2 = [[11.1, 12.2], [13.3, 14.4], [15.5, 16.6], [17.7, 18.8]]

    m.tag1 = IMeshTag(2, float)
    m.tag2 = IMeshTag(2, float)

    for i, mat, ve in m:
        assert_array_equal(m.tag1[i], exp_tag1[i])
        assert_array_equal(m.tag2[i], exp_tag2[i])
