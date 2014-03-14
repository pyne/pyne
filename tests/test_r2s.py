import os
from nose.tools import assert_equal, assert_almost_equal

from pyne.r2s import irradiation_setup
from pyne.material import Material
from pyne.mesh import Mesh, IMeshTag

thisdir = os.path.dirname(__file__)

def test_irradiation_setup():

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
    reverse = False,
    alara_inp = os.path.join(os.getcwd(), "alara_inp")
    alara_matlib= os.path.join(os.getcwd(), "alara_matlib")
    output_mesh= os.path.join(os.getcwd(), "r2s_step1.h5m")
 
    irradiation_setup(meshtal, tally_num, cell_mats, alara_params, geom,
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

    m = Mesh(structured=True, mesh_file=output_mesh)

    print m.mats
    m.flux = IMeshTag(2, float)
    m.err = IMeshTag(2, float)
    m.flux_tot = IMeshTag(1, float)
    m.flux_err_tot = IMeshTag(1, float)

    for i, mat, _ in m:
        assert_almost_equal(mat.density, 1.962963E+00)
        assert_almost_equal(mat.comp, {20040000: 1.886792E-02, 
                                       30060000: 5.886792E-01, 
                                       20070000: 3.924528E-01})
        assert_array_equal(m.flux[i], fluxes[i])
        assert_array_equal(m.err[i], errs[i])
        assert_almost_equal(m.flux_tot[i], tot_fluxes[i])
        assert_almost_equal(m.flux_err_tot[i], tot_errs[i])

    os.remove(alara_inp)
    os.remove(alara_matlib)
    os.remove(fluxin)
