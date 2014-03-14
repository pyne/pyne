from os import path.join, getcwd, remove
from nose.tools import assert_equal, assert_almost_equal

from pyne.r2s import irradiation_setup

thisdir = os.path.dirname(__file__)

def test_irradiation_setup():

    meshtal = join(thisdir, "files_test_r2s", "simple_meshtal")
    tally_num = 4
    cell_mats = {2: Material({2004: 1.0}, density=1.0)
                 3: Material({3007: 0.4, 3006: 0.6}, density=1.0)}
    alara_params = "Bogus line for testing" 
    geom = join(thisdir, "unitbox.h5m")
    num_rays = 10
    grid = True 
    flux_tag = "n_flux"
    fluxin = join(getcwd(), "alara_fluxin")
    reverse = False,
    alara_inp = join(getcwd(), "alara_geom")
    alara_matlib= join(getcwd(), "alara_matlib")
    output_mesh= join(getcwd(), "r2s_step1.h5m")
 
    irradiation_setup(meshtal, tally_num, cell_mats, alara_params, geom,
                      num_rays, grid, flux_tag, fluxin, reverse, alara_inp,
                      alara_matlib, output_mesh)

   remove(fluxin)
   remove(alara_inp)
   remove(alara_matlib)
   remove(output_mesh)
