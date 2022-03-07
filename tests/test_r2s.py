import os
import warnings
from nose.tools import assert_equal, assert_almost_equal
from nose.plugins.skip import SkipTest
import numpy as np
import tables as tb
from numpy.testing import assert_array_equal
import multiprocessing
import filecmp
import sys
from shutil import copyfile

from pyne.mcnp import Meshtal
from pyne.material import Material
from pyne.r2s import (
    irradiation_setup,
    photon_sampling_setup,
    total_photon_source_intensity,
    tag_e_bounds,
    tag_source_intensity,
    tag_decay_time,
    tag_version,
)
from pyne.utils import QAWarning, file_almost_same, file_block_almost_same
from pyne.mesh import Mesh, NativeMeshTag, HAVE_PYMOAB

if not HAVE_PYMOAB:
    raise SkipTest

if sys.version_info[0] > 2:
    izip = zip
else:
    from itertools import izip

warnings.simplefilter("ignore", QAWarning)

thisdir = os.path.dirname(__file__)


def irradiation_setup_structured(flux_tag="n_flux", meshtal_file="meshtal_2x2x1"):

    meshtal = os.path.join(thisdir, "files_test_r2s", meshtal_file)
    tally_num = 4
    cell_mats = {
        2: Material({2004: 1.0}, density=1.0, metadata={"name": "mat_11"}),
        3: Material({3007: 0.4, 3006: 0.6}, density=2.0, metadata={"name": "mat_12"}),
    }
    alara_params = "Bogus line for testing\n"
    geom = os.path.join(thisdir, "unitbox.h5m")
    num_rays = 9
    grid = True
    fluxin = os.path.join(os.getcwd(), "alara_fluxin")
    reverse = True
    alara_inp = os.path.join(os.getcwd(), "alara_inp")
    alara_matlib = os.path.join(os.getcwd(), "alara_matlib")
    output_mesh = os.path.join(os.getcwd(), "r2s_step1.h5m")
    output_material = True
    cell_fracs = np.zeros(
        8,
        dtype=[
            ("idx", np.int64),
            ("cell", np.int64),
            ("vol_frac", np.float64),
            ("rel_error", np.float64),
        ],
    )
    cell_fracs[:] = [
        (0, 2, 0.037037037037037035, 0.5443310539518174),
        (0, 3, 0.9629629629629631, 0.010467904883688123),
        (1, 2, 0.037037037037037035, 0.5443310539518174),
        (1, 3, 0.9629629629629629, 0.010467904883688454),
        (2, 2, 0.037037037037037035, 0.5443310539518174),
        (2, 3, 0.9629629629629629, 0.010467904883688454),
        (3, 2, 0.037037037037037035, 0.5443310539518174),
        (3, 3, 0.9629629629629629, 0.010467904883688454),
    ]

    irradiation_setup(
        meshtal,
        cell_mats,
        cell_fracs,
        alara_params,
        tally_num,
        num_rays,
        grid,
        flux_tag,
        fluxin,
        reverse,
        alara_inp,
        alara_matlib,
        output_mesh,
        output_material,
    )

    #  expected output files
    exp_alara_inp = os.path.join(thisdir, "files_test_r2s", "exp_alara_inp")
    exp_alara_matlib = os.path.join(thisdir, "files_test_r2s", "exp_alara_matlib")
    exp_fluxin = os.path.join(thisdir, "files_test_r2s", "exp_fluxin")

    # test files
    f1 = filecmp.cmp(alara_inp, exp_alara_inp)
    f2 = filecmp.cmp(alara_matlib, exp_alara_matlib)
    f3 = filecmp.cmp(fluxin, exp_fluxin)

    m = Mesh(structured=True, mesh=output_mesh, mats=output_mesh)

    out = [
        m.n_flux[:].tolist(),
        m.n_flux_err[:].tolist(),
        m.n_flux_total[:].tolist(),
        m.n_flux_err_total[:].tolist(),
        [x.comp.items() for y, x, z in m],
        [x.density for y, x, z in m],
    ]

    n_flux = out[0]
    n_flux_err = out[1]
    n_flux_total = out[2]
    n_flux_err_total = out[3]
    densities = out[5]

    comps = np.zeros(shape=(len(out[4])), dtype=dict)
    for i, comp in enumerate(out[4]):
        comps[i] = {}
        for nucid in comp:
            comps[i][nucid[0]] = nucid[1]

    # test r2s step 1 output mesh
    fluxes = [
        [6.93088e-07, 1.04838e-06],
        [6.36368e-07, 9.78475e-07],
        [5.16309e-07, 9.86586e-07],
        [6.36887e-07, 9.29879e-07],
    ]
    errs = [
        [9.67452e-02, 7.55950e-02],
        [9.88806e-02, 7.61482e-02],
        [1.04090e-01, 7.69284e-02],
        [9.75826e-02, 7.54181e-02],
    ]
    tot_fluxes = [1.74147e-06, 1.61484e-06, 1.50290e-06, 1.56677e-06]
    tot_errs = [6.01522e-02, 6.13336e-02, 6.19920e-02, 5.98742e-02]

    i = 0
    for nf, nfe, nft, nfte, comp, density in izip(
        n_flux, n_flux_err, n_flux_total, n_flux_err_total, comps, densities
    ):
        assert_almost_equal(density, 1.962963e00)
        assert_equal(len(comp), 3)
        assert_almost_equal(comp[20040000], 1.886792e-02)
        assert_almost_equal(comp[30060000], 5.886792e-01)
        assert_almost_equal(comp[30070000], 3.924528e-01)
        assert_array_equal(nf, fluxes[i])
        assert_array_equal(nfe, errs[i])
        assert_almost_equal(nft, tot_fluxes[i])
        assert_almost_equal(nfte, tot_errs[i])
        i += 1

    os.remove(alara_inp)
    os.remove(alara_matlib)
    os.remove(fluxin)
    os.remove(output_mesh)

    return [f1, f2, f3]


def test_irradiation_setup_structured():
    p = multiprocessing.Pool()
    r = p.apply_async(irradiation_setup_structured)
    p.close()
    p.join()
    results = r.get()

    # unpack return values
    f1 = results[0]
    f2 = results[1]
    f3 = results[2]

    # test that files match
    assert f1 is True
    assert f2 is True
    assert f3 is True


def test_photon_sampling_setup_structured():

    phtn_src = os.path.join(thisdir, "files_test_r2s", "phtn_src")
    coords = [[0, 1, 2], [0, 1, 2], [0, 1]]
    m = Mesh(structured=True, structured_coords=coords)
    tags = {(10010000, "1 h"): "tag1", ("TOTAL", "shutdown"): "tag2"}
    photon_sampling_setup(m, phtn_src, tags)
    os.remove(phtn_src + ".h5")

    exp_tag1 = [[1.1, 2.2], [3.3, 4.4], [5.5, 6.6], [7.7, 8.8]]
    exp_tag2 = [[11.1, 12.2], [13.3, 14.4], [15.5, 16.6], [17.7, 18.8]]

    m.tag1 = NativeMeshTag(2, float)
    m.tag2 = NativeMeshTag(2, float)

    for i, mat, ve in m:
        assert_array_equal(m.tag1[i], exp_tag1[i])
        assert_array_equal(m.tag2[i], exp_tag2[i])


def irradiation_setup_unstructured(flux_tag="n_flux"):
    meshtal_filename = "meshtal_2x2x1"
    meshtal_file = os.path.join(thisdir, "files_test_r2s", meshtal_filename)

    meshtal = Meshtal(
        meshtal_file,
        {
            4: (
                flux_tag,
                flux_tag + "_err",
                flux_tag + "_total",
                flux_tag + "_err_total",
            )
        },
    )
    #  Explicitly make this mesh unstructured, it will now iterate in yxz
    #  order which is MOAB structured mesh creation order.
    meshtal = Mesh(structured=False, mesh=meshtal.tally[4].mesh)
    meshtal_mesh_file = os.path.join(thisdir, "meshtal.h5m")
    meshtal.write_hdf5(meshtal_mesh_file, write_mats=False)

    if flux_tag != "n_flux":
        # if not using n_flux makes a mesh containing n_flux tag, and then
        # makes a new tag called flux_tag, to use later in the test
        flux_tag_name = "n_flux"
        meshtal = Meshtal(
            meshtal_file,
            {
                4: (
                    flux_tag_name,
                    flux_tag_name + "_err",
                    flux_tag_name + "_total",
                    flux_tag_name + "_err_total",
                )
            },
        )
        #  Explicitly make this mesh unstructured, it will now iterate in yxz
        #  order which is MOAB structured mesh creation order.
        meshtal = Mesh(structured=False, mesh=meshtal.tally[4].mesh)
        meshtal_mesh_file = os.path.join(thisdir, "meshtal.h5m")
        meshtal.write_hdf5(meshtal_mesh_file, write_mats=False)
        new_mesh = Mesh(structured=False, mesh=meshtal_mesh_file)
        new_mesh.TALLY_TAG = NativeMeshTag(2, float)  # 2 egroups
        new_mesh.TALLY_TAG = meshtal.n_flux[:]

        # overwrite the mesh file
        new_mesh.write_hdf5(meshtal_mesh_file, write_mats=False)

    cell_mats = {
        2: Material({2004: 1.0}, density=1.0, metadata={"name": "mat_11"}),
        3: Material({3007: 0.4, 3006: 0.6}, density=2.0, metadata={"name": "mat_12"}),
    }
    alara_params = "Bogus line for testing\n"
    geom = os.path.join(thisdir, "unitbox.h5m")
    fluxin = os.path.join(os.getcwd(), "alara_fluxin")
    reverse = True
    alara_inp = os.path.join(os.getcwd(), "alara_inp")
    alara_matlib = os.path.join(os.getcwd(), "alara_matlib")
    output_mesh = os.path.join(os.getcwd(), "r2s_step1.h5m")
    output_material = True
    cell_fracs = np.zeros(
        4,
        dtype=[
            ("idx", np.int64),
            ("cell", np.int64),
            ("vol_frac", np.float64),
            ("rel_error", np.float64),
        ],
    )
    cell_fracs[:] = [
        (0, 3, 1.0, 1.0),
        (1, 3, 1.0, 1.0),
        (2, 3, 1.0, 1.0),
        (3, 3, 1.0, 1.0),
    ]
    irradiation_setup(
        flux_mesh=meshtal_mesh_file,
        cell_mats=cell_mats,
        cell_fracs=cell_fracs,
        alara_params=alara_params,
        flux_tag=flux_tag,
        fluxin=fluxin,
        reverse=reverse,
        alara_inp=alara_inp,
        alara_matlib=alara_matlib,
        output_mesh=output_mesh,
        output_material=output_material,
    )

    #  expected output files
    exp_alara_inp = os.path.join(thisdir, "files_test_r2s", "exp_alara_inp_un")
    exp_alara_matlib = os.path.join(thisdir, "files_test_r2s", "exp_alara_matlib")
    exp_fluxin = os.path.join(thisdir, "files_test_r2s", "exp_fluxin_un")

    # test files
    f1 = filecmp.cmp(alara_inp, exp_alara_inp)
    f2 = filecmp.cmp(alara_matlib, exp_alara_matlib)
    f3 = filecmp.cmp(fluxin, exp_fluxin)

    m = Mesh(structured=True, mesh=output_mesh, mats=output_mesh)

    out = [
        m.n_flux[:].tolist(),
        m.n_flux_err[:].tolist(),
        m.n_flux_total[:].tolist(),
        m.n_flux_err_total[:].tolist(),
        [x.comp.items() for y, x, z in m],
        [x.density for y, x, z in m],
    ]

    n_flux = out[0]
    n_flux_err = out[1]
    n_flux_total = out[2]
    n_flux_err_total = out[3]
    densities = out[5]

    comps = np.zeros(shape=(len(out[4])), dtype=dict)
    for i, comp in enumerate(out[4]):
        comps[i] = {}
        for nucid in comp:
            comps[i][nucid[0]] = nucid[1]

    # test r2s step 1 output mesh
    fluxes = [
        [6.93088e-07, 1.04838e-06],
        [6.36368e-07, 9.78475e-07],
        [5.16309e-07, 9.86586e-07],
        [6.36887e-07, 9.29879e-07],
    ]
    errs = [
        [9.67452e-02, 7.55950e-02],
        [9.88806e-02, 7.61482e-02],
        [1.04090e-01, 7.69284e-02],
        [9.75826e-02, 7.54181e-02],
    ]
    tot_fluxes = [1.74147e-06, 1.61484e-06, 1.50290e-06, 1.56677e-06]
    tot_errs = [6.01522e-02, 6.13336e-02, 6.19920e-02, 5.98742e-02]

    i = 0
    for nf, nfe, nft, nfte, comp, density in izip(
        n_flux, n_flux_err, n_flux_total, n_flux_err_total, comps, densities
    ):
        assert_almost_equal(density, 2.0)
        assert_equal(len(comp), 2)
        assert_almost_equal(comp[30060000], 0.6)
        assert_almost_equal(comp[30070000], 0.4)
        assert_array_equal(nf, fluxes[i])
        assert_array_equal(nfe, errs[i])
        assert_almost_equal(nft, tot_fluxes[i])
        assert_almost_equal(nfte, tot_errs[i])
        i += 1

    os.remove(meshtal_mesh_file)
    os.remove(alara_inp)
    os.remove(alara_matlib)
    os.remove(fluxin)
    os.remove(output_mesh)

    return [f1, f2, f3]


def test_irradiation_setup_unstructured():

    # make new file with non default tag
    p = multiprocessing.Pool()
    r = p.apply_async(irradiation_setup_unstructured)
    p.close()
    p.join()
    results = r.get()

    # unpack return values
    f1 = results[0]
    f2 = results[1]
    f3 = results[2]

    # test that files match
    assert f1 is True
    assert f2 is True
    assert f3 is True


def test_photon_sampling_setup_unstructured():

    phtn_src = os.path.join(thisdir, "files_test_r2s", "phtn_src")
    coords = [[0, 1, 2], [0, 1, 2], [0, 1]]
    m = Mesh(structured=True, structured_coords=coords)
    m.structured = False
    tags = {(10010000, "1 h"): "tag1", ("TOTAL", "shutdown"): "tag2"}
    photon_sampling_setup(m, phtn_src, tags)
    os.remove(phtn_src + ".h5")

    exp_tag1 = [[1.1, 2.2], [3.3, 4.4], [5.5, 6.6], [7.7, 8.8]]
    exp_tag2 = [[11.1, 12.2], [13.3, 14.4], [15.5, 16.6], [17.7, 18.8]]

    m.tag1 = NativeMeshTag(2, float)
    m.tag2 = NativeMeshTag(2, float)

    for i, mat, ve in m:
        assert_array_equal(m.tag1[i], exp_tag1[i])
        assert_array_equal(m.tag2[i], exp_tag2[i])


def test_total_photon_source_intensity():

    m = Mesh(structured=True, structured_coords=[[0, 1, 2], [0, 1, 3], [0, 1]])
    m.source_density = NativeMeshTag(2, float)
    m.source_density[:] = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]

    intensity = total_photon_source_intensity(m, "source_density")
    assert_equal(intensity, 58)


def test_total_photon_source_intensity_subvoxel():
    # In the calculation of the total photon source intensities under subvoxel
    # mode, the volume fractions of each subvoxel should be multiplied

    # Set up 4 voxels with the volume of: 1, 2, 1, 2
    m = Mesh(structured=True, structured_coords=[[0, 1, 2], [0, 1, 3], [0, 1]])
    # 4 voxels, each voxel contains two subvoxels -> 8 subvoxels
    # The volume fraction of each subvoxel is 0.5
    cell_fracs = np.zeros(
        8,
        dtype=[
            ("idx", np.int64),
            ("cell", np.int64),
            ("vol_frac", np.float64),
            ("rel_error", np.float64),
        ],
    )
    cell_fracs[:] = [
        (0, 11, 0.5, 0.0),
        (0, 12, 0.5, 0.0),
        (1, 11, 0.5, 0.0),
        (1, 12, 0.5, 0.0),
        (2, 13, 0.5, 0.0),
        (2, 11, 0.5, 0.0),
        (3, 12, 0.5, 0.0),
        (3, 13, 0.5, 0.0),
    ]
    m.tag_cell_fracs(cell_fracs)
    # Set up the source density with energy group number of 2
    m.source_density = NativeMeshTag(4, float)
    m.source_density[:] = [
        [0.0, 0.0, 1.0, 1.0],
        [2.0, 2.0, 3.0, 3.0],
        [4.0, 4.0, 5.0, 5.0],
        [6.0, 6.0, 7.0, 7.0],
    ]
    intensity = total_photon_source_intensity(m, "source_density", True)
    # expected intensity: each line represents a mesh voxel
    # for each subvoxel: voxel_vol * cell_fracs * photon_intensity
    expected_intensity = 1 * 0.5 * (0.0 + 0.0) + 1 * 0.5 * (1.0 + 1.0)
    expected_intensity += 2 * 0.5 * (2.0 + 2.0) + 2 * 0.5 * (3.0 + 3.0)
    expected_intensity += 1 * 0.5 * (4.0 + 4.0) + 1 * 0.5 * (5.0 + 5.0)
    expected_intensity += 2 * 0.5 * (6.0 + 6.0) + 2 * 0.5 * (7.0 + 7.0)
    assert_equal(intensity, expected_intensity)


def test_irradiation_setup_unstructured_nondef_tag():
    p = multiprocessing.Pool()
    r = p.apply_async(irradiation_setup_unstructured, ("TALLY_TAG",))
    p.close()
    p.join()
    results = r.get()

    # unpack return values
    f1 = results[0]
    f2 = results[1]
    f3 = results[2]


def test_tag_e_bounds():
    m = Mesh(structured=True, structured_coords=[[0, 1, 2], [0, 1, 3], [0, 1]])
    e_bounds = np.array([0.0, 0.1, 20])
    m = tag_e_bounds(m, e_bounds)
    assert_array_equal(m.e_bounds[m], e_bounds)


def test_tag_source_intensity():
    m = Mesh(structured=True, structured_coords=[[0, 1, 2], [0, 1, 3], [0, 1]])
    source_intensity = 1.0
    m = tag_source_intensity(m, source_intensity)
    assert_array_equal(m.source_intensity[m], source_intensity)


def test_tag_decay_time():
    m = Mesh(structured=True, structured_coords=[[0, 1, 2], [0, 1, 3], [0, 1]])
    decay_time = 1.0
    m = tag_decay_time(m, decay_time)
    assert_array_equal(m.decay_time[m], decay_time)


def _r2s_test_step1(r2s_run_dir, remove_step1_out=True):
    os.chdir(thisdir)
    # copy ../scripts/r2s.py to r2s_run_dir/r2s.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(r2s_run_dir, "r2s.py")
    copyfile(os.path.join(folderpath, "scripts", "r2s.py"), dst)

    # run r2s step1
    os.chdir(r2s_run_dir)
    os.system("python r2s.py step1")

    # output files of r2s step1
    alara_inp = os.path.join(r2s_run_dir, "alara_inp")
    alara_matlib = os.path.join(r2s_run_dir, "alara_matlib")
    alara_fluxin = os.path.join(r2s_run_dir, "alara_fluxin")
    blank_mesh = os.path.join(r2s_run_dir, "blank_mesh.h5m")
    step1_file = os.path.join(r2s_run_dir, "r2s_step1.h5m")

    exp_alara_inp = os.path.join(r2s_run_dir, "exp_alara_inp")
    exp_alara_matlib = os.path.join(r2s_run_dir, "exp_alara_matlib")
    exp_alara_fluxin = os.path.join(r2s_run_dir, "exp_alara_fluxin")

    # compare the output file of step1
    f1 = filecmp.cmp(alara_inp, exp_alara_inp)
    f2 = file_block_almost_same(alara_matlib, exp_alara_matlib)
    f3 = filecmp.cmp(alara_fluxin, exp_alara_fluxin)

    # remove test output files
    if remove_step1_out:
        # these files are used in step2
        os.remove(blank_mesh)
        os.remove(alara_inp)
    # these files are not used in step2
    os.remove(alara_fluxin)
    os.remove(alara_matlib)
    os.remove(step1_file)
    os.remove(dst)

    assert_equal(f1, True)
    assert_equal(f2, True)
    assert_equal(f3, True)


def _r2s_test_step2(r2s_run_dir, remove_step1_out=True):
    os.chdir(thisdir)
    # copy ../scripts/r2s.py to r2s_run_dir/r2s.py
    os.chdir("..")
    folderpath = os.getcwd()
    dst = os.path.join(r2s_run_dir, "r2s.py")
    copyfile(os.path.join(folderpath, "scripts", "r2s.py"), dst)

    # output files of r2s step1
    alara_inp = os.path.join(r2s_run_dir, "alara_inp")
    blank_mesh = os.path.join(r2s_run_dir, "blank_mesh.h5m")
    if remove_step1_out:
        copyfile(os.path.join(r2s_run_dir, "exp_alara_inp"), alara_inp)
        copyfile(os.path.join(r2s_run_dir, "exp_blank_mesh.h5m"), blank_mesh)

    # run r2s step2
    os.chdir(r2s_run_dir)
    os.system("python r2s.py step2")
    os.remove(blank_mesh)
    os.remove(alara_inp)
    os.remove(dst)

    # output files of r2s step2
    p_src = os.path.join(r2s_run_dir, "phtn_src.h5")
    os.remove(p_src)

    # compare the total photon source intensities
    t_p_src = os.path.join(r2s_run_dir, "total_photon_source_intensities.txt")
    exp_t_p_src = os.path.join(r2s_run_dir, "exp_total_photon_source_intensities.txt")
    f5 = file_almost_same(t_p_src, exp_t_p_src)
    assert_equal(f5, True)
    os.remove(t_p_src)

    # compare the results
    src_c1 = os.path.join(r2s_run_dir, "source_1.h5m")
    exp_src_c1 = os.path.join(r2s_run_dir, "exp_source_1.h5m")
    # skip test if h5diff not exist
    if "unstructured" in r2s_run_dir:
        ele_type = "Tet4"
    else:
        ele_type = "Hex8"
    is_h5diff = os.system("which h5diff")
    if is_h5diff == 0:
        # compare two h5 files
        f6 = True  # compre source_density
        command = "".join(
            [
                "h5diff --relative=1e-6 ",
                src_c1,
                " ",
                exp_src_c1,
                " /tstt/elements/",
                ele_type,
                "/tags/source_density",
                " /tstt/elements/",
                ele_type,
                "/tags/source_density",
            ]
        )
        diff_flag = os.system(command)
        # return value 0 if no difference, 1 if differences found, 2 if error
        f6 = True if diff_flag == 0 else False
        assert_equal(f6, True)

        # compare e_bounds
        f4 = True
        command = "".join(
            [
                "h5diff --relative=1e-6 ",
                src_c1,
                " ",
                exp_src_c1,
                " /tstt/tags/e_bounds",
                " /tstt/tags/e_bounds",
            ]
        )
        diff_flag = os.system(command)
        f4 = True if diff_flag == 0 else False
        assert_equal(f4, True)
        # compare decay_time
        f7 = True
        command = "".join(
            [
                "h5diff --relative=1e-6 ",
                src_c1,
                " ",
                exp_src_c1,
                " /tstt/tags/decay_time",
                " /tstt/tags/decay_time",
            ]
        )
        diff_flag = os.system(command)
        f7 = True if diff_flag == 0 else False
        assert_equal(f7, True)
        # compare total_photon_source_intensity
        f8 = True
        command = "".join(
            [
                "h5diff --relative=1e-6 ",
                src_c1,
                " ",
                exp_src_c1,
                " /tstt/tags/source_intensity",
                " /tstt/tags/source_intensity",
            ]
        )
        diff_flag = os.system(command)
        f8 = True if diff_flag == 0 else False
        assert_equal(f8, True)
        # compare r2s soruce file version
        f9 = True
        command = "".join(
            [
                "h5diff ",
                src_c1,
                " ",
                exp_src_c1,
                " /tstt/tags/r2s_source_file_version",
                " /tstt/tags/r2s_source_file_version",
            ]
        )
        diff_flag = os.system(command)
        f9 = True if diff_flag == 0 else False
        assert_equal(f9, True)

    os.remove(src_c1)


def test_r2s_script_step_by_step():
    # skip test without dagmc
    try:
        from pyne import dagmc
    except ImportError:
        raise SkipTest

    remove_step1_out = True
    r2s_run_dir = os.path.join(thisdir, "files_test_r2s", "r2s_examples", "r2s_run")
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test sub-voxel r2s
    r2s_run_dir = os.path.join(
        thisdir, "files_test_r2s", "r2s_examples", "subvoxel_r2s_run"
    )
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test unstructured r2s
    r2s_run_dir = os.path.join(
        thisdir, "files_test_r2s", "r2s_examples", "unstructured_r2s_run"
    )
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test openmc r2s
    try:
        import openmc
    except:
        raise SkipTest
    r2s_run_dir = os.path.join(thisdir, "files_test_r2s", "r2s_examples", "openmc_r2s")
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)


def test_r2s_script():
    # skip test without dagmc
    try:
        from pyne import dagmc
    except ImportError:
        raise SkipTest

    remove_step1_out = False
    # test voxel r2s
    r2s_run_dir = os.path.join(thisdir, "files_test_r2s", "r2s_examples", "r2s_run")
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test sub-voxel r2s
    r2s_run_dir = os.path.join(
        thisdir, "files_test_r2s", "r2s_examples", "subvoxel_r2s_run"
    )
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test unstructured r2s
    r2s_run_dir = os.path.join(
        thisdir, "files_test_r2s", "r2s_examples", "unstructured_r2s_run"
    )
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
    # test openmc r2s
    try:
        import openmc
    except ImportError:
        raise SkipTest
    r2s_run_dir = os.path.join(thisdir, "files_test_r2s", "r2s_examples", "openmc_r2s")
    _r2s_test_step1(r2s_run_dir, remove_step1_out)
    _r2s_test_step2(r2s_run_dir, remove_step1_out)
