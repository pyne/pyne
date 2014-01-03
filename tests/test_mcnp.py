"""PyNE MCNP tools tests"""
import os
import unittest
import nose

import nose.tools
from nose.tools import assert_almost_equal
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.tools import assert_raises

import tables

try:
    from pyne import mcnp
    from pyne.mcnp import read_mcnp_inp
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from pyne.material import Material
from pyne.material import MultiMaterial

thisdir = os.path.dirname(__file__)
ssrname1 = os.path.join(thisdir, "mcnp5_surfsrc.w")
sswname1 = os.path.join(thisdir, "copy_mcnp5_surfsrc.w")
ssrname2 = os.path.join(thisdir, "mcnp6_surfsrc.w")
sswname2 = os.path.join(thisdir, "copy_mcnp6_surfsrc.w")
ssrname3 = os.path.join(thisdir, "mcnpx_surfsrc.w")
sswname3 = os.path.join(thisdir, "copy_mcnpx_surfsrc.w")
ssrnames = [ssrname1, ssrname2, ssrname3]
sswnames = [sswname1, sswname2, sswname3]
ssrname_onetrack = os.path.join(thisdir, "mcnp_surfsrc_onetrack.w")

# mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest
    HAVE_PYTAPS = False
    pass

from pyne.mesh import Mesh, StatMesh, MeshError
from numpy.testing import assert_array_equal


# Test SurfSrc class
def test_read_header_block():
    """Test the read_header() method in the SurfSrc class
    We compare the SurfSrc object variables with expected values from the
    multiple write files
    'mcnp_surfsrc.w', 'mcnpx_surfsrc.w', and 'mcnp6_surfsrc.w'.
    """

    for ssrname in ssrnames:
        yield check_read_header_block, ssrname


def check_read_header_block(ssrname):
        if 'mcnp_surfsrc.w' in ssrname:
            ssr = mcnp.SurfSrc(ssrname, 'rb')
            ssr.read_header()

            # header record values
            assert_equal(ssr.kod, "mcnp    ")
            assert_equal(ssr.ver, "5    ")
            assert_equal(ssr.loddat, "01232009")
            assert_equal(ssr.idtm, " 10/31/11 13:52:39 ")
            assert_equal(ssr.probid, " 10/31/11 13:52:35 ")
            assert_equal(ssr.aid,
                         "c Test deck with H20 cube, point n source, "
                         "SSW of top surface interactions      ")
            assert_equal(ssr.knod, 2)
            # table 1 record values
            assert_equal(ssr.np1, 1000)
            assert_equal(ssr.nrss, 173)
            assert_equal(ssr.ncrd, 11)
            assert_equal(ssr.njsw, 1)
            assert_equal(ssr.niss, 173)
            # table 2 record values
            assert_equal(ssr.niwr, 0)
            assert_equal(ssr.mipts, 3)
            assert_equal(ssr.kjaq, 0)

        elif 'mcnp6_surfsrc.w' in ssrname:
            ssr = mcnp.SurfSrc(ssrname, 'rb')
            ssr.read_header()

            # header record values
            assert_equal(ssr.kod, "SF_00001")
            assert_equal(ssr.ver, "mcnp    6   ")
            assert_equal(ssr.loddat, " 05/08/13")
            assert_equal(ssr.idtm, " 11/18/13 17:50:49 ")
            assert_equal(ssr.probid, " 11/18/13 17:50:43 ")
            assert_equal(ssr.aid, "Simple MCNP Example that uses SSW"
                         "                       "
                         "                        ")
            assert_equal(ssr.knod, 2)
            # table 1 record values
            assert_equal(ssr.np1, 10000)
            assert_equal(ssr.nrss, 1710)
            assert_equal(ssr.ncrd, -11)
        #assert_equal(ssrA, ssrB)
            assert_equal(ssr.njsw, 1)
            assert_equal(ssr.niss, 1701)
            # table 2 record values
            assert_equal(ssr.niwr, 0)
            assert_equal(ssr.mipts, 37)
            assert_equal(ssr.kjaq, 0)

        elif 'mcnpx_surfsrc.w' in ssrname:
            ssr = mcnp.SurfSrc(ssrname, 'rb')
            ssr.read_header()

            # header record values
            assert_equal(ssr.kod, "mcnpx   ")
            assert_equal(ssr.ver, "2.6.0")
            assert_equal(ssr.loddat, "Wed Apr 09 08:00:00 MST 2008")
            assert_equal(ssr.idtm, "  10/28/13 02:16:22")
            assert_equal(ssr.probid, "  10/28/13 02:16:16")
            assert_equal(ssr.aid, "Simple MCNP Example that uses SSW"
                         "                                               ")
            assert_equal(ssr.knod, 2)
            # table 1 record values
            assert_equal(ssr.np1, 10000)
            assert_equal(ssr.nrss, 1658)
            assert_equal(ssr.ncrd, 11)
            assert_equal(ssr.njsw, 1)
            assert_equal(ssr.niss, 1652)
            # table 2 record values
            assert_equal(ssr.niwr, 0)
            assert_equal(ssr.mipts, 35)
            assert_equal(ssr.kjaq, 0)


def test_compare():
    """Test the __cmp__() method in the SurfSrc class
    Tricky to test... this just verifies that comparisons are done right.
    """
    for ssrname in ssrnames:
        yield check_compare, ssrname


def check_compare(ssrname):
        ssrA = mcnp.SurfSrc(ssrname, 'rb')
        ssrB = mcnp.SurfSrc(ssrname, 'rb')
        ssrA.read_header()
        ssrB.read_header()
        assert_true(ssrA == ssrB)
        ssrA.close()
        ssrB.close()


def test_put_header_block():
    """We copy the header block, write to new file, re-read, and compare.
    This tests that information is preserved correctly when written.
    """
    for ssrname, sswname in zip(ssrnames, sswnames):
        yield check_put_header_block, ssrname, sswname


def check_put_header_block(ssrname, sswname):
        ssr = mcnp.SurfSrc(ssrname, "rb")
        ssw = mcnp.SurfSrc(sswname, "wb")
        ssr.read_header()

        # header record values
        ssw.kod = ssr.kod
        ssw.ver = ssr.ver
        ssw.loddat = ssr.loddat
        ssw.idtm = ssr.idtm
        ssw.probid = ssr.probid
        ssw.aid = ssr.aid
        ssw.knod = ssr.knod
        # table 1 record values
        ssw.np1 = ssr.orignp1  # ssr.np1
        ssw.nrss = ssr.nrss
        ssw.ncrd = ssr.ncrd
        ssw.njsw = ssr.njsw
        ssw.niss = ssr.niss
        # table 2 record values
        ssw.niwr = ssr.niwr
        ssw.mipts = ssr.mipts
        ssw.kjaq = ssr.kjaq
        ssw.table2extra = ssr.table2extra
        # surface info record list
        ssw.surflist = ssr.surflist
        # summary table record values
        ssw.summary_table = ssr.summary_table
        ssw.summary_extra = ssr.summary_extra

        ssw.write_header()
        ssw.close()

        sswr = mcnp.SurfSrc(sswname, "rb")
        sswr.read_header()

        assert_equal(ssr.print_header(), sswr.print_header())

        ssr.close()
        sswr.close()

        os.system("rm -f " + sswname)


def test_read_tracklist():
    """We read in tracklists and compare with known values.
    We use a file with a single track for this test.
    """
    ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr.read_header()
    ssr.read_tracklist()

    # print "Length: " + str(len(ssr.tracklist))
    for trackData in ssr.tracklist:
        # Should only be one trackData in tracklist
        # trackData.record is skipped; contains the below components.
        # self.assertEqual(trackData.record  , 0)
        assert_equal(trackData.nps, 1)
        assert_almost_equal(trackData.bitarray, 8.000048e+06)
        assert_almost_equal(trackData.wgt, 0.99995639)
        assert_almost_equal(trackData.erg, 5.54203947)
        assert_almost_equal(trackData.tme, 0.17144023)
        assert_almost_equal(trackData.x, -8.05902e-02)
        assert_almost_equal(trackData.y, 3.122666098e+00)
        assert_almost_equal(trackData.z, 5.00000e+00)
        assert_almost_equal(trackData.u, -0.35133163)
        assert_almost_equal(trackData.v, 0.48465036)
        assert_almost_equal(trackData.cs, 0.80104937)
        assert_almost_equal(trackData.w, 0.80104937)
    return


def test_read_tracklist_into_different_surface():
    """We read in tracklists and compare with known values.
    We use a file with a single track for this test.
    """
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    # Update ssr1 with ssr2's tracklist
    ssr1.update_tracklist(ssr2)

    for trackData in ssr1.tracklist:
        assert_equal(trackData.nps, 1)
        assert_almost_equal(trackData.bitarray, 8.000048e+06)
        assert_almost_equal(trackData.wgt, 0.99995639)
        assert_almost_equal(trackData.erg, 5.54203947)
        assert_almost_equal(trackData.tme, 0.17144023)
        assert_almost_equal(trackData.x, -8.05902e-02)
        assert_almost_equal(trackData.y, 3.122666098e+00)
        assert_almost_equal(trackData.z, 5.00000e+00)
        assert_almost_equal(trackData.u, -0.35133163)
        assert_almost_equal(trackData.v, 0.48465036)
        assert_almost_equal(trackData.cs, 0.80104937)
        assert_almost_equal(trackData.w, 0.80104937)
    return


def test_read_tracklist_into_different_surface_errors():
    """ 6 Exceptions that are handled by update_tracklist
    We iterate through each type of error and try match each exception
    We try to confirm each error caught by update_tracklist
    """
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    # TypeError #1: Test with integer '1' in argument
    def wrong_type():
        ssr1.update_tracklist(1)

    assert_raises(TypeError, wrong_type)

    # AttributeError #2: If there is no header variables in surf_src argument
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")

    def surf_src_arg_no_header():
        ssr1.update_tracklist(ssr2)

    assert_raises(AttributeError, surf_src_arg_no_header)

    # AttributeError #3: If there are no header variables in surf_src
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()

    def surf_src_no_header():
        ssr1.update_tracklist(ssr2)

    assert_raises(AttributeError, surf_src_no_header)

    # AttributeError #4: If there is no tracklist in surf_src argument
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()

    def surf_src_arg_no_tracklist():
        ssr1.update_tracklist(ssr2)

    assert_raises(AttributeError, surf_src_arg_no_tracklist)

    # AttributeError #5: If there is no tracklist in surf_src
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    def surf_src_no_tracklist():
        ssr1.update_tracklist(ssr2)

    assert_raises(AttributeError, surf_src_no_tracklist)

    # ValueError #6: Update ssr1 with ssr1's tracklist
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()
    ssr1.read_tracklist()

    def update_with_self():
        ssr1.update_tracklist(ssr1)

    assert_raises(ValueError, update_with_self)

    return


def test_print_header():
    """Check SurfSrc.print_header() against expected resulting string.
    We use a file with a single track for this test, but only use the
    header of this file.
    """
    ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr.read_header()
    # If comparison output needs to be updated, uncomment the below
    #  and do: nosetests test_mcnp.py --nocapture
    #print ssr.print_header()
    assert_equal(ssr.print_header(),
                 "Code: mcnp     (version: 5    ) [01232009]\n"
                 "Problem info: ( 07/05/12 17:50:19 )"
                 "  07/05/12 17:50:16 \n"
                 "c Test deck with H20 cube, point n source,"
                 " SSW of top surface interactions      \n"
                 "Showing dump #2\n"
                 "1 histories, 1 tracks, 11 record size, "
                 "1 surfaces, 1 histories\n"
                 "0 cells, source particle: 3, macrobody facet flag: 0\n"
                 "Surface [6]: facet -1, type [4]"
                 " with 1 parameters: ( [5.0])\n"
                 "Summary Table: [0, 0, 1, 1, 1, 1,"
                 " 0, 0, 0, 0, 0, 0, 0, 0, 0]")

    return


def test_print_tracklist():
    """Check SurfSrc.print_tracklist() against expected resulting string.
    We use a file with a single track for this test.
    """
    ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr.read_header()
    ssr.read_tracklist()
    # If comparison output needs to be updated, uncomment the below
    #  and do: nosetests test_mcnp.py --nocapture
    #print ssr.print_tracklist()
    assert_equal(ssr.print_tracklist(),
                 'Track Data\n       nps   BITARRAY        WGT        ERG'
                 '        TME             X             Y             Z  '
                 '        U          V     COSINE  |       W\n         '
                 '1 8.00005e+06    0.99996      5.542    0.17144  '
                 '-8.05902e-02   3.12267e+00   5.00000e+00   '
                 '-0.35133    0.48465    0.80105  |    0.80105 \n')

    return


def test_read_mcnp():

    expected_material = Material(nucvec={922350000: 0.04, 922380000: 0.96},
                                 mass=-1.0,
                                 density=19.1,
                                 attrs={"comments": (
                                     " first line of comments second line of "
                                     "comments third line of comments forth "
                                     "line of comments"),
                                     "mat_number": "1",
                                     "name": " leu",
                                     "source": " Some http://URL.com",
                                     "table_ids": {'922350': "15c"}})
    expected_material.mass = -1.0  # to avoid reassignment to +1.0

    expected_multimaterial = MultiMaterial({
        Material(
            {10000000: 0.11189838783149784, 80000000: 0.8881016121685023},
            -1.0, 0.9, 3,
            {"comments":
                (" Here are comments the comments "
                    "continue here are more even more"),
                "mat_number": "2",
                "name": " water",
                "source": " internet",
                "table_ids": {'10000': "05c"}}): 1,
        Material(
            {10000000: 0.11189838783149784, 80000000: 0.8881016121685023},
            -1.0,
            1.0021552889223864, 3,
            {"comments":
                (" Here are comments the comments "
                    "continue here are more even more"),
                "mat_number": "2",
                "name": " water",
                "source": " internet",
                "table_ids": {'10000': "05c"}}): 1})

    read_materials = read_mcnp_inp('mcnp_inp.txt')
    assert_equal(expected_material, read_materials[0])
    assert_equal(
        expected_multimaterial._mats.keys()[0].comp,
        read_materials[1]._mats.keys()[0].comp)
    assert_equal(
        expected_multimaterial._mats.keys()[0].mass,
        read_materials[1]._mats.keys()[0].mass)
    assert_equal(
        expected_multimaterial._mats.keys()[0].density,
        read_materials[1]._mats.keys()[0].density)
    assert_equal(
        expected_multimaterial._mats.keys()[0].atoms_per_mol,
        read_materials[1]._mats.keys()[0].atoms_per_mol)
    assert_equal(
        expected_multimaterial._mats.keys()[0].attrs,
        read_materials[1]._mats.keys()[0].attrs)
    assert_equal(
        expected_multimaterial._mats.keys()[1].comp,
        read_materials[1]._mats.keys()[1].comp)
    assert_equal(
        expected_multimaterial._mats.keys()[1].mass,
        read_materials[1]._mats.keys()[1].mass)
    assert_equal(
        expected_multimaterial._mats.keys()[1].density,
        read_materials[1]._mats.keys()[1].density)
    assert_equal(
        expected_multimaterial._mats.keys()[1].atoms_per_mol,
        read_materials[1]._mats.keys()[1].atoms_per_mol)
    assert_equal(
        expected_multimaterial._mats.keys()[1].attrs,
        read_materials[1]._mats.keys()[1].attrs)


# Test PtracReader class
def test_read_headers():
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
    assert_equal(
        p.problem_title,
        "Generate a well-defined PTRAC file for PyNE test cases")
    del p

    # 8-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
    assert_equal(
        p.problem_title,
        "Generate a well-defined PTRAC file for PyNE test cases")
    del p


def test_determine_format():
    # 4-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
    assert_equal(p.endianness, "<")
    del p

    # 8-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
    assert_equal(p.endianness, "<")
    del p


def test_read_events():
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")

    evt = {}

    p.read_nps_line()
    assert_equal(p.next_event, 1000)

    p.read_event_line(evt)
    assert_equal(evt["xxx"], 0.0)
    assert_equal(evt["yyy"], 0.0)
    assert_equal(evt["zzz"], 0.0)
    del p
    del evt


def test_write_to_hdf5():
    test_files = ["mcnp_ptrac_i4_little.ptrac",
                  "mcnp_ptrac_i8_little.ptrac"]

    for test_file in test_files:
        p = mcnp.PtracReader(test_file)
        h5file = tables.openFile("mcnp_ptrac_hdf5_file.h5", "w")
        tab = h5file.createTable("/", "t", mcnp.PtracEvent, "test")
        p.write_to_hdf5_table(tab)
        tab.flush()
        h5file.close()
        del h5file
        del tab
        del p

        # now check if the data was correctly written.
        # there should be 5 events of type 1000 (src)
        h5file = tables.openFile("mcnp_ptrac_hdf5_file.h5")
        tab = h5file.getNode("/t")
        selected = [1 for x in tab.iterrows() if x["event_type"] == 1000]
        assert_equal(len(selected), 5)
        h5file.close()
        del tab
        del h5file

        # clean up
        if os.path.exists("mcnp_ptrac_hdf5_file.h5"):
            os.unlink("mcnp_ptrac_hdf5_file.h5")


# Test Wwinp class. All three function are tested at once because their inputs
# and ouputs are easily strung together.
def test_wwinp_n():
    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_n.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_n.h5m')
    expected_sm = Mesh(mesh_file=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), 'test_wwinp')

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert_equal(ww1.ni, 1)
    assert_equal(ww1.nr, 10)
    assert_equal(ww1.ne,  [7])
    assert_equal(ww1.nf, [15, 8, 6])
    assert_equal(ww1.origin, [-100, -100, -100])
    assert_equal(ww1.nc, [5, 3, 1])
    assert_equal(ww1.nwg, 1)
    assert_equal(ww1.cm, [[-99, -97,  97, 99, 100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1, 1, 11, 1, 1], [1, 3, 4], [6]])
    assert_array_equal(
        ww1.e,
        [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]])
    assert_equal(
        ww1.bounds,
        [[-100.0, -99.0, -97.0, -79.36363636363636,
            -61.727272727272727, -44.090909090909093,
            -26.454545454545453, -8.818181818181813,
            8.818181818181813, 26.454545454545453,
            44.090909090909093, 61.72727272727272,
            79.363636363636374, 97.0, 99.0, 100.0],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 15):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww1.structured_get_hex(x, y, z)
                    written = ww1.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert_equal(ww2.ni, 1)
    assert_equal(ww2.nr, 10)
    assert_equal(ww2.ne,  [7])
    assert_equal(ww2.nf, [15, 8, 6])
    assert_equal(ww2.origin, [-100, -100, -100])
    assert_equal(ww2.nc, [5, 3, 1])
    assert_equal(ww2.nwg, 1)
    assert_equal(ww2.cm, [[-99, -97,  97, 99, 100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1, 1, 11, 1, 1], [1, 3, 4], [6]])
    assert_array_equal(
        ww2.e,
        [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]])
    assert_equal(
        ww2.bounds,
        [[-100.0, -99.0, -97.0, -79.36363636363636,
            -61.727272727272727, -44.090909090909093,
            -26.454545454545453, -8.818181818181813,
            8.818181818181813, 26.454545454545453,
            44.090909090909093, 61.72727272727272,
            79.363636363636374, 97.0, 99.0, 100.0],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 15):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww2.structured_get_hex(x, y, z)
                    written = ww2.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # write a new wwinp file and verify that is same wwinp file used as an
    # input to this test
    ww2.write_wwinp(output)
    expected_output = wwinp_file

    with open(output) as f:
        written = f.readlines()

    with open(expected_output) as f:
        expected = f.readlines()

    # check to make sure file are the same except for the data/time info
    # on line 1
    assert_equal(written[0].split()[:-2], expected[0].split()[:-2])
    assert_equal(len(written), len(expected))
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert_equal(
                float(written[i].split()[j]),
                float(expected[i].split()[j]))

    os.remove(output)


def test_wwinp_p():
    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_p.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_p.h5m')
    expected_sm = Mesh(mesh_file=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), 'test_wwinp')

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert_equal(ww1.ni, 2)
    assert_equal(ww1.nr, 10)
    assert_equal(ww1.ne,  [0, 7])
    assert_equal(ww1.nf, [1, 8, 6])
    assert_equal(ww1.origin, [-100, -100, -100])
    assert_equal(ww1.nc, [1, 3, 1])
    assert_equal(ww1.nwg, 1)
    assert_equal(ww1.cm, [[100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1], [1, 3, 4], [6]])
    assert_equal(ww1.e[0], [])
    assert_array_equal(
        ww1.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_equal(
        ww1.bounds,
        [[-100.0, 100],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww1.structured_get_hex(x, y, z)
                    written = ww1.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert_equal(ww2.ni, 2)
    assert_equal(ww2.nr, 10)
    assert_equal(ww2.ne,  [0, 7])
    assert_equal(ww2.nf, [1, 8, 6])
    assert_equal(ww2.origin, [-100, -100, -100])
    assert_equal(ww2.nc, [1, 3, 1])
    assert_equal(ww2.nwg, 1)
    assert_equal(ww2.cm, [[100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1], [1, 3, 4], [6]])
    assert_equal(ww2.e[0], [])
    assert_array_equal(
        ww2.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_equal(
        ww2.bounds,
        [[-100.0, 100],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww2.structured_get_hex(x, y, z)
                    written = ww2.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # write a new wwinp file and verify that is same wwinp file used as an
    # input to this test
    ww2.write_wwinp(output)
    expected_output = wwinp_file

    with open(output) as f:
        written = f.readlines()

    with open(expected_output) as f:
        expected = f.readlines()

    # check to make sure file are the same except for the data/time info
    # on line 1
    assert_equal(written[0].split()[:-2], expected[0].split()[:-2])
    assert_equal(len(written), len(expected))
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert_equal(
                float(written[i].split()[j]), float(expected[i].split()[j]))

    os.remove(output)


def test_wwinp_np():
    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_np.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_np.h5m')
    expected_sm = Mesh(mesh_file=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), 'test_wwinp')

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert_equal(ww1.ni, 2)
    assert_equal(ww1.nr, 10)
    assert_equal(ww1.ne,  [7, 1])
    assert_equal(ww1.nf, [1, 8, 6])
    assert_equal(ww1.origin, [-100, -100, -100])
    assert_equal(ww1.nc, [1, 3, 1])
    assert_equal(ww1.nwg, 1)
    assert_equal(ww1.cm, [[100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1], [1, 3, 4], [6]])
    assert_equal(
        ww1.e[0], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_array_equal(ww1.e[1], [100])
    assert_equal(
        ww1.bounds,
        [[-100.0, 100],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww1.structured_get_hex(x, y, z)
                    written = ww1.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 2):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww1.structured_get_hex(x, y, z)
                    written = ww1.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert_equal(ww2.ni, 2)
    assert_equal(ww2.nr, 10)
    assert_equal(ww2.ne,  [7, 1])
    assert_equal(ww2.nf, [1, 8, 6])
    assert_equal(ww2.origin, [-100, -100, -100])
    assert_equal(ww2.nc, [1, 3, 1])
    assert_equal(ww2.nwg, 1)
    assert_equal(ww2.cm, [[100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1], [1, 3, 4], [6]])
    assert_array_equal(
        ww2.e[0], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_array_equal(ww2.e[1], [100])
    assert_equal(
        ww2.bounds,
        [[-100.0, 100],
         [-100.0, -50.0, -13.333333333333336,
          23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0],
         [-100.0, -66.666666666666657, -33.333333333333329,
          0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww2.structured_get_hex(x, y, z)
                    written = ww2.mesh.getTagHandle(
                        'ww_n_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    for x in range(0, 1):
        for y in range(0, 8):
            for z in range(0, 6):
                for e_group in range(1, 2):
                    expected_voxel = expected_sm.structured_get_hex(x, y, z)
                    expected = expected_sm.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[expected_voxel]
                    written_voxel = ww2.structured_get_hex(x, y, z)
                    written = ww2.mesh.getTagHandle(
                        'ww_p_group_00{0}'.format(e_group))[written_voxel]
                    assert_equal(written, expected)

    # write a new wwinp file and verify that is same wwinp file used as an
    # input to this test
    ww2.write_wwinp(output)
    expected_output = wwinp_file

    with open(output) as f:
        written = f.readlines()

    with open(expected_output) as f:
        expected = f.readlines()

    # check to make sure file are the same except for the data/time info
    # on line 1
    assert_equal(written[0].split()[:-2], expected[0].split()[:-2])
    assert_equal(len(written), len(expected))
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert_equal(
                float(written[i].split()[j]), float(expected[i].split()[j]))

    os.remove(output)


# Test Meshtal and Meshtally classes
def test_single_meshtally_meshtal():
    """Test a meshtal file containing a single mesh tally.
    """

    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    meshtal_file = os.path.join(thisdir, "mcnp_meshtal_single_meshtal.txt")
    expected_h5m = os.path.join(thisdir, "tally_single.h5m")
    expected_sm = Mesh(mesh_file=expected_h5m, structured=True)

    meshtal_object = mcnp.Meshtal(meshtal_file)

    # test Meshtal attributes
    assert_equal(meshtal_object.version, 5)

    assert_equal(meshtal_object.ld, "09282010")

    assert_equal(meshtal_object.title,
                 "Input file to general test meshtal file")

    assert_equal(meshtal_object.histories, 100000)

    # test MeshTally attributes
    assert_equal(meshtal_object.tally[4].tally_number, 4)
    assert_equal(meshtal_object.tally[4].particle, "n")
    assert_equal(meshtal_object.tally[4].dose_response, True)
    assert_equal(
        meshtal_object.tally[4].x_bounds,
        [-200.00, -66.67, 66.67, 200.00])
    assert_equal(
        meshtal_object.tally[4].y_bounds,
        [-200.00, -120.00, -40.00, 40.00, 120.00, 200.00])
    assert_equal(
        meshtal_object.tally[4].z_bounds,
        [-200.00, -50.00, 100.00, 200.00])
    assert_equal(
        meshtal_object.tally[4].e_bounds,
        [0.00E+00, 1.00E-01, 2.00E-01, 1.00E+00])

    # test mesh attributes
    for e_group in range(1, 4):
        for v_e, expected_v_e in zip(
                meshtal_object.tally[4].structured_iterate_hex("xyz"),
                expected_sm.structured_iterate_hex("xyz")):
            written = meshtal_object.tally[4].mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[v_e]
            expected = expected_sm.mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[expected_v_e]
            assert_equal(written, expected)


def test_multiple_meshtally_meshtal():
    """Test a meshtal file containing 4 mesh tallies including neutron and
    photon, single energy group and multiple energy group.
    """

    if not HAVE_PYTAPS:
        raise SkipTest

    thisdir = os.path.dirname(__file__)
    meshtal_file = os.path.join(thisdir, "mcnp_meshtal_multiple_meshtal.txt")

    expected_h5m_4 = os.path.join(thisdir, "mcnp_meshtal_tally_4.h5m")
    expected_sm_4 = Mesh(mesh_file=expected_h5m_4, structured=True)

    expected_h5m_14 = os.path.join(thisdir, "mcnp_meshtal_tally_14.h5m")
    expected_sm_14 = Mesh(mesh_file=expected_h5m_14, structured=True)

    expected_h5m_24 = os.path.join(thisdir, "mcnp_meshtal_tally_24.h5m")
    expected_sm_24 = Mesh(mesh_file=expected_h5m_24, structured=True)

    expected_h5m_34 = os.path.join(thisdir, "mcnp_meshtal_tally_34.h5m")
    expected_sm_34 = Mesh(mesh_file=expected_h5m_34, structured=True)

    meshtal_object = mcnp.Meshtal(meshtal_file)

    for e_group in range(1, 7):
        for v_e, expected_v_e in zip(
                meshtal_object.tally[4].structured_iterate_hex("xyz"),
                expected_sm_4.structured_iterate_hex("xyz")):
            written = meshtal_object.tally[4].mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[v_e]
            expected = expected_sm_4.mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[expected_v_e]
            assert_equal(written, expected)

    for e_group in range(1, 2):
        for v_e, expected_v_e in zip(
                meshtal_object.tally[14].structured_iterate_hex("xyz"),
                expected_sm_14.structured_iterate_hex("xyz")):
            written = meshtal_object.tally[14].mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[v_e]
            expected = expected_sm_14.mesh.getTagHandle(
                "n_group_00{0}".format(e_group))[expected_v_e]
            assert_equal(written, expected)

    for e_group in range(1, 7):
        for v_e, expected_v_e in zip(
                meshtal_object.tally[24].structured_iterate_hex("xyz"),
                expected_sm_24.structured_iterate_hex("xyz")):
            written = meshtal_object.tally[24].mesh.getTagHandle(
                "p_group_00{0}".format(e_group))[v_e]
            expected = expected_sm_24.mesh.getTagHandle(
                "p_group_00{0}".format(e_group))[expected_v_e]
            assert_equal(written, expected)

    for e_group in range(1, 2):
        for v_e, expected_v_e in zip(
                meshtal_object.tally[34].structured_iterate_hex("xyz"),
                expected_sm_34.structured_iterate_hex("xyz")):
            written = meshtal_object.tally[34].mesh.getTagHandle(
                "p_group_00{0}".format(e_group))[v_e]
            expected = expected_sm_34.mesh.getTagHandle(
                "p_group_00{0}".format(e_group))[expected_v_e]
            assert_equal(written, expected)
