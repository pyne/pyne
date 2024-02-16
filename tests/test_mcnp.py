"""PyNE MCNP tools tests"""
import os
import unittest
import struct
import warnings

import pytest

import tables
from numpy.testing import assert_array_equal

from pyne.mesh import Mesh, StatMesh, MeshError, HAVE_PYMOAB
from pyne.material import MultiMaterial, Material
from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
try:
    from pyne import mcnp
    from pyne.mcnp import mats_from_inp
except ImportError:
    pytest.skip()


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

# Test SurfSrc class


@pytest.mark.parametrize("ssrname", ssrnames)
def test_read_header_block(ssrname):
    """Test the read_header() method in the SurfSrc class
    We compare the SurfSrc object variables with expected values from the
    multiple write files
    'mcnp_surfsrc.w', 'mcnpx_surfsrc.w', and 'mcnp6_surfsrc.w'.
    """

    check_read_header_block(ssrname)


def check_read_header_block(ssrname):
    if "mcnp_surfsrc.w" in ssrname:
        ssr = mcnp.SurfSrc(ssrname, "rb")

        try:
            ssr.read_header()
        except:
            pytest.skip()

        # header record values
        assert ssr.kod == "mcnp    "
        assert ssr.ver == "5    "
        assert ssr.loddat == "01232009"
        assert ssr.idtm == " 10/31/11 13:52:39 "
        assert ssr.probid == " 10/31/11 13:52:35 "
        assert (
            ssr.aid ==
            "c Test deck with H20 cube, point n source, "
            "SSW of top surface interactions      ")
        assert ssr.knod == 2
        # table 1 record values
        assert ssr.np1 == 1000
        assert ssr.nrss == 173
        assert ssr.ncrd == 11
        assert ssr.njsw == 1
        assert ssr.niss == 173
        # table 2 record values
        assert ssr.niwr == 0
        assert ssr.mipts == 3
        assert ssr.kjaq == 0

    elif "mcnp6_surfsrc.w" in ssrname:
        ssr = mcnp.SurfSrc(ssrname, "rb")
        try:
            ssr.read_header()
        except:
            pytest.skip()
        # header record values
        assert ssr.kod == "SF_00001"
        assert ssr.ver == "mcnp    6   "
        assert ssr.loddat == " 05/08/13"
        assert ssr.idtm == " 11/18/13 17:50:49 "
        assert ssr.probid == " 11/18/13 17:50:43 "
        assert (
            ssr.aid ==
            "Simple MCNP Example that uses SSW"
            "                       "
            "                        ")
        assert ssr.knod == 2
        # table 1 record values
        assert ssr.np1 == 10000
        assert ssr.nrss == 1710
        assert ssr.ncrd == -11
        # assert_equal(ssrA, ssrB)
        assert ssr.njsw == 1
        assert ssr.niss == 1701
        # table 2 record values
        assert ssr.niwr == 0
        assert ssr.mipts == 37
        assert ssr.kjaq == 0

    elif "mcnpx_surfsrc.w" in ssrname:
        ssr = mcnp.SurfSrc(ssrname, "rb")
        try:
            ssr.read_header()
        except:
            pytest.skip()
        # header record values
        assert ssr.kod == "mcnpx   "
        assert ssr.ver == "2.6.0"
        assert ssr.loddat == "Wed Apr 09 08:00:00 MST 2008"
        assert ssr.idtm == "  10/28/13 02:16:22"
        assert ssr.probid == "  10/28/13 02:16:16"
        assert (
            ssr.aid ==
            "Simple MCNP Example that uses SSW"
            "                                               ")
        assert ssr.knod == 2
        # table 1 record values
        assert ssr.np1 == 10000
        assert ssr.nrss == 1658
        assert ssr.ncrd == 11
        assert ssr.njsw == 1
        assert ssr.niss == 1652
        # table 2 record values
        assert ssr.niwr == 0
        assert ssr.mipts == 35
        assert ssr.kjaq == 0


@pytest.mark.parametrize("ssrname", ssrnames)
def test_compare(ssrname):
    """Test the __cmp__() method in the SurfSrc class
    Tricky to test... this just verifies that comparisons are done right.
    """
    check_compare(ssrname)


def check_compare(ssrname):
    ssrA = mcnp.SurfSrc(ssrname, "rb")
    ssrB = mcnp.SurfSrc(ssrname, "rb")
    try:
        ssrA.read_header()
    except:
        pytest.skip()
    ssrB.read_header()
    assert ssrA == ssrB
    ssrA.close()
    ssrB.close()


@pytest.mark.parametrize("ssrname, sswname", zip(ssrnames,sswnames))
def test_put_header_block(ssrname,sswname):
    """We copy the header block, write to new file, re-read, and compare.
    This tests that information is preserved correctly when written.
    """
    check_put_header_block(ssrname,sswname)


def check_put_header_block(ssrname, sswname):
    ssr = mcnp.SurfSrc(ssrname, "rb")
    ssw = mcnp.SurfSrc(sswname, "wb")
    try:
        ssr.read_header()
    except:
        pytest.skip()
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
    ssw.table1extra = ssr.table1extra
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

    assert ssr.print_header() == sswr.print_header()

    ssr.close()
    sswr.close()

    os.system("rm -f " + sswname)


def test_read_tracklist():
    """We read in tracklists and compare with known values.
    We use a file with a single track for this test.
    """
    ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
    try:
        ssr.read_header()
    except:
        pytest.skip()
    ssr.read_tracklist()

    # print "Length: " + str(len(ssr.tracklist))
    for trackData in ssr.tracklist:
        # Should only be one trackData in tracklist
        # trackData.record is skipped; contains the below components.
        # self.assertEqual(trackData.record  , 0)
        assert trackData.nps == 1
        assert trackData.bitarray == pytest.approx(8.000048e06)
        assert trackData.wgt == pytest.approx(0.99995639)
        assert trackData.erg == pytest.approx(5.54203947)
        assert trackData.tme == pytest.approx(0.17144023)
        assert trackData.x == pytest.approx(-8.05902e-02)
        assert trackData.y == pytest.approx(3.122666098e00)
        assert trackData.z == pytest.approx(5.00000e00)
        assert trackData.u == pytest.approx(-0.35133163)
        assert trackData.v == pytest.approx(0.48465036)
        assert trackData.cs == pytest.approx(0.80104937)
        assert trackData.w == pytest.approx(0.80104937)
    return


def test_read_tracklist_into_different_surface():
    """We read in tracklists and compare with known values.
    We use a file with a single track for this test.
    """
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    try:
        ssr1.read_header()
    except:
        pytest.skip()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    # Update ssr1 with ssr2's tracklist
    ssr1.update_tracklist(ssr2)

    for trackData in ssr1.tracklist:
        assert trackData.nps == 1
        assert trackData.bitarray == pytest.approx(8.000048e06)
        assert trackData.wgt == pytest.approx(0.99995639)
        assert trackData.erg == pytest.approx(5.54203947)
        assert trackData.tme == pytest.approx(0.17144023)
        assert trackData.x == pytest.approx(-8.05902e-02)
        assert trackData.y == pytest.approx(3.122666098e00)
        assert trackData.z == pytest.approx(5.00000e00)
        assert trackData.u == pytest.approx(-0.35133163)
        assert trackData.v == pytest.approx(0.48465036)
        assert trackData.cs == pytest.approx(0.80104937)
        assert trackData.w == pytest.approx(0.80104937)
    return


def test_read_tracklist_into_different_surface_errors():
    """6 Exceptions that are handled by update_tracklist
    We iterate through each type of error and try match each exception
    We try to confirm each error caught by update_tracklist
    """
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    try:
        ssr1.read_header()
    except:
        pytest.skip()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    # TypeError #1: Test with integer '1' in argument
    def wrong_type():
        ssr1.update_tracklist(1)

    pytest.raises(TypeError, wrong_type)

    # AttributeError #2: If there is no header variables in surf_src argument
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")

    def surf_src_arg_no_header():
        ssr1.update_tracklist(ssr2)

    pytest.raises(AttributeError, surf_src_arg_no_header)

    # AttributeError #3: If there are no header variables in surf_src
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()

    def surf_src_no_header():
        ssr1.update_tracklist(ssr2)

    pytest.raises(AttributeError, surf_src_no_header)

    # AttributeError #4: If there is no tracklist in surf_src argument
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()
    ssr1.read_tracklist()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()

    def surf_src_arg_no_tracklist():
        ssr1.update_tracklist(ssr2)

    pytest.raises(AttributeError, surf_src_arg_no_tracklist)

    # AttributeError #5: If there is no tracklist in surf_src
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    ssr1.read_header()

    ssr2 = mcnp.SurfSrc(ssrname_onetrack, "rb")
    ssr2.read_header()
    ssr2.read_tracklist()

    def surf_src_no_tracklist():
        ssr1.update_tracklist(ssr2)

    pytest.raises(AttributeError, surf_src_no_tracklist)

    # ValueError #6: Update ssr1 with ssr1's tracklist
    ssrname = "mcnp5_surfsrc.w"
    ssr1 = mcnp.SurfSrc(ssrname, "rb")
    try:
        ssr1.read_header()
    except:
        pytest.skip()

    ssr1.read_tracklist()

    def update_with_self():
        ssr1.update_tracklist(ssr1)

    pytest.raises(ValueError, update_with_self)

    return


def test_print_header():
    """Check SurfSrc.print_header() against expected resulting string.
    We use a file with a single track for this test, but only use the
    header of this file.
    """
    ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
    try:
        ssr.read_header()
    except:
        pytest.skip()
    # If comparison output needs to be updated, uncomment the below
    #  and do: nosetests test_mcnp.py --nocapture
    # print ssr.print_header()
    assert (
        ssr.print_header() ==
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
    try:
        ssr.read_header()
    except struct.error:
        pytest.skip()
    ssr.read_tracklist()
    # If comparison output needs to be updated, uncomment the below
    #  and do: nosetests test_mcnp.py --nocapture
    try:
        observed = ssr.print_tracklist()
    except struct.error:
        pytest.skip()
    assert (
        observed ==
        "Track Data\n       nps   BITARRAY        WGT        ERG"
        "        TME             X             Y             Z  "
        "        U          V     COSINE  |       W\n         "
        "1 8.00005e+06    0.99996      5.542    0.17144  "
        "-8.05902e-02   3.12267e+00   5.00000e+00   "
        "-0.35133    0.48465    0.80105  |    0.80105 \n")

    return


def _gen_xsdir():
    thisdir = os.path.dirname(__file__)
    xsdir_file = os.path.join(thisdir, "files_test_mcnp", "dummy_xsdir")
    return mcnp.Xsdir(xsdir_file)


def test_xsdir():
    xsdir = _gen_xsdir()

    #  test atomic mass ratio tables
    exp_awr = {
        "1000": "0.99931697",
        "3000": "6.88131188",
        "3003": "3.11111111",
        "3004": "4.11111111",
        "3005": "5.111111111",
        "3009": "9.11111111",
        "0001": "1.000000",
    }
    assert xsdir.awr == exp_awr

    # test xs tables
    assert xsdir.tables[0].name == "1001.44c"
    assert xsdir.tables[0].awr == 1.111111
    assert xsdir.tables[0].filename == "many_xs/1001.555nc"
    assert xsdir.tables[0].access == "0"
    assert xsdir.tables[0].filetype == 1
    assert xsdir.tables[0].address == 4
    assert xsdir.tables[0].tablelength == 55555
    assert xsdir.tables[0].recordlength == 0
    assert xsdir.tables[0].entries == 0
    assert xsdir.tables[0].temperature == 5.5555e05
    assert not xsdir.tables[0].ptable
    assert xsdir.tables[1].filename == "such_data/1001.777nc"
    assert xsdir.tables[1].ptable
    assert xsdir.tables[2].filename == "more_data/1001.999nc"
    assert xsdir.tables[2].ptable


def test_xsdir_find_table():
    xsdir = _gen_xsdir()
    table = xsdir.find_table("1001")
    assert table[0].name == "1001.44c"
    assert table[1].name == "1001.66c"
    assert table[2].name == "1001.70c"


def test_xsdir_to_serpent():
    xsdir = _gen_xsdir()
    output = os.path.join(os.getcwd(), "test_output")
    xsdir.to_xsdata(output)

    with open(output, "r") as f:
        lines = f.readlines()

    exp = [
        (
            "1001.44c 1001.44c 1 1001 0 1.111111 6.44688328094e+15 0"
            " many_xs/1001.555nc\n"
        ),
        (
            "1001.66c 1001.66c 1 1001 0 1.111111 6.44688328094e+15 0"
            " such_data/1001.777nc\n"
        ),
        (
            "1001.70c 1001.70c 1 1001 0 1.111111 6.44688328094e+15 0"
            " more_data/1001.999nc\n"
        ),
    ]

    assert lines == exp
    os.remove(output)


def test_xsdir_nucs():
    xsdir = _gen_xsdir()
    assert xsdir.nucs() == set([10010000])


def test_xsdirtable_to_serpent():
    xsdir = _gen_xsdir()
    line = xsdir.tables[0].to_serpent(".")
    exp_line = (
        "1001.44c 1001.44c 1 1001 0 1.111111 6.44688328094e+15 0"
        " ./many_xs/1001.555nc"
    )
    assert line == exp_line


def test_read_mcnp():
    expected_material = Material(
        nucvec={922350000: 0.04, 922380000: 0.96},
        mass=-1.0,
        density=19.1,
        metadata={
            "comments": (
                " first line of comments second line of "
                "comments third line of comments forth "
                "line of comments"
            ),
            "mat_number": "1",
            "name": " leu",
            "source": " Some http://URL.com",
            "table_ids": {"922350": "15c"},
        },
    )
    expected_material.mass = -1.0  # to avoid reassignment to +1.0
    expected_material_default_lib = Material(
        {
            10000000: 0.037298334378933776,
            60000000: 0.6666767493126631,
            80000000: 0.29602491630840305,
        },
        54.04749412269001,
        1.1,
        6.0,
        {
            "mat_number": "3",
            "HLIB": "42h",
            "NLIB": "60c",
            "PLIB": "01p",
            "table_ids": {},
        },
    )
    expected_multimaterial = MultiMaterial(
        {
            Material(
                {10000000: 0.1118983878322976, 80000000: 0.8881016121677024},
                -1.0,
                0.9,
                3,
                {
                    "comments": (
                        " Here are comments the comments "
                        "continue here are more even more"
                    ),
                    "mat_number": "2",
                    "name": " water",
                    "source": " internet",
                    "table_ids": {"10000": "05c"},
                },
            ): 1,
            Material(
                {10000000: 0.1118983878322976, 80000000: 0.8881016121677024},
                -1.0,
                1.0021552889251644,
                3,
                {
                    "comments": (
                        " Here are comments the comments "
                        "continue here are more even more"
                    ),
                    "mat_number": "2",
                    "name": " water",
                    "source": " internet",
                    "table_ids": {"10000": "05c"},
                },
            ): 1,
            Material(
                {10000000: 0.1118983878322976, 80000000: 0.8881016121677024},
                -1.0,
                1.1,
                3,
                {
                    "comments": (
                        " Here are comments the comments "
                        "continue here are more even more"
                    ),
                    "mat_number": "2",
                    "name": " water",
                    "source": " internet",
                    "table_ids": {"10000": "05c"},
                },
            ): 1,
        }
    )

    read_materials = mats_from_inp("mcnp_inp.txt")
    assert expected_material == pytest.approx(read_materials[1])
    assert expected_material_default_lib == read_materials[3]
    assert (
        list(expected_multimaterial._mats.keys())[0].comp.keys() ==
        list(read_materials[2]._mats.keys())[0].comp.keys())
    for i in range(2):
        assert (
            list(list(expected_multimaterial._mats.keys())[0].comp.values())[i] ==
            pytest.approx(list(list(read_materials[2]._mats.keys())[0].comp.values())[i]))
    assert (
        list(expected_multimaterial._mats.keys())[0].mass ==
        pytest.approx(list(read_materials[2]._mats.keys())[0].mass))
    assert (
        list(expected_multimaterial._mats.keys())[0].density ==
        pytest.approx(list(read_materials[2]._mats.keys())[0].density))
    assert (
        list(expected_multimaterial._mats.keys())[0].atoms_per_molecule ==
        list(read_materials[2]._mats.keys())[0].atoms_per_molecule)
    assert (
        list(expected_multimaterial._mats.keys())[0].metadata ==
        list(read_materials[2]._mats.keys())[0].metadata)
    assert (
        list(expected_multimaterial._mats.keys())[1].comp.keys() ==
        list(read_materials[2]._mats.keys())[1].comp.keys())
    for i in range(2):
        assert (
            list(list(expected_multimaterial._mats.keys())[1].comp.values())[i] ==
            pytest.approx(list(list(read_materials[2]._mats.keys())[1].comp.values())[i]))
    assert (
        list(expected_multimaterial._mats.keys())[1].mass ==
        list(read_materials[2]._mats.keys())[1].mass)
    assert (
        list(expected_multimaterial._mats.keys())[1].density ==
        pytest.approx(list(read_materials[2]._mats.keys())[1].density))
    assert (
        list(expected_multimaterial._mats.keys())[1].atoms_per_molecule ==
        list(read_materials[2]._mats.keys())[1].atoms_per_molecule)
    assert (
        list(expected_multimaterial._mats.keys())[1].metadata ==
        list(read_materials[2]._mats.keys())[1].metadata)
    assert (
        list(expected_multimaterial._mats.keys())[2].density ==
        pytest.approx(list(read_materials[2]._mats.keys())[2].density))


# test to ensure the mats_from_inp function can read repeated mcnp
# materials like
# m1  1001.21c 0.1
#     1002.21c 0.3
#     1001.21c 0.5


def test_read_mcnp_wcomments():
    expected_material = Material(
        nucvec={922350000: 0.04, 922380000: 0.96},
        mass=-1.0,
        density=19.1,
        metadata={
            "comments": (
                " first line of comments second line of "
                "comments third line of comments forth "
                "line of comments"
            ),
            "mat_number": "1",
            "name": " leu",
            "source": " Some http://URL.com",
            "table_ids": {"922350": "15c"},
        },
    )
    expected_material.mass = -1.0  # to avoid reassignment to +1.0

    read_materials = mats_from_inp("mcnp_inp_comments.txt")
    assert expected_material == read_materials[1]


# Test PtracReader class


def test_read_headers():
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
    assert (
        p.problem_title == "Generate a well-defined PTRAC file for PyNE test cases")
    del p

    # 8-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
    assert (
        p.problem_title == "Generate a well-defined PTRAC file for PyNE test cases")
    del p


def test_determine_format():
    # 4-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
    assert p.endianness == "<"
    del p

    # 8-byte ints, little endian
    p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
    assert p.endianness == "<"
    del p


def test_read_events():
    p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")

    evt = {}

    p.read_nps_line()
    assert p.next_event == 1000

    p.read_event_line(evt)
    assert evt["xxx"] == 0.0
    assert evt["yyy"] == 0.0
    assert evt["zzz"] == 0.0
    del p
    del evt


def test_write_to_hdf5():
    test_files = [
        "mcnp_ptrac_i4_little.ptrac",
        "mcnp_ptrac_i8_little.ptrac",
        "mcnp6_ptrac_i4_little.ptrac",
    ]

    for test_file in test_files:
        p = mcnp.PtracReader(test_file)
        h5file = tables.open_file("mcnp_ptrac_hdf5_file.h5", "w")
        tab = h5file.create_table("/", "t", mcnp.PtracEvent, "test")
        p.write_to_hdf5_table(tab)
        tab.flush()
        h5file.close()
        del h5file
        del tab
        del p

        # now check if the data was correctly written.
        # there should be 5 events of type 1000 (src)
        h5file = tables.open_file("mcnp_ptrac_hdf5_file.h5")
        tab = h5file.get_node("/t")
        selected = [1 for x in tab.iterrows() if x["event_type"] == 1000]
        assert len(selected) == 5
        h5file.close()
        del tab
        del h5file

        # clean up
        if os.path.exists("mcnp_ptrac_hdf5_file.h5"):
            os.unlink("mcnp_ptrac_hdf5_file.h5")


# Test Wwinp class. All three function are tested at once because their inputs
# and ouputs are easily strung together.
def test_wwinp_n():
    if not HAVE_PYMOAB:
        pytest.skip()

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, "mcnp_wwinp_wwinp_n.txt")
    expected_h5m = os.path.join(thisdir, "mcnp_wwinp_mesh_n.h5m")
    expected_sm = Mesh(mesh=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), "test_wwinp")

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert ww1.ni == 1
    assert ww1.nr == 10
    assert ww1.ne == [7]
    assert ww1.nf == [15, 8, 6]
    assert ww1.origin == [-100, -100, -100]
    assert ww1.nc == [5, 3, 1]
    assert ww1.nwg == 1
    assert ww1.cm == [[-99, -97, 97, 99, 100], [-50, 60, 100], [100]]
    assert ww1.fm == [[1, 1, 11, 1, 1], [1, 3, 4], [6]]
    assert_array_equal(
        ww1.e, [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]]
    )
    assert (
        ww1.bounds ==
        [
            [
                -100.0,
                -99.0,
                -97.0,
                -79.36363636363636,
                -61.727272727272727,
                -44.090909090909093,
                -26.454545454545453,
                -8.818181818181813,
                8.818181818181813,
                26.454545454545453,
                44.090909090909093,
                61.72727272727272,
                79.363636363636374,
                97.0,
                99.0,
                100.0,
            ],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_n[expected_ve]
        written = ww1.ww_n[written_ve]
        assert_array_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert ww2.ni == 1
    assert ww2.nr == 10
    assert ww2.ne == [7]
    assert ww2.nf == [15, 8, 6]
    assert ww2.origin == [-100, -100, -100]
    assert ww2.nc == [5, 3, 1]
    assert ww2.nwg == 1
    assert ww2.cm == [[-99, -97, 97, 99, 100], [-50, 60, 100], [100]]
    assert ww2.fm == [[1, 1, 11, 1, 1], [1, 3, 4], [6]]
    assert_array_equal(
        ww2.e, [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]]
    )
    assert (
        ww2.bounds ==
        [
            [
                -100.0,
                -99.0,
                -97.0,
                -79.36363636363636,
                -61.727272727272727,
                -44.090909090909093,
                -26.454545454545453,
                -8.818181818181813,
                8.818181818181813,
                26.454545454545453,
                44.090909090909093,
                61.72727272727272,
                79.363636363636374,
                97.0,
                99.0,
                100.0,
            ],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_n[expected_ve]
        written = ww2.ww_n[written_ve]
        assert_array_equal(written, expected)

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
    assert written[0].split()[:-2] == expected[0].split()[:-2]
    assert len(written) == len(expected)
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert float(written[i].split()[j]) == float(expected[i].split()[j])

    os.remove(output)


def test_wwinp_p():
    if not HAVE_PYMOAB:
        pytest.skip()

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, "mcnp_wwinp_wwinp_p.txt")
    expected_h5m = os.path.join(thisdir, "mcnp_wwinp_mesh_p.h5m")
    expected_sm = Mesh(mesh=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), "test_wwinp")

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert ww1.ni == 2
    assert ww1.nr == 10
    assert ww1.ne == [0, 7]
    assert ww1.nf == [1, 8, 6]
    assert ww1.origin == [-100, -100, -100]
    assert ww1.nc == [1, 3, 1]
    assert ww1.nwg == 1
    assert ww1.cm == [[100], [-50, 60, 100], [100]]
    assert ww1.fm == [[1], [1, 3, 4], [6]]
    assert ww1.e[0] == []
    assert_array_equal(
        ww1.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]
    )
    assert (
        ww1.bounds ==
        [
            [-100.0, 100],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_p[expected_ve]
        written = ww1.ww_p[written_ve]
        assert_array_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert ww2.ni == 2
    assert ww2.nr == 10
    assert ww2.ne == [0, 7]
    assert ww2.nf == [1, 8, 6]
    assert ww2.origin == [-100, -100, -100]
    assert ww2.nc == [1, 3, 1]
    assert ww2.nwg == 1
    assert ww2.cm == [[100], [-50, 60, 100], [100]]
    assert ww2.fm == [[1], [1, 3, 4], [6]]
    assert ww2.e[0] == []
    assert_array_equal(
        ww2.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]
    )
    assert (
        ww2.bounds ==
        [
            [-100.0, 100],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_p[expected_ve]
        written = ww2.ww_p[written_ve]
        assert_array_equal(written, expected)

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
    assert written[0].split()[:-2] == expected[0].split()[:-2]
    assert len(written) == len(expected)
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert float(written[i].split()[j]) == float(expected[i].split()[j])

    os.remove(output)


def test_wwinp_np():
    if not HAVE_PYMOAB:
        pytest.skip()

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, "mcnp_wwinp_wwinp_np.txt")
    expected_h5m = os.path.join(thisdir, "mcnp_wwinp_mesh_np.h5m")
    expected_sm = Mesh(mesh=expected_h5m, structured=True)
    output = os.path.join(os.getcwd(), "test_wwinp")

    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert ww1.ni == 2
    assert ww1.nr == 10
    assert ww1.ne == [7, 1]
    assert ww1.nf == [1, 8, 6]
    assert ww1.origin == [-100, -100, -100]
    assert ww1.nc == [1, 3, 1]
    assert ww1.nwg == 1
    assert ww1.cm == [[100], [-50, 60, 100], [100]]
    assert ww1.fm == [[1], [1, 3, 4], [6]]
    assert ww1.e[0] == [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]
    assert_array_equal(ww1.e[1], [100])
    assert (
        ww1.bounds ==
        [
            [-100.0, 100],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_n[expected_ve]
        written = ww1.ww_n[written_ve]
        assert_array_equal(written, expected)

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_p[expected_ve]
        written = ww1.ww_p[written_ve]
        assert_array_equal(written, expected)

    # Create an new object based off of only the mesh attribute of the first
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert ww2.ni == 2
    assert ww2.nr == 10
    assert ww2.ne == [7, 1]
    assert ww2.nf == [1, 8, 6]
    assert ww2.origin == [-100, -100, -100]
    assert ww2.nc == [1, 3, 1]
    assert ww2.nwg == 1
    assert ww2.cm == [[100], [-50, 60, 100], [100]]
    assert ww2.fm == [[1], [1, 3, 4], [6]]
    assert_array_equal(
        ww2.e[0], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]
    )
    assert_array_equal(ww2.e[1], [100])
    assert (
        ww2.bounds ==
        [
            [-100.0, 100],
            [
                -100.0,
                -50.0,
                -13.333333333333336,
                23.333333333333329,
                60.0,
                70.0,
                80.0,
                90.0,
                100.0,
            ],
            [
                -100.0,
                -66.666666666666657,
                -33.333333333333329,
                0.0,
                33.333333333333343,
                66.666666666666657,
                100.0,
            ],
        ])

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_n[expected_ve]
        written = ww2.ww_n[written_ve]
        assert_array_equal(written, expected)

    expected_ves = list(expected_sm.structured_iterate_hex("zyx"))
    written_ves = list(expected_sm.structured_iterate_hex("zyx"))
    for expected_ve, written_ve in zip(expected_ves, written_ves):
        expected = expected_sm.ww_p[expected_ve]
        written = ww2.ww_p[written_ve]
        assert_array_equal(written, expected)

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
    assert written[0].split()[:-2] == expected[0].split()[:-2]
    assert len(written) == len(expected)
    for i in range(1, len(expected)):
        for j in range(0, len(expected[i].split())):
            assert float(written[i].split()[j]) == float(expected[i].split()[j])

    os.remove(output)


# Test Meshtal and Meshtally classes
def test_single_meshtally_meshtal():
    """Test a meshtal file containing a single mesh tally."""

    if not HAVE_PYMOAB:
        pytest.skip()

    thisdir = os.path.dirname(__file__)
    meshtal_file = os.path.join(thisdir, "mcnp_meshtal_single_meshtal.txt")
    expected_h5m = os.path.join(thisdir, "mcnp_meshtal_single_mesh.h5m")
    expected_sm = Mesh(mesh=expected_h5m, structured=True)

    tags = {4: ["n_result", "n_rel_error", "n_total_result", "n_total_rel_error"]}

    meshtal_object = mcnp.Meshtal(meshtal_file, tags, meshes_have_mats=True)
    assert meshtal_object.tally[4].mats != None

    # test Meshtal attributes
    assert meshtal_object.version == "5.mpi"

    assert meshtal_object.ld == "09282010"

    assert meshtal_object.title == "Input file to general test meshtal file"

    assert meshtal_object.histories == 100000

    # test MeshTally attributes
    assert meshtal_object.tally[4].tally_number == 4
    assert meshtal_object.tally[4].particle == "neutron"
    assert meshtal_object.tally[4].dose_response == True
    assert meshtal_object.tally[4].x_bounds == (-200.00, -66.67, 66.67, 200.00)
    assert (
        meshtal_object.tally[4].y_bounds ==
        (-200.00, -120.00, -40.00, 40.00, 120.00, 200.00))
    assert meshtal_object.tally[4].z_bounds == (-200.00, -50.00, 100.00, 200.00)
    assert (
        meshtal_object.tally[4].e_bounds == (0.00e00, 1.00e-01, 2.00e-01, 1.00e00))

    # test vector tags
    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_result[v_e]
        expected = expected_sm.n_result[expected_v_e]
        assert_array_equal(written, expected)
    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_rel_error[v_e]
        expected = expected_sm.n_rel_error[expected_v_e]
        assert_array_equal(written, expected)

    # test total tag
    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_total_result[v_e]
        expected = expected_sm.n_total_result[expected_v_e]
        assert written == expected
    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_total_rel_error[v_e]
        expected = expected_sm.n_total_rel_error[expected_v_e]
        assert written == expected


def test_multiple_meshtally_meshtal():
    """Test a meshtal file containing 4 mesh tallies including neutron and
    photon, single energy group and multiple energy group.
    """

    if not HAVE_PYMOAB:
        pytest.skip()

    thisdir = os.path.dirname(__file__)
    meshtal_file = os.path.join(thisdir, "mcnp_meshtal_multiple_meshtal.txt")

    expected_h5m_4 = os.path.join(thisdir, "mcnp_meshtal_tally_4.h5m")
    expected_sm_4 = Mesh(mesh=expected_h5m_4, structured=True)

    expected_h5m_14 = os.path.join(thisdir, "mcnp_meshtal_tally_14.h5m")
    expected_sm_14 = Mesh(mesh=expected_h5m_14, structured=True)

    expected_h5m_24 = os.path.join(thisdir, "mcnp_meshtal_tally_24.h5m")
    expected_sm_24 = Mesh(mesh=expected_h5m_24, structured=True)

    expected_h5m_34 = os.path.join(thisdir, "mcnp_meshtal_tally_34.h5m")
    expected_sm_34 = Mesh(mesh=expected_h5m_34, structured=True)

    tags = {
        4: ["n_result", "n_rel_error", "n_total_result", "n_total_rel_error"],
        14: ["n_result", "n_rel_error", "n_total_result", "n_total_rel_error"],
        24: ["p_result", "p_rel_error", "p_total_result", "p_total_rel_error"],
        34: ["p_result", "p_rel_error", "p_total_result", "p_total_rel_error"],
    }
    meshtal_object = mcnp.Meshtal(meshtal_file, tags)
    assert meshtal_object.version == "5"
    assert meshtal_object.tally[4].mats == None

    # test meshtally 4
    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm_4.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_result[v_e]
        expected = expected_sm_4.n_result[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm_4.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_rel_error[v_e]
        expected = expected_sm_4.n_rel_error[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm_4.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_total_result[v_e]
        expected = expected_sm_4.n_total_result[expected_v_e]
        assert written == expected

    for v_e, expected_v_e in zip(
        meshtal_object.tally[4].structured_iterate_hex("xyz"),
        expected_sm_4.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[4].n_total_rel_error[v_e]
        expected = expected_sm_4.n_total_rel_error[expected_v_e]
        assert written == expected

    # test meshtally 14
    for v_e, expected_v_e in zip(
        meshtal_object.tally[14].structured_iterate_hex("xyz"),
        expected_sm_14.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[14].n_result[v_e]
        expected = expected_sm_14.n_result[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[14].structured_iterate_hex("xyz"),
        expected_sm_14.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[14].n_rel_error[v_e]
        expected = expected_sm_14.n_rel_error[expected_v_e]
        assert_array_equal(written, expected)

    # test meshtally 24
    for v_e, expected_v_e in zip(
        meshtal_object.tally[24].structured_iterate_hex("xyz"),
        expected_sm_24.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[24].p_result[v_e]
        expected = expected_sm_24.p_result[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[24].structured_iterate_hex("xyz"),
        expected_sm_24.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[24].p_rel_error[v_e]
        expected = expected_sm_24.p_rel_error[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[24].structured_iterate_hex("xyz"),
        expected_sm_24.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[24].p_total_result[v_e]
        expected = expected_sm_24.p_total_result[expected_v_e]
        assert written == expected

    for v_e, expected_v_e in zip(
        meshtal_object.tally[24].structured_iterate_hex("xyz"),
        expected_sm_24.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[24].p_total_rel_error[v_e]
        expected = expected_sm_24.p_total_rel_error[expected_v_e]
        assert written == expected

    # test meshtally 34
    for v_e, expected_v_e in zip(
        meshtal_object.tally[34].structured_iterate_hex("xyz"),
        expected_sm_34.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[34].p_result[v_e]
        expected = expected_sm_34.p_result[expected_v_e]
        assert_array_equal(written, expected)

    for v_e, expected_v_e in zip(
        meshtal_object.tally[34].structured_iterate_hex("xyz"),
        expected_sm_34.structured_iterate_hex("xyz"),
    ):
        written = meshtal_object.tally[34].p_rel_error[v_e]
        expected = expected_sm_34.p_rel_error[expected_v_e]
        assert_array_equal(written, expected)


def test_mesh_to_geom():
    if not HAVE_PYMOAB:
        pytest.skip()

    mats = {
        0: Material({"H1": 1.0, "K39": 1.0}, density=42.0),
        1: Material({"H1": 0.1, "O16": 1.0}, density=43.0),
        2: Material({"He4": 42.0}, density=44.0),
        3: Material({"Tm171": 171.0}, density=45.0),
        4: Material({"C12": 1.0}, density=47.0),
        5: Material({"1002": 1.0}, density=5.0),
    }

    m = Mesh(
        structured_coords=[[0, 1, 2, 3], [0, 1, 2], [0, 1]], mats=mats, structured=True
    )

    geom = mcnp.mesh_to_geom(m)

    exp_geom = (
        "Generated from PyNE Mesh\n"
        "1 1 42.0 1 -2 5 -6 8 -9\n"
        "2 2 43.0 1 -2 6 -7 8 -9\n"
        "3 3 44.0 2 -3 5 -6 8 -9\n"
        "4 4 45.0 2 -3 6 -7 8 -9\n"
        "5 5 47.0 3 -4 5 -6 8 -9\n"
        "6 6 5.0 3 -4 6 -7 8 -9\n"
        "7 0 -1:4:-5:7:-8:9\n"
        "\n"
        "1 px 0.0\n"
        "2 px 1.0\n"
        "3 px 2.0\n"
        "4 px 3.0\n"
        "5 py 0.0\n"
        "6 py 1.0\n"
        "7 py 2.0\n"
        "8 pz 0.0\n"
        "9 pz 1.0\n"
        "\n"
        "C name: 0\n"
        "C density = 42.00000\n"
        "m1\n"
        "     1001 -2.1000e+01\n"
        "     19039 -2.1000e+01\n"
        "C name: 1\n"
        "C density = 43.00000\n"
        "m2\n"
        "     1001 -3.9091e+00\n"
        "     8016 -3.9091e+01\n"
        "C name: 2\n"
        "C density = 44.00000\n"
        "m3\n"
        "     2004 -4.4000e+01\n"
        "C name: 3\n"
        "C density = 45.00000\n"
        "m4\n"
        "     69171 -4.5000e+01\n"
        "C name: 4\n"
        "C density = 47.00000\n"
        "m5\n"
        "     6012 -4.7000e+01\n"
        "C name: 5\n"
        "C density = 5.00000\n"
        "m6\n"
        "     1002 -5.0000e+00\n"
    )

    assert geom == exp_geom
