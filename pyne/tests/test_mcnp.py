"""PyNE MCNP tools tests"""
import os
import unittest
import nose

import nose.tools
from nose.tools import assert_equal

import tables

from pyne import mcnp
from pyne.mcnp import read_mcnp_inp
from pyne.material import Material
from pyne.material import MultiMaterial

thisdir = os.path.dirname(__file__)
ssrname = os.path.join(thisdir,"mcnp_surfsrc.w")
sswname = os.path.join(thisdir,"copy_mcnp_surfsrc.w")
ssrname_onetrack = os.path.join(thisdir,"mcnp_surfsrc_onetrack.w")

from itaps import iMesh
from r2s.scdmesh import ScdMesh

from numpy.testing import assert_array_equal

# Test methods for the SurfSrc class
class TestSurfSrc(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_header_block(self):
        """Test the read_header() method in the SurfSrc class
        We compare the SurfSrc object variables with expected values from the
        file 'mcnp_surfsrc.w'.
        """
        ssr = mcnp.SurfSrc(ssrname, 'rb')
        ssr.read_header()

        # header record values
        self.assertEqual(ssr.kod   , "mcnp    ")
        self.assertEqual(ssr.ver   , "5    ")
        self.assertEqual(ssr.loddat, "01232009")
        self.assertEqual(ssr.idtm  , " 10/31/11 13:52:39 ")
        self.assertEqual(ssr.probid, " 10/31/11 13:52:35 ")
        self.assertEqual(ssr.aid   , "c Test deck with H20 cube, point n source, SSW of top surface interactions      ")
        self.assertEqual(ssr.knod  , 2)
        # table 1 record values
        self.assertEqual(ssr.np1   , 1000)
        self.assertEqual(ssr.nrss  , 173)
        self.assertEqual(ssr.ncrd  , 11)
        self.assertEqual(ssr.njsw  , 1)
        self.assertEqual(ssr.niss  , 173)
        # table 2 record values
        self.assertEqual(ssr.niwr  , 0)
        self.assertEqual(ssr.mipts , 3)
        self.assertEqual(ssr.kjaq  , 0)
        
    def test_compare(self):
        """Test the compare() method in the SurfSrc class
        Tricky to test... this just verifies that comparisons are done right.
        """
        ssrA = mcnp.SurfSrc(ssrname, 'rb')
        ssrB = mcnp.SurfSrc(ssrname, 'rb')
        ssrA.read_header()
        ssrB.read_header()
        self.assertTrue(ssrA.compare(ssrB))
        ssrA.close()
        ssrB.close()

    def test_put_header_block(self):
        """We copy the header block, write to new file, re-read, and compare.
        This tests that information is preserved correctly when written.
        """
        ssr = mcnp.SurfSrc(ssrname, "rb")
        ssw = mcnp.SurfSrc(sswname, "wb")
        ssr.read_header()

        # header record values
        ssw.kod    = ssr.kod
        ssw.ver    = ssr.ver
        ssw.loddat = ssr.loddat
        ssw.idtm   = ssr.idtm
        ssw.probid = ssr.probid
        ssw.aid    = ssr.aid
        ssw.knod   = ssr.knod
        # table 1 record values
        ssw.np1    = ssr.orignp1  # ssr.np1
        ssw.nrss   = ssr.nrss  
        ssw.ncrd   = ssr.ncrd
        ssw.njsw   = ssr.njsw
        ssw.niss   = ssr.niss
        # table 2 record values
        ssw.niwr   = ssr.niwr
        ssw.mipts  = ssr.mipts
        ssw.kjaq   = ssr.kjaq
        ssw.table2extra = ssr.table2extra
        # surface info record list
        ssw.surflist     = ssr.surflist
        # summary table record values
        ssw.summaryTable = ssr.summaryTable
        ssw.summaryExtra = ssr.summaryExtra

        ssw.put_header()
        ssw.put_table_1()
        ssw.put_table_2()
        ssw.put_surface_info()
        ssw.put_summary()
        ssw.close()

        sswr = mcnp.SurfSrc(sswname, "rb")
        sswr.read_header()
        
        self.assertEqual(ssr.print_header(), sswr.print_header())
        
        ssr.close()
        sswr.close()
        
        os.system("rm -f " + sswname)

        return

    def test_read_tracklist(self):
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
            self.assertEqual(trackData.nps     , 1) 
            self.assertAlmostEqual(trackData.bitarray, 8.000048e+06) 
            self.assertAlmostEqual(trackData.wgt     , 0.99995639) 
            self.assertAlmostEqual(trackData.erg     , 5.54203947) 
            self.assertAlmostEqual(trackData.tme     , 0.17144023) 
            self.assertAlmostEqual(trackData.x       , -8.05902e-02) 
            self.assertAlmostEqual(trackData.y       , 3.122666098e+00) 
            self.assertAlmostEqual(trackData.z       , 5.00000e+00) 
            self.assertAlmostEqual(trackData.u       , -0.35133163) 
            self.assertAlmostEqual(trackData.v       , 0.48465036) 
            self.assertAlmostEqual(trackData.cs      , 0.80104937) 
            self.assertAlmostEqual(trackData.w       , 0.80104937) 
        return

    def test_print_header(self):
        """Check SurfSrc.print_header() against expected resulting string.
        We use a file with a single track for this test, but only use the
        header of this file.
        """
        ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
        ssr.read_header()
        # If comparison output needs to be updated, uncomment the below
        #  and do: nosetests test_mcnp.py --nocapture
        #print ssr.print_header()
        self.assertEqual(ssr.print_header(),
                "Code: mcnp     (version: 5    ) [01232009]\n" \
                "Problem info: ( 07/05/12 17:50:19 )  07/05/12 17:50:16 \n" \
                "c Test deck with H20 cube, point n source, SSW of top surface interactions      \n" \
                "Showing dump #2\n" \
                "1 histories, 1 tracks, 11 record size, 1 surfaces, 1 histories\n" \
                "0 cells, source particle: 3, macrobody facet flag: 0\n" \
                "Surface [6]: facet -1, type [4] with 1 parameters: ( [5.0])\n" \
                "Summary Table: [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]"
            )
 
        return

    def test_print_tracklist(self):
        """Check SurfSrc.print_tracklist() against expected resulting string.
        We use a file with a single track for this test.
        """
        ssr = mcnp.SurfSrc(ssrname_onetrack, "rb")
        ssr.read_header()
        ssr.read_tracklist()
        # If comparison output needs to be updated, uncomment the below
        #  and do: nosetests test_mcnp.py --nocapture
        #print ssr.print_tracklist()
        self.assertEqual(ssr.print_tracklist(), 'Track Data\n       nps   BITARRAY        WGT        ERG        TME             X             Y             Z          U          V     COSINE  |       W\n         1 8.00005e+06    0.99996      5.542    0.17144  -8.05902e-02   3.12267e+00   5.00000e+00   -0.35133    0.48465    0.80105  |    0.80105 \n')

        return

def test_read_mcnp():

    expected_material = Material(nucvec={922350: 0.04, 922380: 0.96}, 
                                 mass=-1.0, 
                                 density=19.1, 
                                 attrs={"comments":(
                                            " first line of comments second line of "
                                            "comments third line of comments forth "
                                            "line of comments"),
                                        "mat_number": "1",
                                        "name":" leu", 
                                        "source":" Some http://URL.com",
                                        "table_ids":{'922350':"15c"}})
    expected_material.mass = -1.0  # to avoid reassignment to +1.0

    expected_multimaterial = MultiMaterial({
        Material({10000: 0.11190248274452597, 80000: 0.888097517255474}, 
                    -1.0, 0.9, 3, {"comments": (" Here are comments the comments "
                                                "continue here are more even more"),
                                   "mat_number": "2", 
                                   "name":" water",
                                   "source":" internet",
                                   "table_ids":{'10000':"05c"}}): 1,
        Material({10000: 0.11190248274452597, 80000: 0.888097517255474}, -1.0, 
                 1.002153632715214, 3, {"comments": (" Here are comments the comments "
                                            "continue here are more even more"),
                                        "mat_number": "2", 
                                        "name": " water",
                                        "source": " internet",
                                        "table_ids":{'10000':"05c"}}): 1})

    read_materials = read_mcnp_inp('mcnp_inp.txt')
    assert_equal(expected_material, read_materials[0])
    assert_equal(expected_multimaterial._mats.keys()[0].comp,
                                        read_materials[1]._mats.keys()[0].comp)
    assert_equal(expected_multimaterial._mats.keys()[0].mass,
                                        read_materials[1]._mats.keys()[0].mass)
    assert_equal(expected_multimaterial._mats.keys()[0].density,
                                      read_materials[1]._mats.keys()[0].density)
    assert_equal(expected_multimaterial._mats.keys()[0].atoms_per_mol,
                                read_materials[1]._mats.keys()[0].atoms_per_mol)
    assert_equal(expected_multimaterial._mats.keys()[0].attrs,
                                      read_materials[1]._mats.keys()[0].attrs)
    assert_equal(expected_multimaterial._mats.keys()[1].comp,
                                        read_materials[1]._mats.keys()[1].comp)
    assert_equal(expected_multimaterial._mats.keys()[1].mass,
                                        read_materials[1]._mats.keys()[1].mass)
    assert_equal(expected_multimaterial._mats.keys()[1].density,
                                      read_materials[1]._mats.keys()[1].density)
    assert_equal(expected_multimaterial._mats.keys()[1].atoms_per_mol,
                                read_materials[1]._mats.keys()[1].atoms_per_mol)
    assert_equal(expected_multimaterial._mats.keys()[1].attrs,
                                      read_materials[1]._mats.keys()[1].attrs)


# test PtracReader class
class TestPtrac(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_headers(self):
        p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
        assert_equal(p.problem_title,
                "Generate a well-defined PTRAC file for PyNE test cases")
        del p

        # 8-byte ints, little endian
        p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
        assert_equal(p.problem_title,
                "Generate a well-defined PTRAC file for PyNE test cases")
        del p

    def test_determine_format(self):
        # 4-byte ints, little endian
        p = mcnp.PtracReader("mcnp_ptrac_i4_little.ptrac")
        assert_equal(p.endianness, "<")
        del p

        # 8-byte ints, little endian
        p = mcnp.PtracReader("mcnp_ptrac_i8_little.ptrac")
        assert_equal(p.endianness, "<")
        del p

    def test_read_events(self):
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

    def test_write_to_hdf5(self):
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


# Test Wwinp class. All three function are tested at once because there inputs 
# and ouputs are easily strung together.
def test_Wwinp_n():

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_n.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_n.h5m')
    expected_sm = ScdMesh.fromFile(expected_h5m)
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
    assert_equal(ww1.nwg , 1)
    assert_equal(ww1.cm , [[-99, -97,  97, 99, 100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1, 1, 11, 1, 1], [1, 3, 4], [6]])
    assert_array_equal(ww1.e, [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]])
    assert_equal(ww1.bounds, [[-100.0, -99.0, -97.0, -79.36363636363636, 
                             -61.727272727272727, -44.090909090909093, 
                             -26.454545454545453, -8.818181818181813, 
                               8.818181818181813, 26.454545454545453, 
                              44.090909090909093, 61.72727272727272, 
                              79.363636363636374, 97.0, 99.0, 100.0],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]])    
    for x in range(0,15):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_n_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww1.mesh.getHex(x,y,z)
                    written = ww1.mesh.imesh\
                                 .getTagHandle('ww_n_group_00{0}'\
                                 .format(e_group))[written_voxel]
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
    assert_equal(ww2.nwg , 1)
    assert_equal(ww2.cm , [[-99, -97,  97, 99, 100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1, 1, 11, 1, 1], [1, 3, 4], [6]])
    assert_array_equal(ww2.e, [[0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000]])
    assert_equal(ww2.bounds, [[-100.0, -99.0, -97.0, -79.36363636363636, 
                             -61.727272727272727, -44.090909090909093, 
                             -26.454545454545453, -8.818181818181813, 
                               8.818181818181813, 26.454545454545453, 
                              44.090909090909093, 61.72727272727272, 
                              79.363636363636374, 97.0, 99.0, 100.0],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0,15):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_n_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww2.mesh.getHex(x,y,z)
                    written = ww2.mesh.imesh\
                              .getTagHandle('ww_n_group_00{0}'\
                              .format(e_group))[written_voxel]
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
            assert_equal(float(written[i].split()[j]), float(expected[i].split()[j]))

    os.remove(output)


def test_Wwinp_p():

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_p.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_p.h5m')
    expected_sm = ScdMesh.fromFile(expected_h5m)
    output = os.path.join(os.getcwd(), 'test_wwinp')

    
    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert_equal(ww1.ni, 2)
    assert_equal(ww1.nr, 10)
    assert_equal(ww1.ne,  [0,7])
    assert_equal(ww1.nf, [1, 8, 6])
    assert_equal(ww1.origin, [-100, -100, -100])
    assert_equal(ww1.nc, [1, 3, 1])
    assert_equal(ww1.nwg , 1)
    assert_equal(ww1.cm , [[100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1], [1, 3, 4], [6]])
    assert_equal(ww1.e[0], [])
    assert_array_equal(ww1.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_equal(ww1.bounds, [[-100.0, 100],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]]) 
    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_p_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww1.mesh.getHex(x,y,z)
                    written = ww1.mesh.imesh\
                                 .getTagHandle('ww_p_group_00{0}'\
                                 .format(e_group))[written_voxel]
                    assert_equal(written, expected)


    # Create an new object based off of only the mesh attribute of the first 
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert_equal(ww2.ni, 2)
    assert_equal(ww2.nr, 10)
    assert_equal(ww2.ne,  [0,7])
    assert_equal(ww2.nf, [1, 8, 6])
    assert_equal(ww2.origin, [-100, -100, -100])
    assert_equal(ww2.nc, [1, 3, 1])
    assert_equal(ww2.nwg , 1)
    assert_equal(ww2.cm , [[100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1], [1, 3, 4], [6]])
    assert_equal(ww2.e[0], [])
    assert_array_equal(ww2.e[1], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_equal(ww2.bounds, [[-100.0, 100],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_p_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww2.mesh.getHex(x,y,z)
                    written = ww2.mesh.imesh\
                              .getTagHandle('ww_p_group_00{0}'\
                              .format(e_group))[written_voxel]
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
            assert_equal(float(written[i].split()[j]), float(expected[i].split()[j]))

    os.remove(output)


def test_Wwinp_np():

    thisdir = os.path.dirname(__file__)
    wwinp_file = os.path.join(thisdir, 'mcnp_wwinp_wwinp_np.txt')
    expected_h5m = os.path.join(thisdir, 'mcnp_wwinp_mesh_np.h5m')
    expected_sm = ScdMesh.fromFile(expected_h5m)
    output = os.path.join(os.getcwd(), 'test_wwinp')

    
    # Read in the wwinp file to an object and check resulting attributes.
    ww1 = mcnp.Wwinp()
    ww1.read_wwinp(wwinp_file)
    assert_equal(ww1.ni, 2)
    assert_equal(ww1.nr, 10)
    assert_equal(ww1.ne,  [7,1])
    assert_equal(ww1.nf, [1, 8, 6])
    assert_equal(ww1.origin, [-100, -100, -100])
    assert_equal(ww1.nc, [1, 3, 1])
    assert_equal(ww1.nwg , 1)
    assert_equal(ww1.cm , [[100],  [-50, 60, 100], [100]])
    assert_equal(ww1.fm, [[1], [1, 3, 4], [6]])
    assert_equal(ww1.e[0], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_array_equal(ww1.e[1], [100])
    assert_equal(ww1.bounds, [[-100.0, 100],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]]) 
    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_n_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww1.mesh.getHex(x,y,z)
                    written = ww1.mesh.imesh\
                                 .getTagHandle('ww_n_group_00{0}'\
                                 .format(e_group))[written_voxel]
                    assert_equal(written, expected)

    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 2):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_p_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww1.mesh.getHex(x,y,z)
                    written = ww1.mesh.imesh\
                                 .getTagHandle('ww_p_group_00{0}'\
                                 .format(e_group))[written_voxel]
                    assert_equal(written, expected)


    # Create an new object based off of only the mesh attribute of the first 
    # object and check resutling attributes.
    ww2 = mcnp.Wwinp()
    ww2.read_mesh(ww1.mesh)
    assert_equal(ww2.ni, 2)
    assert_equal(ww2.nr, 10)
    assert_equal(ww2.ne,  [7,1])
    assert_equal(ww2.nf, [1, 8, 6])
    assert_equal(ww2.origin, [-100, -100, -100])
    assert_equal(ww2.nc, [1, 3, 1])
    assert_equal(ww2.nwg , 1)
    assert_equal(ww2.cm , [[100],  [-50, 60, 100], [100]])
    assert_equal(ww2.fm, [[1], [1, 3, 4], [6]])
    assert_array_equal(ww2.e[0], [0.1, 0.14678, 0.21544, 0.31623, 0.46416, 0.68129, 1.0000])
    assert_array_equal(ww2.e[1], [100])
    assert_equal(ww2.bounds, [[-100.0, 100],
                             [-100.0, -50.0, -13.333333333333336, 
                              23.333333333333329, 60.0, 70.0, 80.0, 90.0, 100.0], 
                             [-100.0, -66.666666666666657, -33.333333333333329, 
                              0.0, 33.333333333333343, 66.666666666666657, 100.0]])
    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 8):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_n_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww2.mesh.getHex(x,y,z)
                    written = ww2.mesh.imesh\
                              .getTagHandle('ww_n_group_00{0}'\
                              .format(e_group))[written_voxel]
                    assert_equal(written, expected)

    for x in range(0,1):
        for y in range(0,8):
            for z in range(0,6):
                for e_group in range(1, 2):
                    expected_voxel = expected_sm.getHex(x,y,z)
                    expected = expected_sm.imesh\
                               .getTagHandle('ww_p_group_00{0}'\
                               .format(e_group))[expected_voxel]
                    written_voxel = ww2.mesh.getHex(x,y,z)
                    written = ww2.mesh.imesh\
                              .getTagHandle('ww_p_group_00{0}'\
                              .format(e_group))[written_voxel]
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
            assert_equal(float(written[i].split()[j]), float(expected[i].split()[j]))

    os.remove(output)

