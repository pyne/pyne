"""PyNE MCNP tools tests"""
import os
import unittest
import nose

import nose.tools

from pyne import mcnp

thisdir = os.path.dirname(__file__)
ssrname = os.path.join(thisdir,"mcnp_surfsrc.w")
sswname = os.path.join(thisdir,"copy_mcnp_surfsrc.w")

class TestSurfSrc(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_read_header_block(self):
        """
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
        # surface info recod list
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
