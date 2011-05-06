"""nucname tests"""
from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises

from pyne import nucname

def test_LLzz():
    assert_equal(nucname.LLzz["HE"], 2)
    assert_equal(nucname.LLzz["U"], 92)

def test_zzLL():
    assert_equal(nucname.zzLL[1],  "H")
    assert_equal(nucname.zzLL[94], "PU")


class TestElementGroupsLL(TestCase):
    """Tests the LL element groups are right."""

    def test_LAN(self):
        assert_equal(nucname.LAN, set(['CE', 'DY', 'ER', 'EU', 'GD', 'HO', 'LA', 
               'LU', 'ND', 'PM', 'PR', 'SM', 'TB', 'TM', 'YB']))

    def test_ACT(self):
        assert_equal(nucname.ACT, set(['AC', 'AM', 'BH', 'BK', 'CF', 'CM', 'DB', 
               'DS', 'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 'PA', 'PU', 
               'RF', 'RG', 'SG', 'TH', 'U']))
            
    def test_TRU(self):
        assert_equal(nucname.TRU, set(['AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 
               'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 'PU', 'RF', 'RG', 
               'SG']))

    def test_MA(self):
        assert_equal(nucname.MA, set(['AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 
              'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 'RF', 'RG', 'SG']))

    def test_FP(self):
        assert_equal(nucname.FP, set(['AG', 'AL', 'AR', 'AS', 'AT', 'AU', 'B',  
              'BA', 'BE', 'BI', 'BR', 'C',  'CA', 'CD', 'CE', 'CL', 'CO', 'CR', 
              'CS', 'CU', 'DY', 'ER', 'EU', 'F',  'FE', 'FR', 'GA', 'GD', 'GE',
              'H',  'HE', 'HF', 'HG', 'HO', 'I',  'IN', 'IR', 'K',  'KR', 'LA', 
              'LI', 'LU', 'MG', 'MN', 'MO', 'N',  'NA', 'NB', 'ND', 'NE', 'NI', 
              'O',  'OS', 'P',  'PB', 'PD', 'PM', 'PO', 'PR', 'PT', 'RA', 'RB', 
              'RE', 'RH', 'RN', 'RU', 'S',  'SB', 'SC', 'SE', 'SI', 'SM', 'SN', 
              'SR', 'TA', 'TB', 'TC', 'TE', 'TI', 'TL', 'TM', 'V',  'W',  'XE',  
              'Y', 'YB', 'ZN', 'ZR']))


class TestElementGroupsZz(TestCase):
    """Tests the zz element groups are right."""

    def test_lan(self):
        assert_equal(nucname.lan, set([57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 
            68, 69, 70, 71]))

    def test_act(self):
        assert_equal(nucname.act, set([89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 
            100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111]))

    def test_tru(self):
        assert_equal(nucname.tru, set([93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 
            103, 104, 105, 106, 107, 108, 109, 110, 111]))

    def test_ma(self):
        assert_equal(nucname.ma, set([93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 
            104, 105, 106, 107, 108, 109, 110, 111]))

    def test_fp(self):
        assert_equal(nucname.fp, set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
            14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
            31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 
            65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
            82, 83, 84, 85, 86, 87, 88]))


class TestZzaaamFunctions(TestCase):
    """Tests the zzaaam_2_* conversion functions."""

    def test_zzaaam_2_LLAAAM(self):
        assert_equal(nucname.zzaaam_2_LLAAAM(942390), "PU239")
        assert_equal(nucname.zzaaam_2_LLAAAM(952421), "AM242M")

    # Need to supress sdtout from C to work nicely
    #@raises(RuntimeError)
    #def test_zzaaam_2_LLAAAM_2(self):
    #    nucname.zzaaam_2_LLAAAM(94239)

    def test_zzaaam_2_MCNP(self):
        assert_equal(nucname.zzaaam_2_MCNP(10010),  1001)
        assert_equal(nucname.zzaaam_2_MCNP(952421), 95642)

    def test_zzaaam_2_LLAAAM_List(self):
        assert_equal(nucname.zzaaam_2_LLAAAM_List([942390, 80160]), 
            ["PU239", "O16"])

    def test_zzaaam_2_MCNP_List(self):
        assert_equal(nucname.zzaaam_2_MCNP_List([10010, 932370]), 
            [1001, 93237])


class TestLLAAAMFunctions(TestCase):
    """Tests the LLAAAM_2_* conversion functions."""

    def test_LLAAAM_2_zzaaam(self):
        assert_equal(nucname.LLAAAM_2_zzaaam("he4"),    20040)
        assert_equal(nucname.LLAAAM_2_zzaaam("Cm-244"), 962440)
        assert_equal(nucname.LLAAAM_2_zzaaam("PU239"),  942390)
        assert_equal(nucname.LLAAAM_2_zzaaam("AM242M"), 952421)

    def test_LLAAAM_2_MCNP(self):
        assert_equal(nucname.LLAAAM_2_MCNP("H1"),     1001)
        assert_equal(nucname.LLAAAM_2_MCNP("AM242M"), 95642)

    def test_LLAAAM_2_zzaaam_List(self):
        assert_equal(nucname.LLAAAM_2_zzaaam_List(["PU239", "O16"]), 
            [942390, 80160])

    def test_LLAAAM_2_MCNP_List(self):
        assert_equal(nucname.LLAAAM_2_MCNP_List(["H1", "NP237"]), 
            [1001, 93237])


class TestMCNPFunctions(TestCase):
    """Tests the MCNP_2_* conversion functions."""

    def test_MCNP_2_zzaaam(self):
        assert_equal(nucname.MCNP_2_zzaaam(2004),  20040)
        assert_equal(nucname.MCNP_2_zzaaam(95642), 952421)

    def test_MCNP_2_LLAAAM(self):
        assert_equal(nucname.MCNP_2_LLAAAM(94239), "PU239")
        assert_equal(nucname.MCNP_2_LLAAAM(95642), "AM242M")

    def test_MCNP_2_zzaaam_List(self):
        assert_equal(nucname.MCNP_2_zzaaam_List([94239, 8016]), 
            [942390, 80160])

    def test_MCNP_2_LLAAAM_List(self):
        assert_equal(nucname.MCNP_2_LLAAAM_List([94239, 8016]), 
            ["PU239", "O16"])

class TestCurrentForm(TestCase):
    """Tests the CurrentForm flag."""

    def test_CurrentForm_zzaaam(self):
        assert_equal(nucname.CurrentForm(922380),   "zzaaam")
        assert_equal(nucname.CurrentForm("922350"), "zzaaam")

    def test_CurrentForm_LLAAAM(self):
        assert_equal(nucname.CurrentForm("U235"),  "LLAAAM")
        assert_equal(nucname.CurrentForm("h-004"), "LLAAAM")

    def test_CurrentForm_MCNP(self):
        assert_equal(nucname.CurrentForm(92238),   "MCNP")
        assert_equal(nucname.CurrentForm("92235"), "MCNP")

class TestRearRemoveDuplicates(TestCase):
    """Tests the RearRemoveDuplicates function."""

    def test_RearRemoveDuplicates(self):
        assert_equal(nucname.RearRemoveDuplicates([6, 1, 2, 2, 5, 6]), 
            [6, 1, 2, 5])

class TestMixedFunctions(TestCase):
    """Tests the mixed_2_* conversion functions."""

    def test_mixed_2_zzaaam(self):
        assert_equal(nucname.mixed_2_zzaaam(20040),   20040)
        assert_equal(nucname.mixed_2_zzaaam("PU239"), 942390)
        assert_equal(nucname.mixed_2_zzaaam(95642),   952421)

    def test_mixed_2_LLAAAM(self):
        assert_equal(nucname.mixed_2_LLAAAM(952421),  "AM242M")
        assert_equal(nucname.mixed_2_LLAAAM("PU239"), "PU239")
        assert_equal(nucname.mixed_2_LLAAAM(94239),   "PU239")

    def test_mixed_2_MCNP(self):
        assert_equal(nucname.mixed_2_MCNP(10010),  1001)
        assert_equal(nucname.mixed_2_MCNP("H1"),   1001)
        assert_equal(nucname.mixed_2_MCNP(95642),  95642)

    def test_mixed_2_zzaaam_List(self):
        assert_equal(nucname.mixed_2_zzaaam_List([20040, "PU239", 95642]),
            [20040, 942390, 952421])

    def test_mixed_2_LLAAAM_List(self):
        assert_equal(nucname.mixed_2_LLAAAM_List([952421, "PU239", 93237]),
            ["AM242M", "PU239", "NP237"])

    def test_mixed_2_MCNP_List(self):
        assert_equal(nucname.mixed_2_MCNP_List([10010, "HE4", 95642]),
            [1001, 2004, 95642])

class TestIsovecKeysFunctions(TestCase):
    """Tests the isovec_keys_2_* conversion functions."""

    d = {
        922350:  0.5, 
        "U-238": 0.5, 
        95642:   1.0, 
        }

    def test_isovec_keys_2_zzaaam(self):
        assert_equal(nucname.isovec_keys_2_zzaaam(self.d), {
            922350: 0.5, 
            922380: 0.5, 
            952421: 1.0,             
            })

    def test_isovec_keys_2_LLAAAM(self):
        assert_equal(nucname.isovec_keys_2_LLAAAM(self.d), {
            "U235":   0.5, 
            "U238":   0.5, 
            "AM242M": 1.0,             
            })

    def test_isovec_keys_2_MCNP(self):
        assert_equal(nucname.isovec_keys_2_MCNP(self.d), {
            92235: 0.5, 
            92238: 0.5, 
            95642: 1.0,             
            })

class TestHelperFunctions(TestCase):
    """Tests the nucname module helper functions."""

    def test_nuc_weight_zzaaam(self):
        assert_equal(nucname.nuc_weight_zzaaam(80160),  16.0)
        assert_equal(nucname.nuc_weight_zzaaam(922350), 235.0)
        assert_equal(nucname.nuc_weight_zzaaam(952421), 242.0)

    def test_nuc_weight_int(self):
        # zzaam form
        assert_equal(nucname.nuc_weight(80160),  16.0)
        assert_equal(nucname.nuc_weight(922350), 235.0)
        assert_equal(nucname.nuc_weight(952421), 242.0)

        # MCNP form
        assert_equal(nucname.nuc_weight(8016),  16.0)
        assert_equal(nucname.nuc_weight(92235), 235.0)
        assert_equal(nucname.nuc_weight(95642), 242.0)

    def test_nuc_weight_string(self):
        # zzaam form
        assert_equal(nucname.nuc_weight("80160"),  16.0)
        assert_equal(nucname.nuc_weight("922350"), 235.0)
        assert_equal(nucname.nuc_weight("952421"), 242.0)

        # zzaam form
        assert_equal(nucname.nuc_weight("O-16"),   16.0)
        assert_equal(nucname.nuc_weight("U235"),   235.0)
        assert_equal(nucname.nuc_weight("Am242M"), 242.0)

        # MCNP form
        assert_equal(nucname.nuc_weight("8016"),  16.0)
        assert_equal(nucname.nuc_weight("92235"), 235.0)
        assert_equal(nucname.nuc_weight("95642"), 242.0)


if __name__ == "__main__":
    nose.main()

