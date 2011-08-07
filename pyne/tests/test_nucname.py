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


def test_LAN():
    assert_equal(nucname.LAN, set(['CE', 'DY', 'ER', 'EU', 'GD', 'HO', 'LA', 
                                   'LU', 'ND', 'PM', 'PR', 'SM', 'TB', 'TM', 'YB']))

def test_ACT():
    assert_equal(nucname.ACT, set(["AC", "TH", "PA", "U",  "NP", "PU", "AM", "CM", \
                                   "BK", "CF", "ES", "FM", "MD", "NO", "LR"]))
            
def test_TRU():
    assert_equal(nucname.TRU, set(['AM', 'BH', 'BK', 'CF', 'CM', 'DB', 'DS', 
     'ES', 'FM', 'HS', 'LR', 'MD', 'MT', 'NO', 'NP', 'PU', 'RF', 'RG', 'SG']))

def test_MA():
    assert_equal(nucname.MA, set(["NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", \
                                  "NO", "LR"]))

def test_FP():
    assert_equal(nucname.FP, set(['AG', 'AL', 'AR', 'AS', 'AT', 'AU', 'B',  
          'BA', 'BE', 'BI', 'BR', 'C',  'CA', 'CD', 'CE', 'CL', 'CO', 'CR', 
          'CS', 'CU', 'DY', 'ER', 'EU', 'F',  'FE', 'FR', 'GA', 'GD', 'GE',
          'H',  'HE', 'HF', 'HG', 'HO', 'I',  'IN', 'IR', 'K',  'KR', 'LA', 
          'LI', 'LU', 'MG', 'MN', 'MO', 'N',  'NA', 'NB', 'ND', 'NE', 'NI', 
          'O',  'OS', 'P',  'PB', 'PD', 'PM', 'PO', 'PR', 'PT', 'RA', 'RB', 
          'RE', 'RH', 'RN', 'RU', 'S',  'SB', 'SC', 'SE', 'SI', 'SM', 'SN', 
          'SR', 'TA', 'TB', 'TC', 'TE', 'TI', 'TL', 'TM', 'V',  'W',  'XE',  
          'Y',  'YB', 'ZN', 'ZR']))


def test_lan():
    assert_equal(nucname.lan, set([57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 
                                   68, 69, 70, 71]))

def test_act():
    assert_equal(nucname.act, set([89,  90,  91,  92, 93, 94, 95, 96, 97, 98, 99, 
                                   100, 101, 102, 103]))

def test_tru():
    assert_equal(nucname.tru, set([93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 
                                   103, 104, 105, 106, 107, 108, 109, 110, 111]))

def test_ma():
    assert_equal(nucname.ma, set([93, 95, 96, 97, 98, 99, 100, 101, 102, 103]))

def test_fp():
    assert_equal(nucname.fp, set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
            14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
            31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 
            65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
            82, 83, 84, 85, 86, 87, 88]))


def test_current_form():
    assert_equal(nucname.current_form(922380),   "zzaaam")
    assert_equal(nucname.current_form("922350"), "zzaaam")

    assert_equal(nucname.current_form("U235"),  "LLAAAM")
    assert_equal(nucname.current_form("h-004"), "LLAAAM")

    assert_equal(nucname.current_form(92238),   "MCNP")
    assert_equal(nucname.current_form("92235"), "MCNP")



def test_zzaaam():
    assert_equal(nucname.zzaaam(20040), 20040)

    assert_equal(nucname.zzaaam("he4"),    20040)
    assert_equal(nucname.zzaaam("Cm-244"), 962440)
    assert_equal(nucname.zzaaam("PU239"),  942390)
    assert_equal(nucname.zzaaam("AM242M"), 952421)

    assert_equal(nucname.zzaaam(2004),  20040)
    assert_equal(nucname.zzaaam(95642), 952421)
    assert_equal(nucname.zzaaam(95942), 952424)

    assert_equal(nucname.zzaaam("Am-242m"), 952421)

    assert_equal(nucname.zzaaam("he"), 20000)
    assert_equal(nucname.zzaaam("U"), 920000)
    assert_equal(nucname.zzaaam("Np"), 930000)


def test_LLAAAM():
    assert_equal(nucname.LLAAAM(942390), "PU239")
    assert_equal(nucname.LLAAAM(952421), "AM242M")

    assert_equal(nucname.LLAAAM("PU239"), "PU239")

    assert_equal(nucname.LLAAAM(94239), "PU239")
    assert_equal(nucname.LLAAAM(95642), "AM242M")

    assert_equal(nucname.LLAAAM("Am-242m"), "AM242M")


def test_mcnp():
    assert_equal(nucname.mcnp(10010),  1001)
    assert_equal(nucname.mcnp(952421), 95642)

    assert_equal(nucname.mcnp("H1"),     1001)
    assert_equal(nucname.mcnp("AM242M"), 95642)

    assert_equal(nucname.mcnp(1001),  1001)
    assert_equal(nucname.mcnp(95642), 95642)

    assert_equal(nucname.mcnp("Am-242m"), 95642)



def test_serpent():
    assert_equal(nucname.serpent(942390), "Pu-239")
    assert_equal(nucname.serpent(952421), "Am-242m")

    assert_equal(nucname.serpent("Pu-239"), "Pu-239")

    assert_equal(nucname.serpent(94239), "Pu-239")
    assert_equal(nucname.serpent(95642), "Am-242m")



def test_nuc_weight():
    # zzaam form
    assert_equal(nucname.nuc_weight(80160),  16.0)
    assert_equal(nucname.nuc_weight(922350), 235.0)
    assert_equal(nucname.nuc_weight(952421), 242.0)

    # MCNP form
    assert_equal(nucname.nuc_weight(8016),  16.0)
    assert_equal(nucname.nuc_weight(92235), 235.0)
    assert_equal(nucname.nuc_weight(95642), 242.0)

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

