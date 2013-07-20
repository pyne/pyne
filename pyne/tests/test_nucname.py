"""nucname tests"""
from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in

from pyne import nucname

def test_name_zz():
    assert_equal(nucname.name_zz["He"], 2)
    assert_equal(nucname.name_zz["U"], 92)

def test_zz_name():
    assert_equal(nucname.zz_name[1],  "H")
    assert_equal(nucname.zz_name[94], "Pu")


def test_LAN():
    assert_equal(nucname.LAN, set(['Ce', 'Dy', 'Er', 'Eu', 'Gd', 'Ho', 'La', 
                                   'Lu', 'Nd', 'Pm', 'Pr', 'Sm', 'Tb', 'Tm', 'Yb']))

def test_ACT():
    assert_equal(nucname.ACT, set(["Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", 
                                   "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]))
            
def test_TRU():
    assert_equal(nucname.TRU, set(['Am', 'Bh', 'Bk', 'Cf', 'Cm', 'Db', 'Ds', 
     'Es', 'Fm', 'Hs', 'Lr', 'Md', 'Mt', 'No', 'Np', 'Pu', 'Rf', 'Rg', 'Sg', 
     'Fl', 'Lv', 'Cn']))

def test_MA():
    assert_equal(nucname.MA, set(["Np", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", 
                                  "No", "Lr"]))

def test_FP():
    assert_equal(nucname.FP, set(['Ag', 'Al', 'Ar', 'As', 'At', 'Au', 'B',  
          'Ba', 'Be', 'Bi', 'Br', 'C',  'Ca', 'Cd', 'Ce', 'Cl', 'Co', 'Cr', 
          'Cs', 'Cu', 'Dy', 'Er', 'Eu', 'F',  'Fe', 'Fr', 'Ga', 'Gd', 'Ge',
          'H',  'He', 'Hf', 'Hg', 'Ho', 'I',  'In', 'Ir', 'K',  'Kr', 'La', 
          'Li', 'Lu', 'Mg', 'Mn', 'Mo', 'N',  'Na', 'Nb', 'Nd', 'Ne', 'Ni', 
          'O',  'Os', 'P',  'Pb', 'Pd', 'Pm', 'Po', 'Pr', 'Pt', 'Ra', 'Rb', 
          'Re', 'Rh', 'Rn', 'Ru', 'S',  'Sb', 'Sc', 'Se', 'Si', 'Sm', 'Sn', 
          'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Ti', 'Tl', 'Tm', 'V',  'W',  'Xe',  
          'Y',  'Yb', 'Zn', 'Zr']))


def test_lan():
    assert_equal(nucname.lan, set([57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 
                                   68, 69, 70, 71]))

def test_act():
    assert_equal(nucname.act, set([89,  90,  91,  92, 93, 94, 95, 96, 97, 98, 99, 
                                   100, 101, 102, 103]))

def test_tru():
    assert_equal(nucname.tru, set([93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 
                         103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 114, 116]))

def test_ma():
    assert_equal(nucname.ma, set([93, 95, 96, 97, 98, 99, 100, 101, 102, 103]))

def test_fp():
    assert_equal(nucname.fp, set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
            14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
            31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 
            65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
            82, 83, 84, 85, 86, 87, 88]))



def test_id():
    assert_equal(nucname.id(20040), 20040000)

    assert_equal(nucname.id("he4"),    20040000)
    assert_equal(nucname.id("Cm-244"), 962440000)
    assert_equal(nucname.id("PU239"),  942390000)
    assert_equal(nucname.id("AM242M"), 952420001)

    assert_equal(nucname.id(2004),  20040000)
    assert_equal(nucname.id(95642), 952420000)
    assert_equal(nucname.id(95242), 952420001)
    assert_equal(nucname.id(92636), 922360001)
    assert_equal(nucname.id(95942), 952420004)

    assert_equal(nucname.id("Am-242m"), 952420001)

    assert_equal(nucname.id("he"), 20000000)
    assert_equal(nucname.id("U"), 920000000)
    assert_equal(nucname.id("Np"), 930000000)

    assert_equal(nucname.id("4he"),   20040000)
    assert_equal(nucname.id("244CM"), 962440000)
    assert_equal(nucname.id("239Pu"), 942390000)
    assert_equal(nucname.id("242AM"), 952420000)

    assert_equal(nucname.id(40020),   20040000)
    assert_equal(nucname.id(2440961), 962440001)
    assert_equal(nucname.id(2390940), 942390000)
    assert_equal(nucname.id(2420950), 952420000)


def test_name():
    assert_equal(nucname.name(942390), "Pu239")
    assert_equal(nucname.name(952421), "Am242M")

    assert_equal(nucname.name("PU239"), "Pu239")

    assert_equal(nucname.name(94239), "Pu239")
    assert_equal(nucname.name(95242), "Am242M")
    assert_equal(nucname.name(95642), "Am242")
    assert_equal(nucname.name(92636), "U236M")

    assert_equal(nucname.name("Am-242m"), "Am242M")

    assert_equal(nucname.name(40020),   "He4")
    assert_equal(nucname.name(2440961), "Cm244M")
    assert_equal(nucname.name(2390940), "Pu239")
    assert_equal(nucname.name(2420950), "Am242")


def test_zzaaam():
    assert_equal(nucname.zzaaam(20040), 20040)

    assert_equal(nucname.zzaaam("he4"),    20040)
    assert_equal(nucname.zzaaam("Cm-244"), 962440)
    assert_equal(nucname.zzaaam("PU239"),  942390)
    assert_equal(nucname.zzaaam("AM242M"), 952421)

    assert_equal(nucname.zzaaam(2004),  20040)
    assert_equal(nucname.zzaaam(95642), 952420)
    assert_equal(nucname.zzaaam(95242), 952421)
    assert_equal(nucname.zzaaam(92636), 922361)
    assert_equal(nucname.zzaaam(95942), 952424)

    assert_equal(nucname.zzaaam("Am-242m"), 952421)

    assert_equal(nucname.zzaaam("he"), 20000)
    assert_equal(nucname.zzaaam("U"), 920000)
    assert_equal(nucname.zzaaam("Np"), 930000)

    assert_equal(nucname.zzaaam("4he"),   20040)
    assert_equal(nucname.zzaaam("244CM"), 962440)
    assert_equal(nucname.zzaaam("239Pu"), 942390)
    assert_equal(nucname.zzaaam("242AM"), 952420)

    assert_equal(nucname.zzaaam(40020),   20040)
    assert_equal(nucname.zzaaam(2440961), 962441)
    assert_equal(nucname.zzaaam(2390940), 942390)
    assert_equal(nucname.zzaaam(2420950), 952420)


def test_mcnp():
    assert_equal(nucname.mcnp(10010),  1001)
    assert_equal(nucname.mcnp(952421), 95242)
    assert_equal(nucname.mcnp(952420), 95642)
    assert_equal(nucname.mcnp(922361), 92636)

    assert_equal(nucname.mcnp("H1"),     1001)
    assert_equal(nucname.mcnp("AM242M"), 95242)
    assert_equal(nucname.mcnp("AM242"), 95642)
    assert_equal(nucname.mcnp("U236M"), 92636)

    assert_equal(nucname.mcnp(1001),  1001)
    assert_equal(nucname.mcnp(95642), 95642)

    assert_equal(nucname.mcnp("Am-242"),  95642)
    assert_equal(nucname.mcnp("Am-242m"), 95242)
    assert_equal(nucname.mcnp("U-236m"),  92636)

    assert_equal(nucname.mcnp(10010),  1001)
    assert_equal(nucname.mcnp(2420951), 95242)
    assert_equal(nucname.mcnp(2420950), 95642)
    assert_equal(nucname.mcnp(2360921), 92636)


def test_serpent():
    assert_equal(nucname.serpent(942390), "Pu-239")
    assert_equal(nucname.serpent(952421), "Am-242m")

    assert_equal(nucname.serpent("Pu-239"), "Pu-239")

    assert_equal(nucname.serpent(94239), "Pu-239")
    assert_equal(nucname.serpent(95642), "Am-242")
    assert_equal(nucname.serpent(95242), "Am-242m")
    assert_equal(nucname.serpent(92636), "U-236m")

    assert_equal(nucname.serpent(2390940), "Pu-239")
    assert_equal(nucname.serpent(2420951), "Am-242m")



def test_nist():
    assert_equal(nucname.nist(942390), "239Pu")
    assert_equal(nucname.nist(952421), "242Am")

    assert_equal(nucname.nist("Pu-239"), "239Pu")

    assert_equal(nucname.nist(94239), "239Pu")
    assert_equal(nucname.nist(95642), "242Am")

    assert_equal(nucname.nist("10000"), "H")
    assert_equal(nucname.nist("920000"), "U")
    assert_equal(nucname.nist("940000"), "Pu")

    assert_equal(nucname.nist(2390940), "239Pu")
    assert_equal(nucname.nist(2420951), "242Am")


def test_cinder():
    assert_equal(nucname.cinder(10020),  20010)
    assert_equal(nucname.cinder(952421), 2420951)

    assert_equal(nucname.cinder("H2"),     20010)
    assert_equal(nucname.cinder("AM242M"), 2420951)

    assert_equal(nucname.cinder(1002),  20010)
    assert_equal(nucname.cinder(95642), 2420950)
    assert_equal(nucname.cinder(95242), 2420951)
    assert_equal(nucname.cinder(92636), 2360921)

    assert_equal(nucname.cinder("Am-242m"), 2420951)

    assert_equal(nucname.cinder(20010),  20010)
    assert_equal(nucname.cinder(2420951), 2420951)


def test_alara():
    assert_equal(nucname.alara(942390), "pu:239")
    assert_equal(nucname.alara(952421), "am:242")

    assert_equal(nucname.alara("PU239"), "pu:239")

    assert_equal(nucname.alara(94239), "pu:239")
    assert_equal(nucname.alara(95242), "am:242")
    assert_equal(nucname.alara(95642), "am:242")
    assert_equal(nucname.alara(92636), "u:236")
    assert_equal(nucname.alara(2000),  "he")

    assert_equal(nucname.alara("Am-242m"), "am:242")

    assert_equal(nucname.alara(40020),   "he:4")
    assert_equal(nucname.alara(20000),   "he")
    assert_equal(nucname.alara(2440961), "cm:244")
    assert_equal(nucname.alara(2390940), "pu:239")
    assert_equal(nucname.alara(2420950), "am:242")



if __name__ == "__main__":
    nose.main()

