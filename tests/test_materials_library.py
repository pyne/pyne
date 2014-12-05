from pyne.dbgen.materials_library import make_elements
from pyne.dbgen.materials_library import grab_materials_compendium
from pyne.dbgen.materials_library import make_materials_compendium
import pyne.nucname as nucname

import tables as tb


def assert_close(obs, exp, margin=1e-6):
    # default margin is because we only have 1e-6 of precision
    assert(abs(obs-exp) < margin)


def test_grab_materials_compendium():
    mats = grab_materials_compendium('../pyne/dbgen/materials_compendium.csv')
    assert(len(mats) == 372)
    a150tep = mats[0]
    # gotta normalize our expectation!
    assert_close(a150tep.comp[nucname.id(1001)], (0.101327 / 1.000001))
    assert_close(a150tep.comp[nucname.id(6000)], (0.775501 / 1.000001))
    assert_close(a150tep.comp[nucname.id(7014)], (0.035057 / 1.000001))
    assert_close(a150tep.comp[nucname.id(8016)], (0.052316 / 1.000001))
    assert_close(a150tep.comp[nucname.id(9019)], (0.017422 / 1.000001))
    assert_close(a150tep.comp[nucname.id(20000)], (0.018378 / 1.000001))

    pubr = [mat for mat in mats if mat.metadata["name"] == "Plutonium Bromide"][0]
    bromium = sum((frac for nuc, frac in pubr.comp.items() if nucname.zzzaaa(nuc) // 1000 == 35))
    assert_close(bromium, 0.500617)


# def test_make_materials_compendium():
#     elts = make_elements()
#     mats = grab_materials_compendium(os.path.join(os.path.split(__file__)[0],
#                                                   '../pyne/dbgen/materials_compendium.csv'))
#     make_materials_compendium('test_materials_library.h5', mats, elts)
#     with tb.openFile('test_materials_library.h5', 'r') as f:
#         table = f.root.detector.readout
#         assert(False)
