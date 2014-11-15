#!/usr/bin/python

from pyne import material
from pyne.material import Material, MaterialLibrary
import nose
from nose.tools import assert_equal, assert_raises
import dagmc_get_materials as gtag


"""
Existence/Absence of a graveyard group
"""
def test_graveyard_1():
    # 'mat:graveyard' group exists
    tag_1 = ['mat:graveyard', 'mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_1), [
                 ('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])
    # 'mat:Graveyard' group exists
    tag_2 = ['mat:Graveyard', 'mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_2), [
                 ('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])


def test_graveyard_2():
    # graveyard group is absent
    tag_3 = ['mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception, gtag.check_matname, tag_3)


def test_graveyard_3():
    # graveyard exists as 'graveyard'
    tag_4 = ['graveyard', 'mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_4), [
                 ('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])

'''test check_matname function'''
"""
test groups naming
"""
def test_group_1():
    # ':' is missing in 'matLead/rho:-11.35'
    tag_5 = ['mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'matLead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception, gtag.check_matname, tag_5)


def test_group_2():
    # density is in wrong format in 'mat:Mercury/rho:mercury'
    tag_6 = ['graveyard', 'mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:mercury']
    assert_raises(Exception, gtag.check_matname, tag_6)


def test_group_3():
    # material name is absent in 'mat:/rho:-0.001205'
    tag_7 = ['mat:graveyard', 'mat:/rho:-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception, gtag.check_matname, tag_7)


def test_group_4():
    # an error in the group name; "/" without a density in 'mat:Nitrogen/'
    tag_8 = ['mat:graveyard', 'mat:Nitrogen/', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception, gtag.check_matname, tag_8)


def test_group_5():
    # ':' is absent in the density part in 'mat:Nitrogen/rho-0.001205'
    tag_9 = ['mat:graveyard', 'mat:Nitrogen/rho-0.001205', 'tally_4.cell.flux.p',
             'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception, gtag.check_matname, tag_9)


def test_group_6():
    # no desnity provided
    tag_10 = ['mat:graveyard', 'mat:Nitrogen', 'tally_4.cell.flux.p',
              'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_equal(gtag.check_matname(tag_10), [('Nitrogen', ''), (
        'Steel, Stainless 321', ''), ('Lead', ''), ('Mercury', '')])


def test_group_7():
    # some densities are provided
    tag_11 = ['mat:graveyard', 'mat:Nitrogen', 'tally_4.cell.flux.p',
              'mat:Steel, Stainless 321', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_11), [('Nitrogen', ''), (
        'Steel, Stainless 321', ''), ('Lead', '-11.35'), ('Mercury', '-7.874')])


def test_group_8():
    # no density provided and material name is absent
    tag_12 = ['mat:graveyard', 'mat:', 'tally_4.cell.flux.p',
              'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_raises(Exception, gtag.check_matname, tag_12)


def test_group_9():
    # no density is provided and ':' is absent in 'matNitrogen'
    tag_13 = ['mat:graveyard', 'matNitrogen', 'tally_4.cell.flux.p',
              'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_raises(Exception, gtag.check_matname, tag_13)


def test_group_10():
    # 'mat' is absent from the group names
    tag_14 = ['graveyard', 'Nitrogen', 'tally_4.cell.flux.p',
              'Steel, Stainless 321', 'Lead', 'Mercury']
    assert_raises(Exception, gtag.check_matname, tag_14)


'''test fluka_material_naming function'''
"""
test attributes
"""


def test_fluka_1():
    # material doesn't exist before in the fluka materials list
    steel = material.Material({60120000: 0.10992222222222224, 60130000: 0.0011888888888888893, 140280000: 0.10247000000000002, 140290000: 0.005205555555555556, 140300000: 0.0034355555555555563, 150310000: 0.11111111111111112, 160320000: 0.10554444444444445, 160330000: 0.0008333333333333334, 160340000: 0.004722222222222223, 160360000: 1.1111111111111113e-05, 220460000: 0.009166666666666668, 220470000: 0.008266666666666669, 220480000: 0.08191111111111112, 220490000: 0.006011111111111112, 220500000: 0.005755555555555556,
                              240500000: 0.004827777777777778, 240520000: 0.09309888888888891, 240530000: 0.010556666666666667, 240540000: 0.002627777777777779, 250550000: 0.11111111111111112, 260540000: 0.006494444444444446, 260560000: 0.10194888888888891, 260570000: 0.0023544444444444455, 260580000: 0.0003133333333333333, 280580000: 0.07564111111111112, 280600000: 0.029136666666666672, 280610000: 0.0012665555555555557, 280620000: 0.004038444444444445, 280640000: 0.0010283333333333336}, 1.0, -7.0, -1.0, {"mat_number": "4", "name": "Steel, Stainless 321"})
    flukamat_list = ['NITROGEN', 'LEAD']
    mat = gtag.fluka_material_naming(steel, flukamat_list)
    original_name = mat.metadata['name']
    name = mat.metadata['fluka_name']
    assert_equal(name, 'STEELSTA')
    assert_equal(original_name, 'Steel, Stainless 321')


def test_fluka_2():
    # material exists in the fluka materials list
    nitrogen = material.Material(
        {70140000: 0.99636, 70150000: 0.00364}, 1.0, 0.001165, -1.0, {"mat_number": "1", "name": "Nitrogen"})
    flukamat_list = ['NITROGEN', 'LEAD']
    mat = gtag.fluka_material_naming(nitrogen, flukamat_list)
    original_name = mat.metadata['name']
    name = mat.metadata['fluka_name']
    assert_equal(name, 'NITROGE1')
    assert_equal(original_name, 'Nitrogen')


def test_fluka_3():
    # material exists twice in the the fluka materials list
    nitrogen = material.Material(
        {70140000: 0.99636, 70150000: 0.00364}, 1.0, 0.001165, -1.0, {"mat_number": "1", "name": "Nitrogen"})
    flukamat_list = ['NITROGEN', 'LEAD', 'NITROGE1']
    mat = gtag.fluka_material_naming(nitrogen, flukamat_list)
    original_name = mat.metadata['name']
    name = mat.metadata['fluka_name']
    assert_equal(name, 'NITROGE2')
    assert_equal(original_name, 'Nitrogen')


''' test print_near_match function'''
"""
test 
"""


def test_match_1():
    air_1 = material.Material({60120000: 0.24732500000000002, 60130000: 0.0026750000000000003, 70140000: 0.24909, 70150000: 0.00091, 80160000: 0.24939250000000002, 80170000:
                              9.5e-05, 80180000: 0.0005124999999999999, 180360000: 0.000834, 180380000: 0.00015725, 180400000: 0.24900875}, 1.0, 0.001205, -1.0, {"name": "Air (dry, near sea level)"})
    air_2 = material.Material({10010000: 0.19997700000000002, 10020000: 2.3e-05, 60120000: 0.19786, 60130000: 0.0021400000000000004, 80160000: 0.19951400000000002, 80170000: 7.6e-05,
                              80180000: 0.00040999999999999994, 90190000: 0.2, 140280000: 0.184446, 140290000: 0.00937, 140300000: 0.006184}, 1.0, 1.76, -1.0, {"name": "C-552 Air-Equivalent Plastic"})
    list_of_matches = [
        'C-552 Air-Equivalent Plastic', 'Air (dry, near sea level)']
    material_library = MaterialLibrary()
    material_library[air_1.metadata['name']] = air_1
    material_library[air_2.metadata['name']] = air_2
    assert_equal(
        gtag.print_near_match('air', material_library), list_of_matches)
