#!/usr/bin/python

from pyne import material
from pyne.material import Material, MaterialLibrary
import nose
from nose.tools import assert_equal, assert_raises
import get_tag_values_matlib as gtag


'''test check_matname function'''
"""
Existence/Absence of a graveyard group
"""
def test_garveyard():
    # 'mat:graveyard' group exists
    tag_1=['mat:graveyard','mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_1),[('graveyard', ' '), ('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])
    # 'mat:Graveyard' group exists
    tag_2=['mat:Graveyard','mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_2),[('Graveyard', ' '), ('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])
    # graveyard group is absent
    tag_3=['mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    #assert_raises(Exception, gtag.check_matname, tag_3)
    #graveyard exists as 'graveyard'
    tag_4=['graveyard','mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_4),[('Nitrogen', '-0.001205'), ('Steel, Stainless 321', '-2'), ('Lead', '-11.35'), ('Mercury', '-7.874')])
    
    
    
"""
test groups naming
"""
def test_group_names():
    # ':' is missing 
    tag_5=['mat:Nitrogen/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'matLead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception,gtag.check_matname,tag_5)    
    #material name is absent
    tag_6=['mat:graveyard','mat:/rho:-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception,gtag.check_matname,tag_6)
    #an error in the group name "/" without a density
    tag_7=['mat:graveyard','mat:Nitrogen/', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception,gtag.check_matname,tag_7)
    # ':' is absent in the density part
    tag_8=['mat:graveyard','mat:Nitrogen/rho-0.001205', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321/rho:-2', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_raises(Exception,gtag.check_matname,tag_8)
    #no desnity provided
    tag_9=['mat:graveyard','mat:Nitrogen', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_equal(gtag.check_matname(tag_9),[('Nitrogen', ' '), ('Steel, Stainless 321', ' '), ('Lead', ' '), ('Mercury', ' ')])
    #some densities are provided
    tag_10=['mat:graveyard','mat:Nitrogen', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321', 'mat:Lead/rho:-11.35', 'mat:Mercury/rho:-7.874']
    assert_equal(gtag.check_matname(tag_10),[('Nitrogen', ' '), ('Steel, Stainless 321', ' '), ('Lead', '-11.35'), ('Mercury', '-7.874')]) 
    # no density provided and material name is absent
    tag_11=['mat:graveyard','mat:', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_raises(Exception,gtag.check_matname,tag_11)       
    # no density is provided and ':' is absent
    tag_12=['mat:graveyard','matNitrogen', 'tally_4.cell.flux.p', 'mat:Steel, Stainless 321', 'mat:Lead', 'mat:Mercury']
    assert_raises(Exception,gtag.check_matname,tag_12)
    # 'mat' is absent from the group names
    tag_13=['graveyard','Nitrogen', 'tally_4.cell.flux.p', 'Steel, Stainless 321', 'Lead', 'Mercury']
    assert_raises(Exception,gtag.check_matname,tag_13)


"""
main
"""
def main():    
    test_garveyard()
    test_group_names()
    
    
    
    
if __name__ == '__main__':
    main()        
