"""openmc tests"""
from __future__ import unicode_literals, division
from io import StringIO
import warnings

import nose 
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in, assert_true, assert_false

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)

from pyne import nucname
from pyne import openmc

sample_xs = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
</cross_sections>
""")

def test_ace_table_init():
    atab = openmc.AceTable(zaid='92235', path='U235.ace', 
                           cross_sections_path='/tmp/cross_sections.xml')
    assert_equal('92235', atab.zaid)
    assert_equal('U235.ace', atab.path)
    assert_equal('/tmp/U235.ace', atab.abspath)
    assert_equal(nucname.id('U235'), atab.nucid)

def test_ace_table_xml():
    atab = openmc.AceTable(zaid='92235', path='U235.ace', 
                           cross_sections_path='/tmp/cross_sections.xml')
    exp = '<ace_table path="U235.ace" zaid="92235"/>'
    obs = atab.xml()
    assert_equal(exp, obs)

def test_cross_sections_read():
    sample_xs.seek(0)
    xs = openmc.CrossSections(sample_xs)
    assert_equal('ascii', xs.filetype)
    assert_true(xs.path is None)

    exp = [openmc.AceTable(alias='H-1.71c', awr='0.999167', location='1', 
                           name='1001.71c', path='293.6K/H_001_293.6K.ace',
                           temperature='2.53e-08', zaid='1001'),
           openmc.AceTable(alias='Am-242m.73c', awr='239.9801', location='1',
                           metastable='1', name='95242.73c', 
                           path='900K/Am_242_900K.ace', temperature='7.756e-08',
                           zaid='95242'),
           openmc.AceTable(awr='89.1324', location='1', name='ZrZrH.71t', 
                           path='tsl/zrzrh.acer', temperature='2.551e-08', 
                           zaid='0')]
    assert_equal(exp, xs.ace_tables)

def test_cross_sections_roundtrip():
    sample_xs.seek(0)
    xs = openmc.CrossSections(sample_xs)
    sample_xs.seek(0)
    exp = sample_xs.read()
    obs = xs.xml()
    assert_equal(exp, obs)


if __name__ == "__main__":
    nose.runmodule()

