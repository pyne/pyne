"""openmc_utils tests"""
from __future__ import unicode_literals, division
from io import StringIO
import warnings

import nose
from nose.tools import assert_equal, assert_true

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne import nucname
from pyne import openmc_utils

sample_xs = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
</cross_sections>
""")

sample_xs_with_dir = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <directory>/</directory>
  <filetype>ascii</filetype>
  <ace_table alias="H-1.71c" awr="0.999167" location="1" name="1001.71c" path="293.6K/H_001_293.6K.ace" temperature="2.53e-08" zaid="1001"/>
  <ace_table alias="Am-242m.73c" awr="239.9801" location="1" metastable="1" name="95242.73c" path="900K/Am_242_900K.ace" temperature="7.756e-08" zaid="95242"/>
  <ace_table awr="89.1324" location="1" name="ZrZrH.71t" path="tsl/zrzrh.acer" temperature="2.551e-08" zaid="0"/>
""")

sample_xs_with_mcnp_id = StringIO("""<?xml version="1.0" ?>
<cross_sections>
  <filetype>ascii</filetype>
  <ace_table alias="Co-58m.70c" awr="57.4381" location="28897" metastable="1" name="27458.70c" path="endf70b" temperature="2.5301e-08" zaid="27458"/>
  <ace_table alias="Co-58.70c" awr="57.4381" location="28699" name="27058.70c" path="endf70b" temperature="2.5301e-08" zaid="27058"/>
</cross_sections>
""")


def test_ace_table_init():
    atab = openmc_utils.AceTable(zaid='92235', path='U235.ace',
                           cross_sections_path='/tmp/cross_sections.xml')
    assert_equal('92235', atab.zaid)
    assert_equal('U235.ace', atab.path)
    assert_equal('/tmp/U235.ace', atab.abspath)
    assert_equal(nucname.id('U235'), atab.nucid)


def test_ace_table_xml():
    atab = openmc_utils.AceTable(zaid='92235', path='U235.ace',
                           cross_sections_path='/tmp/cross_sections.xml')
    exp = '<ace_table path="U235.ace" zaid="92235"/>'
    obs = atab.xml()
    assert_equal(exp, obs)


def test_cross_sections_read():
    sample_xs.seek(0)
    xs = openmc_utils.CrossSections(sample_xs)
    assert_equal('ascii', xs.filetype)
    assert_true(xs.path is None)

    exp = [openmc_utils.AceTable(alias='H-1.71c', awr='0.999167', location='1',
                           name='1001.71c', path='293.6K/H_001_293.6K.ace',
                           temperature='2.53e-08', zaid='1001'),
           openmc_utils.AceTable(alias='Am-242m.73c', awr='239.9801', location='1',
                           metastable='1', name='95242.73c',
                           path='900K/Am_242_900K.ace', temperature='7.756e-08',
                           zaid='95242'),
           openmc_utils.AceTable(awr='89.1324', location='1', name='ZrZrH.71t',
                           path='tsl/zrzrh.acer', temperature='2.551e-08',
                           zaid='0')]
    assert_equal(exp, xs.ace_tables)


def test_cross_sections_abspath_with_dir():
    xs = openmc_utils.CrossSections(sample_xs_with_dir)
    assert_equal('ascii', xs.filetype)
    assert_equal(xs.path, "/")

    exp_abspaths = ["/293.6K/H_001_293.6K.ace",
                    "/900K/Am_242_900K.ace", "/tsl/zrzrh.acer"]
    obs_abspaths = [table.abspath for table in xs.ace_tables]
    assert_equal(exp_abspaths, obs_abspaths)


def test_cross_sections_mcnp_id():
    xstables = openmc_utils.CrossSections(sample_xs_with_mcnp_id).ace_tables
    mcnp_obs = [table.nucid for table in xstables if table.alias == "Co-58m.70c"][0]
    assert_equal(mcnp_obs, 270580001)
    nucid_obs = [table.nucid for table in xstables if table.alias == "Co-58.70c"][0]
    assert_equal(nucid_obs, 270580000)


def test_cross_sections_roundtrip():
    sample_xs.seek(0)
    xs = openmc_utils.CrossSections(sample_xs)
    sample_xs.seek(0)
    exp = sample_xs.read()
    obs = xs.xml()
    assert_equal(exp, obs)


if __name__ == "__main__":
    nose.runmodule()

