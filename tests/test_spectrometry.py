"""Spectrometry tests """
import pytest

import warnings
from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne import spectanalysis as sa
from pyne import gammaspec
gspec1 = gammaspec.read_spe_file('test.spe')
gspec2 = gammaspec.read_dollar_spe_file("gv_format_spect.spe")
gspec3 = gammaspec.read_spec_id_file("gv_format_spect.spe")
gspec4 = gammaspec.read_spec_id_file("gv_format_spect2.spe")
eff_coeff = [-2.818615042612040000, -0.727352820018942000, -0.039579888648190400,
             -0.059230525466409600, 0.023772637347443000, 0.032530647507267100]


def test_read_dollar_spe():
    assert_equal(gspec2.spec_name, "No sample description was entered.")
    assert_equal(gspec2.file_name, "gv_format_spect.spe")
    assert_equal(gspec2.real_time, 209)
    assert_equal(gspec2.live_time, 199)
    assert_equal(gspec2.start_time, "11:43:41")
    assert_equal(gspec2.start_date, "08/01/2014")
    assert_equal(gspec2.dead_time, 10.0)
    assert_equal(gspec2.det_id, "2")
    assert_equal(gspec2.det_descp, "DSPEC1")
    assert_equal(gspec2.start_chan_num, 0)
    assert_equal(gspec2.num_channels, 1024)
    assert_equal(len(gspec2.channels), 1024)
    assert_equal(len(gspec2.counts), 1024)
    assert_equal(len(gspec2.ebin), 1024)
    assert_equal(len(gspec2.energy_channel_fit), 1024)


def test_read_spe():
    assert_equal(gspec1.spec_name, "1K_MIX~1.SPC")
    assert_equal(gspec1.file_name, "test.spe")
    assert_equal(gspec1.real_time, 209.100006)
    assert_equal(gspec1.live_time, 199.800003)
    assert_equal(gspec1.start_time, "11:43:41")
    assert_equal(gspec1.start_date, "01-Aug-2014")
    assert_almost_equal(gspec1.dead_time, 9.300003)
    assert_equal(gspec1.det_id, "2")
    assert_equal(gspec1.det_descp, "DSPEC1")
    assert_equal(gspec1.start_chan_num, 0)
    assert_equal(gspec1.num_channels, 1024)
    assert_equal(len(gspec1.channels), 1024)
    assert_equal(len(gspec1.counts), 1024)
    assert_equal(len(gspec1.ebin), 1024)
    assert_equal(len(gspec1.energy_channel_fit), 1024)
    assert_equal(gspec1.counts[100], gspec2.counts[100])


def test_calib():
    assert gammaspec.calc_e_eff(1, eff_coeff, 1) == 0.059688551591347033
    pytest.raises(ValueError, gammaspec.calc_e_eff, 1, eff_coeff, 10)


def test_str():
    s = str(gspec1)
    assert len(s) > 0


def test_5point_smooth():
    smoothspec = sa.five_point_smooth(gspec1)
    assert smoothspec.start_time == gspec1.start_time
    assert smoothspec.counts[0] == gspec1.counts[0]
    assert smoothspec.counts[-1] == gspec1.counts[-1]
    assert len(smoothspec.counts) == len(gspec1.counts)
    assert gspec1.counts[100] == gspec2.counts[100]


def test_rect_smooth():
    smoothspec = sa.rect_smooth(gspec1, 7)
    assert smoothspec.start_time == gspec1.start_time
    assert smoothspec.counts[0] == gspec1.counts[0]
    assert smoothspec.counts[-1] == gspec1.counts[-1]
    assert len(smoothspec.counts) == len(gspec1.counts)
    assert gspec1.counts[100] == gspec2.counts[100]


def test_calc_bg():
    bg = sa.calc_bg(gspec1, 475, 484, 1)
    pytest.raises(ValueError, sa.calc_bg, gspec1, 500, 484, 1)
    pytest.raises(ValueError, sa.calc_bg, gspec1, -20, 484, 1)
    pytest.raises(ValueError, sa.calc_bg, gspec1, 500, 10484, 1)
    pytest.raises(ValueError, sa.calc_bg, gspec1, 500, 10484, 20)


def test_net_area():
    nc=sa.net_area(gspec1, 475, 484)
    assert_almost_equal(nc, 20454.5)


def test_gross_counts():
    gc = sa.gross_count(gspec1, 475, 484)
    assert gc == sum(gspec1.counts[475:484])
    pytest.raises(ValueError, sa.gross_count, gspec1, 500, 484)
    pytest.raises(ValueError, sa.gross_count, gspec1, -20, 484)
    pytest.raises(ValueError, sa.gross_count, gspec1, 500, 10484)


def test_end_point_average_area():
    nc = sa.end_point_average_area(gspec1,475,484,var=5)
    assert_almost_equal(nc,20355.5) 


def test_fwhm():
    nc = sa.fwhm(gspec4,750,950,3)
    assert_almost_equal(nc, 72.6131488820804)


def test_resolution():
    fwhm = 70.0
    e_0 = 662 
    nc = sa.resolution(fwhm,e_0)
    assert_almost_equal(nc, 0.10574018126888217)
    

def test_read_spec_id_file():
    assert_equal(gspec3.spec_name, "No sample description was entered.")
    assert_equal(gspec3.file_name, "gv_format_spect.spe")
    assert_equal(gspec3.real_time, 209)
    assert_equal(gspec3.live_time, 199)
    assert_equal(gspec3.start_time, "11:43:41")
    assert_equal(gspec3.start_date, "08/01/2014")
    assert_equal(gspec3.dead_time, 10.0)
    assert_equal(gspec3.det_id, "2")
    assert_equal(gspec3.det_descp, "DSPEC1")
    assert_equal(gspec3.start_chan_num, 0)
    assert_equal(gspec3.num_channels, 1024)
    assert_equal(len(gspec3.channels), 1024)
    assert_equal(len(gspec3.counts), 1024)
    assert_equal(len(gspec3.ebin), 1024)
    assert_equal(len(gspec3.energy_channel_fit), 1024)


def test_calc_energy_poly():
    c1,c2,c3 = 860,250,50
    e1,e2,e3 = 662,200,40 
    gspec4.calc_energy_poly(c1,c2,c3,e1,e2,e3)
    assert_almost_equal(gspec4.epoly[2], 0.973600485731636)


def test_net_count():
    nc = sa.net_counts(gspec1, 475, 484, 1)


