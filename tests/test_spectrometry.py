"""Spectrometry tests """
import pytest

import numpy as np
import warnings
from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)

from pyne import gammaspec
from pyne import spectanalysis as sa

gspec1 = gammaspec.read_spe_file("test.spe")
gspec2 = gammaspec.read_dollar_spe_file("gv_format_spect.spe")
eff_coeff = [
    -2.818615042612040000,
    -0.727352820018942000,
    -0.039579888648190400,
    -0.059230525466409600,
    0.023772637347443000,
    0.032530647507267100,
]


def test_read_dollar_spe():
    assert gspec2.spec_name == "No sample description was entered."
    assert gspec2.file_name == "gv_format_spect.spe"
    assert gspec2.real_time == 209
    assert gspec2.live_time == 199
    assert gspec2.start_time == "11:43:41"
    assert gspec2.start_date == "08/01/2014"
    assert gspec2.dead_time == 10.0
    assert gspec2.det_id == "2"
    assert gspec2.det_descp == "DSPEC1"
    assert gspec2.start_chan_num == 0
    assert gspec2.num_channels == 1024
    assert len(gspec2.channels) == 1024
    assert len(gspec2.counts) == 1024
    assert len(gspec2.ebin) == 1024


def test_read_spe():
    assert gspec1.spec_name == "1K_MIX~1.SPC"
    assert gspec1.file_name == "test.spe"
    assert gspec1.real_time == 209.100006
    assert gspec1.live_time == 199.800003
    assert gspec1.start_time == "11:43:41"
    assert gspec1.start_date == "01-Aug-2014"
    assert gspec1.dead_time == pytest.approx(9.300003)
    assert gspec1.det_id == "2"
    assert gspec1.det_descp == "DSPEC1"
    assert gspec1.start_chan_num == 0
    assert gspec1.num_channels == 1024
    assert len(gspec1.channels) == 1024
    assert len(gspec1.counts) == 1024
    assert len(gspec1.ebin) == 1024
    assert gspec1.counts[100] == gspec2.counts[100]


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


def test_gross_counts():
    gc = sa.gross_count(gspec1, 475, 484)
    assert gc == sum(gspec1.counts[475:484])
    pytest.raises(ValueError, sa.gross_count, gspec1, 500, 484)
    pytest.raises(ValueError, sa.gross_count, gspec1, -20, 484)
    pytest.raises(ValueError, sa.gross_count, gspec1, 500, 10484)


def test_net_count():
    nc = sa.net_counts(gspec1, 475, 484, 1)

