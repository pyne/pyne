from nose.tools import assert_equal
import sys
sys.path.insert(1, '/home/zeineddine/pyne/pyne')
from mcnp_mat_cards import get_tag_values
from pyne import mcnp

def test_get_tag_values1():
    assert_equal(get_tag_values('Form14.h5m'), ['m_2_rho_0.09', 'tally_0.cell.flux.n', 'm_1_rho_-8.8'])

def test_get_tag_values2():
    assert_equal(get_tag_values('Form20.h5m'), ['M_5_TorusPrism_9', 'imp_4Beryllium', '91_2Steel_3'])

def test_get_tag_values3():
    assert_equal(get_tag_values('Form10.h5m'), ['MAT_SphereCylinder', '31_2_6S_BORON^Steel', '91_2Steel_3', 'tally_<iit>_MAT.He', '.34Lithium', '31_2*6S*Water@Steel', '31_2Boron', 'M_5_TorusPrism_9'])
