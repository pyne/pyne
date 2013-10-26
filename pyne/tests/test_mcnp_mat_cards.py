from nose.tools import assert_equal
from pyne import mcnp_mat_cards
from pyne.mcnp_mat_cards import get_tag_values
from pyne import mcnp

def test_get_tag_values():
    assert_equal(get_tag_values('Form2.h5m'), ['m_2_rho_0.09', 'tally_0.cell.flux.n', 'm_1_rho_-8.8'])
    assert_equal(get_tag_values('Form3.h5m'), ['M_5_TorusPrism_9', 'imp_4Beryllium', '91_2Steel_3'])
    assert_equal(get_tag_values('Form1.h5m'), ['MAT_SphereCylinder', '31_2_6S_BORON^Steel', '91_2Steel_3', 'tally_<iit>_MAT.He', '.34Lithium', '31_2*6S*Water@Steel', '31_2Boron', 'M_5_TorusPrism_9'])
