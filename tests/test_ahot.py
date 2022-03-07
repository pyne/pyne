import nose
from nose.plugins.skip import Skip, SkipTest

try:
    from spatial_solvers import ahot_script
except:
    raise SkipTest


def test_ahotn_ln():
    ahot_script.test_ahotn_ln()


def test_ahotn_ll():
    ahot_script.test_ahotn_ll()


def test_ahotn_nefd():
    ahot_script.test_ahotn_nefd()


def test_dgfem_ld():
    ahot_script.test_dgfem_ld()


def test_dgfem_dense():
    ahot_script.test_dgfem_dense()


def test_dgfem_lagrange():
    ahot_script.test_dgfem_lagrange()


def test_sct_step():
    ahot_script.test_sct_step()


def test_ahotn_ln_alternating():
    ahot_script.test_ahotn_ln_alternating()


def test_ahotn_ll_alternating():
    ahot_script.test_ahotn_ll_alternating()


def test_ahotn_nefd_alternating():
    ahot_script.test_ahotn_nefd_alternating()


def test_dgfem_ld_alternating():
    ahot_script.test_dgfem_ld_alternating()


def test_dgfem_dense_alternating():
    ahot_script.test_dgfem_dense_alternating()


def test_dgfem_lagrange_alternating():
    ahot_script.test_dgfem_lagrange_alternating()


def test_sct_step_alternating():
    ahot_script.test_sct_step_alternating()


if __name__ == "__main__":
    nose.runmodule()
