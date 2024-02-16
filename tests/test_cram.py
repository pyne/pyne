"""CRAM Tests"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import pytest

from pyne import cram

@pytest.fixture(params=["f8", np.complex128])
def dtype(request):
    return request.param


def test_solve_identity_ones(dtype):
   b = np.ones(cram.N, dtype=dtype)
   mat = sp.eye(cram.N, format="csr", dtype=dtype)
   obs = cram.solve(mat, b)
   exp = spla.spsolve(mat, b)
   assert np.allclose(exp, obs)


def test_solve_identity_range(dtype):
    b = np.arange(cram.N, dtype=dtype)
    mat = sp.eye(cram.N, format="csr", dtype=dtype)
    obs = cram.solve(mat, b)
    exp = spla.spsolve(mat, b)
    assert np.allclose(exp, obs)


def test_solve_ones_ones(dtype):
    b = np.ones(cram.N, dtype=dtype)
    mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
    obs = cram.solve(mat, b)
    exp = spla.spsolve(mat, b)
    assert np.allclose(exp, obs)


def test_solve_ones_range(dtype):
    b = np.arange(cram.N, dtype=dtype)
    mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
    obs = cram.solve(mat, b)
    exp = spla.spsolve(mat, b)
    assert np.allclose(exp, obs)


def test_solve_range_range(dtype):
    b = np.arange(cram.N, dtype=dtype)
    mat = cram.ones(dtype=dtype) + sp.diags(
        [b], offsets=[0], shape=(cram.N, cram.N), format="csr", dtype=dtype
    )
    obs = cram.solve(mat, b)
    exp = spla.spsolve(mat, b)
    assert np.allclose(exp, obs)


def test_diag_add(dtype):
    mat = cram.ones(dtype=dtype)
    res = mat + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
    exp = cram.flatten_sparse_matrix(res)
    obs = cram.diag_add(mat, 9.0)
    assert np.allclose(exp, obs)


def test_dot(dtype):
    x = np.arange(cram.N, dtype=dtype)
    mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
    exp = mat.dot(x)
    obs = cram.dot(mat, x)
    assert np.allclose(exp, obs)
