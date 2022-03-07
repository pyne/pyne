"""CRAM Tests"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from nose.tools import assert_true

from pyne import cram

DTYPES = ["f8", np.complex128]


def test_solve_identity_ones():
    for dtype in DTYPES:
        b = np.ones(cram.N, dtype=dtype)
        mat = sp.eye(cram.N, format="csr", dtype=dtype)
        obs = cram.solve(mat, b)
        exp = spla.spsolve(mat, b)
        yield assert_true, np.allclose(exp, obs)


def test_solve_identity_range():
    for dtype in DTYPES:
        b = np.arange(cram.N, dtype=dtype)
        mat = sp.eye(cram.N, format="csr", dtype=dtype)
        obs = cram.solve(mat, b)
        exp = spla.spsolve(mat, b)
        yield assert_true, np.allclose(exp, obs)


def test_solve_ones_ones():
    for dtype in DTYPES:
        b = np.ones(cram.N, dtype=dtype)
        mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
        obs = cram.solve(mat, b)
        exp = spla.spsolve(mat, b)
        yield assert_true, np.allclose(exp, obs)


def test_solve_ones_range():
    for dtype in DTYPES:
        b = np.arange(cram.N, dtype=dtype)
        mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
        obs = cram.solve(mat, b)
        exp = spla.spsolve(mat, b)
        yield assert_true, np.allclose(exp, obs)


def test_solve_range_range():
    for dtype in DTYPES:
        b = np.arange(cram.N, dtype=dtype)
        mat = cram.ones(dtype=dtype) + sp.diags(
            [b], offsets=[0], shape=(cram.N, cram.N), format="csr", dtype=dtype
        )
        obs = cram.solve(mat, b)
        exp = spla.spsolve(mat, b)
        yield assert_true, np.allclose(exp, obs)


def test_diag_add():
    for dtype in DTYPES:
        mat = cram.ones(dtype=dtype)
        res = mat + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
        exp = cram.flatten_sparse_matrix(res)
        obs = cram.diag_add(mat, 9.0)
        yield assert_true, np.allclose(exp, obs)


def test_dot():
    for dtype in DTYPES:
        x = np.arange(cram.N, dtype=dtype)
        mat = cram.ones(dtype=dtype) + 9 * sp.eye(cram.N, format="csr", dtype=dtype)
        exp = mat.dot(x)
        obs = cram.dot(mat, x)
        yield assert_true, np.allclose(exp, obs)
