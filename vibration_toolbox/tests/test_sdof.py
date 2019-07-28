import numpy as np
import pytest
import vibration_toolbox as vtb
from numpy.testing import assert_allclose


def test_free_response():
    response = np.array([[1.00000000],
                         [0.99591926],
                         [0.99168070],
                         [0.98728508],
                         [0.98273317]])
    assert_allclose(vtb.free_response()[1][:5], response)


def test_fourier_series():
    f = np.hstack((np.arange(-1, 1, .04), np.arange(1, -1, -.04)))
    f += 1
    t = np.arange(0, len(f)) / len(f)
    a, b = vtb.fourier_series(f, t, 5)
    assert_allclose(a[:4],
                    np.array([2.00000000e+00,
                              -8.10836188e-01,
                              7.01998169e-17,
                              -9.03304154e-02]),
                    atol=1e-7,
                    rtol=1e-7)
    assert_allclose(b[:4],
                    np.array([1.13686838e-15,
                              4.45501585e-17,
                              3.37507799e-16,
                              6.09825913e-17]),
                    atol=1e-7,
                    rtol=1e-7)


def test_transmissibility():
    _, D, _ = vtb.transmissibility(zs=[0.05, 0.1, 0.25, 0.5, 0.7],
                                   rmin=0,
                                   rmax=2)
    assert_allclose(D[:10],
                    np.array([1.,
                              1.00001115,
                              1.00004459,
                              1.00010032,
                              1.00017834,
                              1.00027863,
                              1.00040118,
                              1.00054598,
                              1.000713,
                              1.00090222]),
                    rtol=1e-7,
                    atol=1e-7)


def test_rotating_unbalance():
    _, Xn = vtb.rotating_unbalance(m=1, m0=0.5, e=0.1, zs=[0.1], rmin=0,
                                   rmax=3.5, normalized=True)
    assert_allclose(Xn[0, :5],
                    np.array([0. + 0.00000000e+00j,
                              0.01002962 - 2.01187431e-05j,
                              0.02006506 - 8.05225885e-05j,
                              0.03011213 - 1.81354593e-04j,
                              0.04017667 - 3.22853882e-04j]),
                    atol=1e-7)
