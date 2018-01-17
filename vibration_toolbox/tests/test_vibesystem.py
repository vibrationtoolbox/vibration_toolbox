import numpy as np
import pytest
import vibration_toolbox as vtb
from numpy.testing import assert_allclose


@pytest.fixture
def sys0():
    m0, m1 = 1, 1
    c0, c1, c2 = 5, 5, 5
    k0, k1, k2 = 1e3, 1e3, 1e3

    M = np.array([[m0, 0],
                  [0, m1]])
    C = np.array([[c0+c1, -c2],
                  [-c1, c2+c2]])
    K = np.array([[k0+k1, -k2],
                  [-k1, k2+k2]])
    return vtb.VibeSystem(M, C, K)


def test_sys0_freq(sys0):
    assert_allclose(sys0.wn, np.array([31.622777, 54.772256]))
    assert_allclose(sys0.wd, np.array([31.523801, 54.256336]))
