import numpy as np
import pytest
import vibration_toolbox as vtb
from numpy.testing import assert_allclose


@pytest.fixture
def _normalize_vect():
    X = np.array([[0.84 + 0.j,  0.14 - 0.j,  0.84 - 0.j,  0.14 + 0.j],
                  [0.01 - 0.3j,  0.00 + 0.15j,  0.01 + 0.3j,  0.00 - 0.15j],
                  [-0.09 + 0.42j, -0.01 + 0.65j, -0.09 - 0.42j, -0.01 - 0.65j],
                  [0.15 + 0.04j, -0.74 + 0.j,  0.15 - 0.04j, -0.74 - 0.j]])
    Y = np.array([[-0.03 - 0.41j,  0.04 + 0.1j, -0.03 + 0.41j,  0.04 - 0.1j],
                  [0.88 + 0.j,  0.68 + 0.j,  0.88 - 0.j,  0.68 - 0.j],
                  [-0.21 - 0.j,  0.47 + 0.05j, -0.21 + 0.j,  0.47 - 0.05j],
                  [0.00 - 0.08j,  0.05 - 0.54j,  0.00 + 0.08j,  0.05 + 0.54j]])
    Yn = vtb.mdof._normalize(X, Y)
    print(Yn)
    return Yn


def test_normalise(_normalize_vect):
    Yn = np.array([[0.57822773 - 4.69882840e-02j,  0.11696967 - 5.85231773e-02j,
                    0.57822773 + 4.69882840e-02j,  0.11696967 + 5.85231773e-02j],
                   [0.00998912 + 1.24180506e+00j, -0.0687932 - 8.22911025e-01j,
                    0.00998912 - 1.24180506e+00j, -0.0687932 + 8.22911025e-01j],
                   [-0.00238377 - 2.96339843e-01j,  0.01295993 - 5.73835061e-01j,
                    -0.00238377 + 2.96339843e-01j,  0.01295993 + 5.73835061e-01j],
                   [0.11289137 - 9.08101611e-04j, -0.65854649 - 5.87827297e-03j,
                    0.11289137 + 9.08101611e-04j, -0.65854649 + 5.87827297e-03j]])
    assert_allclose(_normalize_vect, Yn)


def modes_sys():
    M = np.array([[1, 0],                [0, 1]])
    K = np.array([[2, -1], [-1, 6]])
    C = np.array([[0.3, -0.02], [-0.02, 0.1]])
    wn, wd, zeta, X, Y = vtb.modes_system(M, K, C)
    return wn, wd, zeta, X, Y


def test_modes_sys_np():
    wn = np.array([1.328707, 2.49613, 1.328707, 2.49613])
    wd = np.array([1.321268, 2.495419, 1.321268, 2.495419])
    zeta = np.array([0.105665, 0.023878, 0.105665, 0.023878])
    X = np.array([[-0.06187175 - 0.58226568j, -0.01117454 + 0.08415145j,
                   -0.06187175 + 0.58226568j,
                   -0.01117454 - 0.08415145j],
                  [-0.0031966 - 0.13686382j, -0.00864532 - 0.36196512j,
                   -0.0031966 + 0.13686382j,
                   -0.00864532 + 0.36196512j],
                  [0.77801575 + 0.0j,  -0.20932709 - 0.03290071j,
                   0.77801575 - 0.0j,
                   -0.20932709 + 0.03290071j],
                  [0.1812826 + 0.0149919j, 0.90376982 + 0.0j,
                   0.1812826 - 0.0149919j,
                   0.90376982 - 0.0j]])
    Y = np.array([[0.01582041 + 0.81605587j, 0.01036355 - 0.30817433j,
                   0.01582041 - 0.81605587j,
                   0.01036355 + 0.30817433j],
                  [-0.05184151 + 0.18429687j, 0.01345577 + 1.31152703j,
                   -0.05184151 - 0.18429687j,
                   0.01345577 - 1.31152703j],
                  [0.61081141 + 0.05967192j, -0.12152971 - 0.02007119j,
                   0.61081141 - 0.05967192j,
                   -0.12152971 + 0.02007119j],
                  [0.14117309 + 0.02567391j, 0.5253469 + 0.00408669j,
                   0.14117309 - 0.02567391j,
                   0.5253469 - 0.00408669j]])

    wn_func, wd_func, zeta_func, X_func, Y_func = modes_sys()

    assert_allclose(wn_func, wn, rtol=1e-05)
    assert_allclose(wd_func, wd, rtol=1e-05)
    assert_allclose(zeta_func, zeta, rtol=1e-04)
    np.set_printoptions(precision=8)
    print(Y_func)
    assert_allclose(X_func, X, rtol=1e-05)
    assert_allclose(Y_func, Y, rtol=1e-05)


def modes_sys_prop():
    M = np.array([[1, 0],                [0, 1]])
    K = np.array([[2, -1], [-1, 6]])
    C = 0.2 * K
    wn, wd, zeta, X, Y = vtb.modes_system(M, K, C)
    return (wn, wd, zeta, X, Y)


def test_modes_sys_prop():
    wn, wd, zeta, X, Y = modes_sys_prop()
    assert_allclose(X, np.array([[-0.973249,  0.229753],
                                 [-0.229753, -0.973249]]), rtol=1e-05)
