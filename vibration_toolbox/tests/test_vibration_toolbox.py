import numpy as np
np.set_printoptions(precision=4, suppress=True)
import vibrationtesting as vt
import numpy.testing as nt

# Examples for future use

'''
M = np.array([[4, 0, 0],
              [0, 4, 0],
              [0, 0, 4]])
K = np.array([[8, -4, 0],
              [-4, 8, -4],
              [0, -4, 4]])
omega, zeta, Psi = vt.sos_modal(M, K)


def test_serep():
    M = np.array([[4, 0, 0],
                  [0, 4, 0],
                  [0, 0, 4]])
    K = np.array([[8, -4, 0],
                  [-4, 8, -4],
                  [0, -4, 4]])
    retained = np.array([[1, 2]])
    Mred, Kred, T, truncated_dofs = vt.serep(M, K, retained)
    Mr_soln = np.array([[16.98791841, -16.19566936],
                        [-16.19566936,  24.19566936]])
    Kr_soln = np.array([[20.98791841, -12.98791841],
                        [-12.98791841,  10.21983253]])
    nt.assert_array_almost_equal(Mred, Mr_soln)
    nt.assert_array_almost_equal(Kred, Kr_soln)

# test_serep()


def test_sos_modal():
    M = np.array([[4, 0, 0],
                  [0, 4, 0],
                  [0, 0, 4]])
    K = np.array([[8, -4, 0],
                  [-4, 8, -4],
                  [0, -4, 4]])
    omega, zeta, Psi = vt.sos_modal(M, K, K / 10)
    nt.assert_array_almost_equal(
        omega, np.array([0.445042,  1.24698,  1.801938]))
    K_diag = np.array([[0.198062,  0., -0.],
                       [0.,  1.554958, -0.],
                       [-0., -0.,  3.24698]])

    nt.assert_array_almost_equal(Psi.T@K@Psi, K_diag)

    K2 = K - np.eye(K.shape[0])@M * (Psi.T@K@Psi)[0, 0]
    omega, zeta, Psi = vt.sos_modal(M, K2)
    nt.assert_array_almost_equal(omega, np.array([0.,  1.164859,  1.746115]))
    Psi_true = np.array([[-0.163993,  0.368488, -0.295505],
                         [-0.295505,  0.163993,  0.368488],
                         [-0.368488, -0.295505, -0.163993]])
    nt.assert_array_almost_equal(Psi_true, Psi)

    # Diagonalizes?
    Psi_true = np.array([[0.,  0., -0.],
                         [-0.,  1.356896,  0.],
                         [-0.,  0.,  3.048917]])
    nt.assert_array_almost_equal(Psi.T@K2@Psi, Psi_true)

    # How about non-proportional damping

    C = K / 10
    C[0, 0] = 2 * C[0, 0]
    omega, zeta, Psi = vt.sos_modal(M, K2, C)

    #  Damping matrix cannot be completely diagonalized.

    nt.assert_array_almost_equal(omega, np.array([0.,  1.164859,  1.746115]))

    nt.assert_array_almost_equal(zeta, np.array([0.,  0.113371,  0.112981]))
    C_diag = np.array([[0.041321, -0.048343,  0.038768],
                       [-0.048343,  0.264123, -0.087112],
                       [0.038768, -0.087112,  0.394556]])
    nt.assert_array_almost_equal(C_diag, Psi.T@C@Psi)
    '''
