"""Multiple Degree of Freedom Analysis Tools."""

import numpy as np
import scipy.linalg as la
import scipy.signal as signal
import matplotlib as mpl

__all__ = ["modes_system",
           "modes_system_undamped",
           "response_system",
           "response_system_undamped"]


mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = (10, 6)


def _eigen(A, B=None):
    """Return sorted eigenvector/eigenvalue pairs.

    e.g. for a given system linalg.eig will return eingenvalues as:
    (array([ 0. +89.4j,  0. -89.4j,  0. +89.4j,  0. -89.4j,  0.+983.2j,
             0.-983.2j,  0. +40.7j,  0. -40.7j])
    This function will sort this eigenvalues as:
    (array([ 0. +40.7j,  0. +89.4j,  0. +89.4j,  0.+983.2j,  0. -40.7j,
             0. -89.4j,  0. -89.4j,  0.-983.2j])
    Correspondent eigenvectors will follow the same order.

    Note: Works fine for moderately sized models. Does not leverage the
    full set of constraints to optimize the solution. See the vibrationtesting
    module for a more advanced solver.

    Parameters
    ----------
    A: array
        A complex or real matrix whose eigenvalues and eigenvectors
        will be computed.
    B: float or str
        Right-hand side matrix in a generalized eigenvalue problem.
        Default is None, identity matrix is assumed.

    Returns
    -------
    evalues: array
        Sorted eigenvalues
    evectors: array
        Sorted eigenvalues

    Examples
    --------
    >>> import vibration_toolbox as vtb
    >>> L = np.array([[2, -1, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> lam, P = vtb.mdof._eigen(L)
    >>> lam
    array([  0.56+0.j,   2.63+0.j,  10.81+0.j])

    """
    if B is None:
        evalues, evectors = la.eig(A)
    else:
        evalues, evectors = la.eig(A, B)

    if all(eigs == 0 for eigs in evalues.imag):
        if all(eigs > 0 for eigs in evalues.real):
            idxp = evalues.real.argsort()  # positive in increasing order
            idxn = np.array([], dtype=int)
        else:
            # positive in increasing order
            idxp = evalues.real.argsort()[int(len(evalues) / 2):]
            # negative in decreasing order
            idxn = evalues.real.argsort()[int(len(evalues) / 2) - 1::-1]

    else:
        # positive in increasing order
        idxp = evalues.imag.argsort()[int(len(evalues) / 2):]
        # negative in decreasing order
        idxn = evalues.imag.argsort()[int(len(evalues) / 2) - 1::-1]

    idx = np.hstack([idxp, idxn])

    return evalues[idx], evectors[:, idx]


def _normalize(X, Y):
    """
    Return normalized left eigenvectors.

    This function is used to normalize vectors of the matrix
    Y with respect to X so that Y.T @ X = I (identity).
    This is used to normalize the matrix with the left eigenvectors.

    Parameters
    ----------
    X: array
        A complex or real matrix
    Y: array
        A complex or real matrix to be normalized

    Returns
    -------
    Yn: array
        Normalized matrix

    Examples
    --------
    >>> # This test has been moved to tests_mdof
    >>> X = np.array([[ 0.84+0.j  ,  0.14-0.j  ,  0.84-0.j  ,  0.14+0.j  ],
    ...               [ 0.01-0.3j ,  0.00+0.15j,  0.01+0.3j ,  0.00-0.15j],
    ...               [-0.09+0.42j, -0.01+0.65j, -0.09-0.42j, -0.01-0.65j],
    ...               [ 0.15+0.04j, -0.74+0.j  ,  0.15-0.04j, -0.74-0.j  ]])
    >>> Y = np.array([[-0.03-0.41j,  0.04+0.1j , -0.03+0.41j,  0.04-0.1j ],
    ...            [ 0.88+0.j  ,  0.68+0.j  ,  0.88-0.j  ,  0.68-0.j  ],
    ...            [-0.21-0.j  ,  0.47+0.05j, -0.21+0.j  ,  0.47-0.05j],
    ...            [ 0.00-0.08j,  0.05-0.54j,  0.00+0.08j,  0.05+0.54j]])
    >>> Yn = _normalize(X, Y)
    >>> Yn # doctest: +SKIP
    array([[ 0.58-0.05j,  0.12-0.06j,  0.58+0.05j,  0.12+0.06j],
           [ 0.01+1.24j, -0.07-0.82j,  0.01-1.24j, -0.07+0.82j],
           [-0.  -0.3j ,  0.01-0.57j, -0.  +0.3j ,  0.01+0.57j],
           [ 0.11-0.j  , -0.66-0.01j,  0.11+0.j  , -0.66+0.01j]])

    """
    Yn = np.zeros_like(X)
    YTX = Y.T @ X  # normalize y so that Y.T @ X will return I
    factors = [1 / a for a in np.diag(YTX)]
    # multiply each column in y by a factor in 'factors'
    for col in enumerate(Y.T):
        Yn[col[0]] = col[1] * factors[col[0]]
    Yn = Yn.T

    return Yn


def modes_system_undamped(M, K):
    r"""Return eigensolution of multiple DOF system.

    Returns the natural frequencies (w),
    eigenvectors (P), mode shapes (S) and the modal transformation
    matrix S for an undamped system.

    See Notes for explanation of the underlying math.

    Parameters
    ----------
    M: float array
        Mass matrix
    K: float array
        Stiffness matrix

    Returns
    -------
    w: float array
        The natural frequencies of the system
    P: float array
        The eigenvectors of the system.
    S: float array
        The mass-normalized mode shapes of the system.
    Sinv: float array
        The modal transformation matrix S^-1(takes x -> r(modal coordinates))

    Notes
    -----
    Given :math:`M\ddot{x}(t)+Kx(t)=0`, with mode shapes :math:`u`, the matrix
    of mode shapes :math:`S=[u_1 u_1 \ldots]` can be created. If the modal
    coordinates are the vector :math:`r(t)`. The modal transformation separates
    space and time from :math:`x(t)` such that :math:`x(t)=S r(t)`.
    Substituting into the governing equation:

    :math:`MS\ddot{r}(t)+KSr(t)=0`

    Premultiplying by :math:`S^T`

    :math:`S^TMS\ddot{r}(t)+S^TKSr(t)=0`

    The matrices :math:`S^TMS` and :math:`S^TKS` will be diagonalized by this
    process (:math:`u_i` are the eigenvectors of :math:`M^{-1}K`).

    If scaled properly (mass normalized so :math:`u_i^TMu_i=1`) then
    :math:`S^TMS=I` and :math:`S^TKS=\Omega^2` where :math:`\Omega^2` is a
    diagonal matrix of the natural frequencies squared in radians per second.

    Further, inverses are unstable so the better way to solve linear equations is with
    Gauss elimination.

    :math:`AB=C` given known :math:`A` and :math:`C`
    is solved using `la.solve(A, C, assume_a='pos')`.

    :math:`BA=C` given known :math:`A` and :math:`C` is solved by first
    transposing the equation to :math:`A^TB^T=C^T`, then solving for
    :math:`C^T`. The resulting command is
    `la.solve(A.T, C.T, assume_a='pos').T`

    Examples
    --------
    >>> M = np.array([[4, 0, 0],
    ...               [0, 4, 0],
    ...               [0, 0, 4]])
    >>> K = np.array([[8, -4, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> w, P, S, Sinv = modes_system_undamped(M, K)
    >>> w # doctest: +SKIP
    array([0.45, 1.25, 1.8 ])
    >>> S
    array([[ 0.16, -0.37, -0.3 ],
           [ 0.3 , -0.16,  0.37],
           [ 0.37,  0.3 , -0.16]])

    """
    L = la.cholesky(M)
    lam, P = _eigen(la.solve(L, la.solve(L, K, assume_a='pos').T,
                             assume_a='pos').T)
    w = np.real(np.sqrt(lam))
    S = la.solve(L, P, assume_a='pos')
    Sinv = la.solve(L.T, P, assume_a='pos').T

    return w, P, S, Sinv


def modes_system(M, K, C=None):
    """Natural frequencies, damping ratios, and mode shapes of MDOF system.
    This function will return the natural frequencies (wn), the
    damped natural frequencies (wd), the damping ratios (zeta),
    the right eigenvectors (X) and the left eigenvectors (Y) for a
    system defined by M, K and C.
    If the dampind matrix 'C' is none or if the damping is proportional,
    wd and zeta will be none and X and Y will be equal.

    Parameters
    ----------
    M: array
        Mass matrix
    K: array
        Stiffness matrix
    C: array
        Damping matrix

    Returns
    -------
    wn: array
        The natural frequencies of the system
    wd: array
        The damped natural frequencies of the system
    zeta: array
        The damping ratios
    X: array
        The right eigenvectors
    Y: array
        The left eigenvectors

    Examples
    --------
    >>> # This test has been moved to tests_mdof
    >>> M = np.array([[1, 0],
    ...               [0, 1]])
    >>> K = np.array([[2, -1],
    ...               [-1, 6]])
    >>> C = np.array([[0.3, -0.02],
    ...               [-0.02, 0.1]])
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C) # doctest: +SKIP
    Damping is non-proportional, eigenvectors are complex.
    >>> wn # doctest: +SKIP
    array([1.33, 2.5 , 1.33, 2.5 ])
    >>> wd # doctest: +SKIP
    array([1.32, 2.5 , 1.32, 2.5 ])
    >>> zeta # doctest: +SKIP
    array([0.11, 0.02, 0.11, 0.02])
    >>> X # doctest: +SKIP
    array([[-0.06-0.58j, -0.01+0.08j, -0.06+0.58j, -0.01-0.08j],
           [-0.  -0.14j, -0.01-0.36j, -0.  +0.14j, -0.01+0.36j],
           [ 0.78+0.j  , -0.21-0.03j,  0.78-0.j  , -0.21+0.03j],
           [ 0.18+0.01j,  0.9 +0.j  ,  0.18-0.01j,  0.9 -0.j  ]])
    >>> Y # doctest: +SKIP
    array([[ 0.02+0.82j,  0.01-0.31j,  0.02-0.82j,  0.01+0.31j],
           [-0.05+0.18j,  0.01+1.31j, -0.05-0.18j,  0.01-1.31j],
           [ 0.61+0.06j, -0.12-0.02j,  0.61-0.06j, -0.12+0.02j],
           [ 0.14+0.03j,  0.53+0.j  ,  0.14-0.03j,  0.53-0.j  ]])
    >>> C = 0.2*K # with proportional damping
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C) # doctest: +SKIP
    Damping is proportional or zero, eigenvectors are real
    >>> X # doctest: +SKIP
    array([[-0.97,  0.23],
           [-0.23, -0.97]])
    """

    n = len(M)

    Z = np.zeros((n, n))
    I = np.eye(n)

    if (C is None or np.all(C == 0) or  # check if C has only zero entries
            la.norm(la.solve(M, C, assume_a='pos') @ K
                    - la.solve(M, K, assume_a='pos') @ C, 2) <
            1e-8 * la.norm(la.solve(M, K, assume_a='pos') @ C, 2)):
        w, P, S, Sinv = modes_system_undamped(M, K)
        wn = w
        wd = w
        #zeta = None
        zeta = np.diag(S.T@C@S) / 2 / wn
        wd = wn * np.sqrt(1 - zeta**2)
        X = P
        Y = P
        print('Damping is proportional or zero, eigenvectors are real')
        return wn, wd, zeta, X, Y

    Z = np.zeros((n, n))
    I = np.eye(n)

    # creates the state space matrix
    A = np.vstack([np.hstack([Z, I]),
                   np.hstack([-la.solve(M, K, assume_a='pos'),
                              -la.solve(M, C, assume_a='pos')])])

    w, X = _eigen(A)
    _, Y = _eigen(A.T)

    wd = abs(np.imag(w))
    wn = np.absolute(w)
    zeta = (-np.real(w) / np.absolute(w))

    Y = _normalize(X, Y)

    print('Damping is non-proportional, eigenvectors are complex.')

    return wn, wd, zeta, X, Y


def response_system_undamped(M, K, x0, v0, max_time):
    """
    This function calculates the time response for an undamped system
    and returns the vector (state-space) X. The n first rows contain the
    displacement (x) and the n last rows contain velocity (v) for each
    coordinate. Each column is related to a time-step.
    The time array is also returned.

    Parameters
    ----------
    M: array
        Mass matrix
    K: array
        Stiffness matrix
    x0: array
        Array with displacement initial conditions
    v0: array
        Array with velocity initial conditions
    max_time: float
        End time

    Returns
    -------
    t: array
        Array with the time
    X: array
        The state-space vector for each time

    Examples
    --------
    >>> M = np.array([[1, 0],
    ...               [0, 4]])
    >>> K = np.array([[12, -2],
    ...               [-2, 12]])
    >>> x0 = np.array([1, 1])
    >>> v0 = np.array([0, 0])
    >>> max_time = 10
    >>> t, X = response_system_undamped(M, K, x0, v0, max_time)
    >>> # first column is the initial conditions [x1, x2, v1, v2]
    >>> X[:, 0] # doctest: +SKIP
    array([1., 1., 0., 0.])
    >>> X[:, 1] # displacement and velocities after delta t
    array([ 1.  ,  1.  , -0.04, -0.01])
    """

    t = np.linspace(0, max_time, int(250 * max_time))
    dt = t[1] - t[0]

    n = len(M)

    Z = np.zeros((n, n))
    I = np.eye(n, n)

    # creates the state space matrix
    A = np.vstack([np.hstack([Z,               I]),
                   np.hstack([-la.solve(M, K, assume_a='pos'), Z])])

    # creates the x array and set the first line according to the initial
    # conditions
    X = np.zeros((2 * n, len(t)))
    X[:, 0] = np.hstack([x0, v0])

    Ad = la.expm(A * dt)
    for i in range(len(t) - 1):
        X[:, i + 1] = Ad @ X[:, i]

    return t, X


def response_system(M, C, K, F, x0, v0, t):
    """
    Returns system response given the initial
    displacement vector 'X0', initial velocity vector 'V0',
    the mass matrix 'M', the stiffness matrix 'M', and the damping
    matrix 'C' and force 'F'.
    T is a row vector of evenly spaced times.
    F is a matrix of forces over time, each column corresponding
    to the corresponding column of T, each row corresponding to
    the same numbered DOF.

    Parameters
    ----------
    M: array
        Mass matrix
    K: array
        Stiffness matrix
    C: array
        Damping matrix
    x0: array
        Array with displacement initial conditions
    v0: array
        Array with velocity initial conditions
    t: array
        Array withe evenly spaced times

    Returns
    -------
    T : array
        Time values for the output.
    yout : array
        System response.
    xout : array
        Time evolution of the state vector.

    Examples
    --------
    >>> M = np.array([[9, 0],
    ...               [0, 1]])
    >>> K = np.array([[27, -3],
    ...               [-3, 3]])
    >>> C = K/10
    >>> x0 = np.array([0, 1])
    >>> v0 = np.array([1, 0])
    >>> t = np.linspace(0, 10, 100)
    >>> F = np.vstack([0*t,
    ...                3*np.cos(2*t)])
    >>> tou, yout, xout = response_system(M, C, K, F, x0, v0, t)
    >>> tou[:10] # doctest: +SKIP
    array([0.  , 0.1 , 0.2 , 0.3 , 0.4 , 0.51, 0.61, 0.71, 0.81, 0.91])

    >>> yout[:10]
    array([[ 0.  ,  1.  ,  1.  ,  0.  ],
           [ 0.1 ,  1.  ,  0.99,  0.04],
           [ 0.2 ,  1.01,  0.95,  0.1 ],
           [ 0.29,  1.02,  0.88,  0.15],
           [ 0.38,  1.04,  0.79,  0.19],
           [ 0.45,  1.06,  0.68,  0.2 ],
           [ 0.51,  1.08,  0.55,  0.17],
           [ 0.56,  1.09,  0.41,  0.09],
           [ 0.59,  1.09,  0.26, -0.04],
           [ 0.61,  1.08,  0.11, -0.22]])

    >>> xout[:10]
    array([[ 0.  ,  1.  ,  1.  ,  0.  ],
           [ 0.1 ,  1.  ,  0.99,  0.04],
           [ 0.2 ,  1.01,  0.95,  0.1 ],
           [ 0.29,  1.02,  0.88,  0.15],
           [ 0.38,  1.04,  0.79,  0.19],
           [ 0.45,  1.06,  0.68,  0.2 ],
           [ 0.51,  1.08,  0.55,  0.17],
           [ 0.56,  1.09,  0.41,  0.09],
           [ 0.59,  1.09,  0.26, -0.04],
           [ 0.61,  1.08,  0.11, -0.22]])
    """

    n = len(M)

    Z = np.zeros((n, n))
    I = np.eye(n)

    # creates the state space matrix
    A = np.vstack([np.hstack([Z,               I]),
                   np.hstack([-la.solve(M, K, assume_a='pos'),
                              -la.solve(M, C, assume_a='pos')])])
    B = np.vstack([Z,
                   la.inv(M)])
    C = np.eye(2 * n)
    D = 0 * B

    sys = signal.lti(A, B, C, D)

    IC = np.hstack([x0, v0])
    F = F.T
    T, yout, xout = signal.lsim(sys, F, t, IC)

    return T, yout, xout
