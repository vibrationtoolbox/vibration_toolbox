import numpy as np
import scipy.linalg as la
import scipy.signal as signal
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = (10, 6)

__all__ = ["modes_system",
           "modes_system_undamped",
           "response_system",
           "response_system_undamped"]


def eigen(A, B=None):
    """
    This function is used to sort eigenvalues and eigenvectors
    e.g. for a given system linalg.eig will return eingenvalues as:
    (array([ 0. +89.4j,  0. -89.4j,  0. +89.4j,  0. -89.4j,  0.+983.2j,  0.-983.2j,  0. +40.7j,  0. -40.7j])
    This function will sort this eigenvalues as:
    (array([ 0. +40.7j,  0. +89.4j,  0. +89.4j,  0.+983.2j,  0. -40.7j,  0. -89.4j,  0. -89.4j,  0.-983.2j])
    Correspondent eigenvectors will follow the same order.

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
    >>> L = np.array([[2, -1, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> lam, P = eigen(L)
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
            idxp = evalues.real.argsort()[int(len(evalues)/2):]  # positive in increasing order
            idxn = evalues.real.argsort()[int(len(evalues)/2) - 1::-1]  # negative in decreasing order

    else:
        idxp = evalues.imag.argsort()[int(len(evalues)/2):]  # positive in increasing order
        idxn = evalues.imag.argsort()[int(len(evalues)/2) - 1::-1]  # negative in decreasing order

    idx = np.hstack([idxp, idxn])

    return evalues[idx], evectors[:, idx]


def normalize(X, Y):
    """
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
    >>> X = np.array([[ 0.84+0.j  ,  0.14-0.j  ,  0.84-0.j  ,  0.14+0.j  ],
    ...               [ 0.01-0.3j ,  0.00+0.15j,  0.01+0.3j ,  0.00-0.15j],
    ...               [-0.09+0.42j, -0.01+0.65j, -0.09-0.42j, -0.01-0.65j],
    ...               [ 0.15+0.04j, -0.74+0.j  ,  0.15-0.04j, -0.74-0.j  ]])
    >>> Y = np.array([[-0.03-0.41j,  0.04+0.1j , -0.03+0.41j,  0.04-0.1j ],
    ...            [ 0.88+0.j  ,  0.68+0.j  ,  0.88-0.j  ,  0.68-0.j  ],
    ...            [-0.21-0.j  ,  0.47+0.05j, -0.21+0.j  ,  0.47-0.05j],
    ...            [ 0.00-0.08j,  0.05-0.54j,  0.00+0.08j,  0.05+0.54j]])
    >>> Yn = normalize(X, Y)
    >>> Yn
    array([[ 0.58-0.05j,  0.12-0.06j,  0.58+0.05j,  0.12+0.06j],
           [ 0.01+1.24j, -0.07-0.82j,  0.01-1.24j, -0.07+0.82j],
           [-0.00-0.3j ,  0.01-0.57j, -0.00+0.3j ,  0.01+0.57j],
           [ 0.11-0.j  , -0.66-0.01j,  0.11+0.j  , -0.66+0.01j]])

    """

    Yn = np.zeros_like(X)
    YTX = Y.T @ X  # normalize y so that Y.T @ X will return I
    factors = [1/a for a in np.diag(YTX)]
    # multiply each column in y by a factor in 'factors'
    for col in enumerate(Y.T):
        Yn[col[0]] = col[1]*factors[col[0]]
    Yn = Yn.T

    return Yn


def modes_system_undamped(M, K):
    """
    This function will return the natural frequencies (w),
    eigenvectors (P), mode shapes (S) abd the modal transformation
    matrix S^-1(takes x -> r(modal coordinates) for an undamped system.

    Parameters
    ----------
    M: array
        Mass matrix
    K: array
        Stiffness matrix

    Returns
    -------
    w: array
        The natural frequencies of the system
    P: array
        The eigenvectors of the system are.
    S: array
        The mode shapes of the system.
    Sinv: array
        The modal transformation matrix S^-1(takes x -> r(modal coordinates))

    Examples
    >>> M = np.array([[4, 0, 0],
    ...               [0, 4, 0],
    ...               [0, 0, 4]])
    >>> K = np.array([[8, -4, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> w, P, S, Sinv = modes_system_undamped(M, K)
    >>> w
    array([ 0.45,  1.25,  1.8 ])

    """

    L = la.cholesky(M)
    Linv = la.inv(L)
    lam, P = eigen(Linv @ K @ Linv.T)
    w = np.real(np.sqrt(lam))
    S = Linv @ P
    Sinv = P.T @ Linv

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
    >>> M = np.array([[1, 0],
    ...               [0, 1]])
    >>> K = np.array([[1, -0.4],
    ...               [0.4, 6]])
    >>> C = np.array([[0.3, -4],
    ...               [4, 0.2]])
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C)
    Damping is non-proportional, eigenvectors are complex.
    >>> wn
    array([ 0.52,  4.77,  0.52,  4.77])
    >>> wd
    array([ 0.51,  4.77, -0.51, -4.77])
    >>> zeta
    array([ 0.22,  0.03,  0.22,  0.03])
    >>> X
    array([[ 0.84+0.j  ,  0.14-0.j  ,  0.84-0.j  ,  0.14+0.j  ],
           [ 0.01-0.3j ,  0.00+0.15j,  0.01+0.3j ,  0.00-0.15j],
           [-0.09+0.42j, -0.01+0.65j, -0.09-0.42j, -0.01-0.65j],
           [ 0.15+0.04j, -0.74+0.j  ,  0.15-0.04j, -0.74-0.j  ]])
    >>> Y
    array([[ 0.58-0.05j,  0.12-0.06j,  0.58+0.05j,  0.12+0.06j],
           [ 0.01+1.26j, -0.06-0.83j,  0.01-1.26j, -0.06+0.83j],
           [-0.00-0.3j ,  0.02-0.58j, -0.00+0.3j ,  0.02+0.58j],
           [ 0.11+0.j  , -0.66-0.01j,  0.11-0.j  , -0.66+0.01j]])
    >>> C = K*2 # with proportional damping
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C)
    Damping is proportional or zero, eigenvectors are real
    >>> X
    array([[-1.  ,  0.08],
           [ 0.08, -1.  ]])
    """

    n = len(M)

    Z = np.zeros((n, n))
    I = np.eye(n)
    Minv = la.inv(M)

    if (C is None or np.all(C == 0) or  # check if C has only zero entries
            la.norm(Minv @ C @ K - Minv @ K @ C, 2) <
            1e-8*la.norm(Minv @ K @ C, 2)):
        w, P, S, Sinv = modes_system_undamped(M, K)
        wn = w
        wd = w
        zeta = None
        X = P
        Y = P
        print('Damping is proportional or zero, eigenvectors are real')
        return wn, wd, zeta, X, Y

    Z = np.zeros((n, n))
    I = np.eye(n)

    # creates the state space matrix
    A = np.vstack([np.hstack([Z, I]),
                   np.hstack([-la.pinv(M) @ K, -la.pinv(M) @ C])])

    w, X = eigen(A)
    _, Y = eigen(A.T)

    wd = np.imag(w)
    wn = np.absolute(w)
    zeta = (-np.real(w)/np.absolute(w))

    Y = normalize(X, Y)

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
    >>> X[:, 0] # first column is the initial conditions [x1, x2, v1, v2]
    array([ 1.,  1.,  0.,  0.])
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
                   np.hstack([-la.pinv(M) @ K, Z])])

    # creates the x array and set the first line according to the initial
    # conditions
    X = np.zeros((2*n, len(t)))
    X[:, 0] = np.hstack([x0, v0])

    Ad = la.expm(A * dt)
    for i in range(len(t) - 1):
        X[:, i + 1] = Ad @ X[:, i]

    return t, X


def response_system(M, C, K, F, x0, v0, t):
    """
    This function solves the system given the initial
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
    >>> tou[:10]
    array([ 0.  ,  0.1 ,  0.2 ,  0.3 ,  0.4 ,  0.51,  0.61,  0.71,  0.81,  0.91])

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
                   np.hstack([-la.pinv(M) @ K, -la.pinv(M) @ C])])
    B = np.vstack([Z,
                   la.inv(M)])
    C = np.eye(2*n)
    D = 0*B

    sys = signal.lti(A, B, C, D)

    IC = np.hstack([x0, v0])
    F = F.T
    T, yout, xout = signal.lsim(sys, F, t, IC)

    return T, yout, xout


if __name__ == "__main__":
    import doctest

    doctest.testmod(optionflags=doctest.ELLIPSIS)
    # doctest.run_docstring_examples(frfest,globals(),optionflags=doctest.ELLIPSIS)
    # doctest.run_docstring_examples(asd,globals(),optionflags=doctest.ELLIPSIS)
    """ What this does.
    python (name of this file)  -v
    will test all of the examples in the help.

    Leaving off -v will run the tests without any output. Success will return
    nothing.

    See the doctest section of the python manual.
    https://docs.python.org/3.5/library/doctest.html
    """
