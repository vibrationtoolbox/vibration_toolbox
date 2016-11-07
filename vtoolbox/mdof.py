import scipy as sp
import scipy.linalg as la
import scipy.signal as signal
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = (10, 6)


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
    ----------
    evalues: array
        Sorted eigenvalues
    evectors: array
        Sorted eigenvalues

    Examples:
    >>> L = sp.array([[2, -1, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> lam, P = eigen(L)
    >>> lam
    array([  0.56258062+0.j,   2.63206172+0.j,  10.80535766+0.j])
    """
    if B is None:
        evalues, evectors = la.eig(A)
    else:
        evalues, evectors = la.eig(A, B)

    if all(eigs == 0 for eigs in evalues.imag):
        if all(eigs > 0 for eigs in evalues.real):
            idxp = evalues.real.argsort()  # positive in increasing order
            idxn = sp.array([], dtype=int)
        else:
            idxp = evalues.real.argsort()[int(len(evalues)/2):]  # positive in increasing order
            idxn = evalues.real.argsort()[int(len(evalues)/2) - 1::-1]  # negative in decreasing order

    else:
        idxp = evalues.imag.argsort()[int(len(evalues)/2):]  # positive in increasing order
        idxn = evalues.imag.argsort()[int(len(evalues)/2) - 1::-1]  # negative in decreasing order

    idx = sp.hstack([idxp, idxn])

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
    ----------
    Yn: array
        Normalized matrix

    Examples:
    >>> X = sp.array([[ 0.84+0.j  ,  0.14-0.j  ,  0.84-0.j  ,  0.14+0.j  ],
    ...               [ 0.01-0.3j ,  0.00+0.15j,  0.01+0.3j ,  0.00-0.15j],
    ...               [-0.09+0.42j, -0.01+0.65j, -0.09-0.42j, -0.01-0.65j],
    ...               [ 0.15+0.04j, -0.74+0.j  ,  0.15-0.04j, -0.74-0.j  ]])
    >>> Y = sp.array([[-0.03-0.41j,  0.04+0.1j , -0.03+0.41j,  0.04-0.1j ],
    ...            [ 0.88+0.j  ,  0.68+0.j  ,  0.88-0.j  ,  0.68-0.j  ],
    ...            [-0.21-0.j  ,  0.47+0.05j, -0.21+0.j  ,  0.47-0.05j],
    ...            [ 0.00-0.08j,  0.05-0.54j,  0.00+0.08j,  0.05+0.54j]])
    >>> Yn = normalize(X, Y)
    >>> Yn
    array([[ 0.57822773 -4.69882840e-02j,  0.11696967 -5.85231773e-02j,
             0.57822773 +4.69882840e-02j,  0.11696967 +5.85231773e-02j],
           [ 0.00998912 +1.24180506e+00j, -0.06879320 -8.22911025e-01j,
             0.00998912 -1.24180506e+00j, -0.06879320 +8.22911025e-01j],
           [-0.00238377 -2.96339843e-01j,  0.01295993 -5.73835061e-01j,
            -0.00238377 +2.96339843e-01j,  0.01295993 +5.73835061e-01j],
           [ 0.11289137 -9.08101611e-04j, -0.65854649 -5.87827297e-03j,
             0.11289137 +9.08101611e-04j, -0.65854649 +5.87827297e-03j]])
    """

    Yn = sp.zeros_like(X)
    YTX = Y.T @ X # normalize y so that Y.T @ X will return I
    factors = [1/a for a in sp.diag(YTX)]
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
    ----------
    w: array
        The natural frequencies of the system
    P: array
        The eigenvectors of the system are.
    S: array
        The mode shapes of the system.
    Sinv: array
        The modal transformation matrix S^-1(takes x -> r(modal coordinates))

    Examples:
    >>> M = sp.array([[4, 0, 0],
    ...               [0, 4, 0],
    ...               [0, 0, 4]])
    >>> K = sp.array([[8, -4, 0],
    ...               [-4, 8, -4],
    ...               [0, -4, 4]])
    >>> w, P, S, Sinv = modes_system_undamped(M, K)
    >>> w
    array([ 0.44504187,  1.2469796 ,  1.80193774])
    """
    L = la.cholesky(M)
    Linv = la.inv(L)
    lam, P = eigen(Linv @ K @ Linv.T)
    w = sp.real(sp.sqrt(lam))
    S = Linv @ P
    Sinv = P.T @ Linv

    return w, P, S, Sinv


def modes_system(M, K, C=None):
    """
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
    ----------
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

    Examples:
    >>> M = sp.array([[1, 0],
    ...               [0, 1]])
    >>> K = sp.array([[1, -0.4],
    ...               [0.4, 6]])
    >>> C = sp.array([[0.3, -4],
    ...               [4, 0.2]])
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C)
    Damping is non-proportional, eigenvectors are complex.
    >>> wn
    array([ 0.5206175 ,  4.76729024,  0.5206175 ,  4.76729024])
    >>> wd
    array([ 0.50825856,  4.76531454, -0.50825856, -4.76531454])
    >>> zeta
    array([ 0.21659744,  0.02878692,  0.21659744,  0.02878692])
    >>> X
    array([[ 0.83593587+0.j        ,  0.13536177-0.00208189j,
             0.83593587-0.j        ,  0.13536177+0.00208189j],
           [ 0.00811744-0.29648107j,  0.00444279+0.15426956j,
             0.00811744+0.29648107j,  0.00444279-0.15426956j],
           [-0.09426382+0.42487157j, -0.00865558+0.6453271j ,
            -0.09426382-0.42487157j, -0.00865558-0.6453271j ],
           [ 0.14977369+0.03755828j, -0.73575270+0.j        ,
             0.14977369-0.03755828j, -0.73575270-0.j        ]])
    >>> Y
    array([[ 0.57880640-0.04653746j,  0.12020049-0.05551985j,
             0.57880640+0.04653746j,  0.12020049+0.05551985j],
           [ 0.01204879+1.25620685j, -0.06167340-0.8279964j ,
             0.01204879-1.25620685j, -0.06167340+0.8279964j ],
           [-0.00101305-0.3004544j ,  0.01513120-0.57704216j,
            -0.00101305+0.3004544j ,  0.01513120+0.57704216j],
           [ 0.10657189+0.0025583j , -0.65801243-0.00842571j,
             0.10657189-0.0025583j , -0.65801243+0.00842571j]])
    >>> C = K*2 # with proportional damping
    >>> wn, wd, zeta, X, Y = modes_system(M, K, C)
    Damping is proportional or zero, eigenvectors are real
    >>> X
    array([[-0.99677405,  0.08025891],
           [ 0.08025891, -0.99677405]])
    """

    n = len(M)

    Z = sp.zeros((n, n))
    I = sp.eye(n)
    Minv = la.inv(M)

    if (C is None or sp.all(C == 0) or # check if C has only zero entries
        la.norm(Minv @ C @ K - Minv @ K @ C, 2) < 1e-8*la.norm(Minv @ K @ C, 2)):
        w, P, S, Sinv = modes_system_undamped(M, K)
        wn = w
        wd = w
        zeta = None
        X = P
        Y = P
        print('Damping is proportional or zero, eigenvectors are real')
        return wn, wd, zeta, X, Y

    Z = sp.zeros((n, n))
    I = sp.eye(n)

    # creates the state space matrix
    A = sp.vstack([sp.hstack([Z, I]),
                   sp.hstack([-la.pinv(M) @ K, -la.pinv(M) @ C])])

    w, X = eigen(A)
    _, Y = eigen(A.T)

    wd = sp.imag(w)
    wn = sp.absolute(w)
    zeta = (-sp.real(w)/sp.absolute(w))

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
    ----------
    t: array
        Array with the time
    X: array
        The state-space vector for each time

    Examples:
    >>> M = sp.array([[1, 0],
    ...               [0, 4]])
    >>> K = sp.array([[12, -2],
    ...               [-2, 12]])
    >>> x0 = sp.array([1, 1])
    >>> v0 = sp.array([0, 0])
    >>> max_time = 10
    >>> t, X = response_system_undamped(M, K, x0, v0, max_time)
    >>> X[:, 0] # first column of X will contain the initial conditions [x1, x2, v1, v2]
    array([ 1.,  1.,  0.,  0.])
    >>> X[:, 1] # displacement and velocities after delta t
    array([ 0.99991994,  0.99997998, -0.04001478, -0.01000397])
    """

    t = sp.linspace(0, max_time, int(250 * max_time))
    dt = t[1] - t[0]

    n = len(M)

    Z = sp.zeros((n, n))
    I = sp.eye(n, n)

    # creates the state space matrix
    A = sp.vstack([sp.hstack([Z,               I]),
                   sp.hstack([-la.pinv(M) @ K, Z])])

    # creates the x array and set the first line according to the initial
    # conditions
    X = sp.zeros((2*n, len(t)))
    X[:, 0] = sp.hstack([x0, v0])

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
    ----------
    T : array
        Time values for the output.
    yout : array
        System response.
    xout : array
        Time evolution of the state vector.

    Examples:
    >>> M = sp.array([[9, 0],
    ...               [0, 1]])
    >>> K = sp.array([[27, -3],
    ...               [-3, 3]])
    >>> C = K/10
    >>> x0 = sp.array([0, 1])
    >>> v0 = sp.array([1, 0])
    >>> t = sp.linspace(0, 10, 100)
    >>> F = sp.vstack([0*t,
    ...                3*sp.cos(2*t)])
    >>> tou, yout, xout = response_system(M, C, K, F, x0, v0, t)
    >>> tou[:10]
    array([ 0.        ,  0.1010101 ,  0.2020202 ,  0.3030303 ,  0.4040404 ,
            0.50505051,  0.60606061,  0.70707071,  0.80808081,  0.90909091])
    >>> yout[:10]
    array([[ 0.        ,  1.        ,  1.        ,  0.        ],
           [ 0.1006699 ,  1.0019166 ,  0.9882704 ,  0.04164362],
           [ 0.19867119,  1.00889664,  0.94748001,  0.09763215],
           [ 0.29117815,  1.02159599,  0.87992707,  0.15230092],
           [ 0.37563287,  1.03911711,  0.78859082,  0.19072496],
           [ 0.44980625,  1.05913319,  0.67697613,  0.19968169],
           [ 0.5118431 ,  1.07810368,  0.54895551,  0.16848849],
           [ 0.56029149,  1.09156627,  0.40861572,  0.08967373],
           [ 0.59411703,  1.09448692,  0.26011522, -0.04055038],
           [ 0.61270358,  1.08164667,  0.10755703, -0.22203433]])
    >>> xout[:10]
    array([[ 0.        ,  1.        ,  1.        ,  0.        ],
           [ 0.1006699 ,  1.0019166 ,  0.9882704 ,  0.04164362],
           [ 0.19867119,  1.00889664,  0.94748001,  0.09763215],
           [ 0.29117815,  1.02159599,  0.87992707,  0.15230092],
           [ 0.37563287,  1.03911711,  0.78859082,  0.19072496],
           [ 0.44980625,  1.05913319,  0.67697613,  0.19968169],
           [ 0.5118431 ,  1.07810368,  0.54895551,  0.16848849],
           [ 0.56029149,  1.09156627,  0.40861572,  0.08967373],
           [ 0.59411703,  1.09448692,  0.26011522, -0.04055038],
           [ 0.61270358,  1.08164667,  0.10755703, -0.22203433]])
    """

    n = len(M)

    Z = sp.zeros((n, n))
    I = sp.eye(n)

    # creates the state space matrix
    A = sp.vstack([sp.hstack([Z,               I]),
                   sp.hstack([-la.pinv(M) @ K, -la.pinv(M) @ C])])
    B = sp.vstack([Z,
                   la.inv(M)])
    C = sp.eye(2*n)
    D = 0*B

    sys = signal.lti(A, B, C, D)

    IC = sp.hstack([x0, v0])
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

    Leaving off -v will run the tests without any output. Success will return nothing.

    See the doctest section of the python manual.
    https://docs.python.org/3.5/library/doctest.html
    """
