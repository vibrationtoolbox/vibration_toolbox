import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.linalg as la
import matplotlib as mpl
from scipy.interpolate import UnivariateSpline
from ipywidgets import interact, interactive, widgets
from IPython.display import clear_output, display, HTML
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = (10, 6)
import IPython.core.display as ipcd


from ipywidgets.widgets.interaction import interact, interactive


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


def modes_system_undamped(M, K):
    """
    This function will return the natural frequencies (w),
    eigenvectors (P), mode shapes (S) abd the modal transformation
    matrix S^-1(takes x -> r(modal coordinates) for an undamped system.

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
    lam: array
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
    >>> lam, P, S, Sinv = modes_system_undamped(M, K)
    >>> lam
    array([ 0.19806226+0.j,  1.55495813+0.j,  3.24697960+0.j])
    """
    L = la.cholesky(M)
    Linv = la.inv(L)
    lam, P = eigen(Linv @ K @ Linv.T)
    S = Linv @ P
    Sinv = P.T @ Linv

    return lam, P, S, Sinv


if __name__ == "__main__":
    import doctest
    import vtoolbox as vtb

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
