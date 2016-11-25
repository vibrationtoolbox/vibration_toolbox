import numpy as np
import scipy.linalg as la
from scipy import signal

__all__ = ['VibSystem']


class VibSystem(object):
    def __init__(self, M, C, K, name=None):
        r"""A multiple degrees of freedom system.

        This class will create a multiple degree of
        freedom system given M, C and K matrices.

        Parameters
        ----------
        M : array
            Mass matrix.
        C : array
            Damping matrix.
        K : array
            Stiffness matrix.
        name : str
            Name of the system.

        Attributes
        ----------
        evalues : array
            System's eigenvalues.
        evectors : array
            System's eigenvectors.
        wn : array
            System's natural frequencies in Hz.
        wd : array
            System's damped natural frequencies in Hz.
        H:
            Continuous-time linear time invariant system

        Examples
        --------
        For a system consisting of two masses connected
        to each other and both connected to a wall we have
        the following matrices:

        >>> m1, m2 = 1, 1
        >>> c1, c2, c3 = 1, 1, 1
        >>> k1, k2, k3 = 1e3, 1e3, 1e3

        >>> M = np.array([[m1, 0],
        ...               [0, m2]])
        >>> C = np.array([[c1+c2, -c2],
        ...               [-c2, c2+c3]])
        >>> K = np.array([[k1+k2, -k2],
        ...               [-k2, k2+k3]])
        >>> sys1 = VibSystem(M, C, K)
        >>> sys1.wn
        array([ 5.03292121,  8.71727525])
        >>> sys1.wd
        array([ 5.03229206,  8.71400566])
        """
        self._M = M
        self._C = C
        self._K = K
        self.name = name
        # Values for evalues and evectors will be calculated by self._calc_system
        self.evalues = None
        self.evectors = None
        self.wn = None
        self.wd = None
        self.H = None

        self.n = len(M)
        self._calc_system()

    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, value):
        self._M = value
        # if the parameter is changed this will update the system
        self._calc_system()

    @property
    def C(self):
        return self._C

    @C.setter
    def C(self, value):
        self._C = value
        # if the parameter is changed this will update the system
        self._calc_system()

    @property
    def K(self):
        return self._K

    @K.setter
    def K(self, value):
        self._K = value
        # if the parameter is changed this will update the system
        self._calc_system()

    def _calc_system(self):
        self.evalues, self.evectors = self._eigen()
        self.wn = (np.absolute(self.evalues)/(2*np.pi))[:self.n]
        self.wd = (np.imag(self.evalues)/(2*np.pi))[:self.n]
        self.H = self._H()

    def A(self):
        """State space matrix"""
        Z = np.zeros((self.n, self.n))
        I = np.eye(self.n)

        A = np.vstack(
            [np.hstack([Z, I]),
             np.hstack([la.solve(-self.M, self.K), la.solve(-self.M, self.C)])])
        return A

    @staticmethod
    def _index(eigenvalues):
        r"""Function used to generate an index that will sort
        eigenvalues and eigenvectors based on the imaginary (wd)
        part of the eigenvalues. Positive eigenvalues will be
        positioned at the first half of the array.
        """
        # avoid float point errors when sorting
        evals_truncated = np.around(eigenvalues, decimals=10)
        a = np.imag(evals_truncated)  # First column
        b = np.absolute(evals_truncated)  # Second column
        ind = np.lexsort((b, a))  # Sort by imag, then by absolute
        # Positive eigenvalues first
        positive = [i for i in ind[len(a) // 2:]]
        negative = [i for i in ind[:len(a) // 2]]

        idx = np.array([positive, negative]).flatten()

        return idx

    def _eigen(self, sorted_=True):
        r"""This method will return the eigenvalues and eigenvectors of
        the state space matrix A, sorted by the index method which
        considers the imaginary part (wd) of the eigenvalues for sorting.
        To avoid sorting use sorted_=False
        """
        evalues, evectors = la.eig(self.A())
        if sorted_ is False:
            return evalues, evectors

        idx = self._index(evalues)

        return evalues[idx], evectors[:, idx]

    def _H(self):
        r"""Continuous-time linear time invariant system.

        This method is used to create a Continuous-time linear
        time invariant system for the mdof system.
        From this system we can obtain poles, impulse response,
        generate a bode, etc.


        """
        Z = np.zeros((self.n, self.n))
        I = np.eye(self.n)

        # x' = Ax + Bu
        B2 = I
        A = self.A()
        B = np.vstack([Z,
                       la.solve(self.M, B2)])

        # y = Cx + Du
        # Observation matrices
        Cd = I
        Cv = Z
        Ca = Z

        C = np.hstack((Cd - Ca @ la.solve(self.M, self.K), Cv - Ca @ la.solve(self.M, self.C)))
        D = Ca @ la.solve(self.M, B2)

        sys = signal.lti(A, B, C, D)

        return sys

    def time_response(self, F, t, ic=None):
        r"""Time response for a mdof system.

        This method returns the time response for a mdof system
        given a force, time and initial conditions.

        Parameters
        ----------
        F : array
            Force array (needs to have the same length as time array).
        t : array
            Time array.
        ic : array
            The initial conditions on the state vector (zero by default).

        Returns
        ----------
        t : array
            Time values for the output.
        yout : array
            System response.
        xout : array
            Time evolution of the state vector.


        Examples
        --------
        >>> m1, m2 = 1, 1
        >>> c1, c2, c3 = 1, 1, 1
        >>> k1, k2, k3 = 1e3, 1e3, 1e3

        >>> M = np.array([[m1, 0],
        ...               [0, m2]])
        >>> C = np.array([[c1+c2, -c2],
        ...               [-c2, c2+c3]])
        >>> K = np.array([[k1+k2, -k2],
        ...               [-k2, k2+k3]])
        >>> sys1 = VibSystem(M, C, K) # create the system
        >>> t = np.linspace(0, 25, 1000) # time array
        >>> F2 = np.zeros((len(t), 2))
        >>> F2[:, 1] = 1000*np.sin(40*t) # force applied on m2
        >>> t, yout, xout = sys1.time_response(F2, t)
        >>> yout[:5, 0] # response on m1
        array([ 0.        ,  0.00304585,  0.07034816,  0.32048951,  0.60748282])
        >>> yout[:5, 1] # response on m2
        array([ 0.        ,  0.08160348,  0.46442657,  0.7885541 ,  0.47808092])
        """
        if ic is not None:
            return signal.lsim(self.H, F, t, ic)
        else:
            return signal.lsim(self.H, F, t)
