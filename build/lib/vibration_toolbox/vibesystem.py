import numpy as np
import scipy.linalg as la
from scipy import signal

__all__ = ['VibeSystem']


class VibeSystem(object):
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
        name : str, optional
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
        >>> sys1 = VibeSystem(M, C, K)
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
        """Mass Matrix"""
        return self._M

    @M.setter
    def M(self, value):
        self._M = value
        # if the parameter is changed this will update the system
        self._calc_system()

    @property
    def C(self):
        """Damping Matrix"""
        return self._C

    @C.setter
    def C(self, value):
        self._C = value
        # if the parameter is changed this will update the system
        self._calc_system()

    @property
    def K(self):
        """Stiffness Matrix"""
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
        """State space matrix
        """
        Z = np.zeros((self.n, self.n))
        I = np.eye(self.n)

        A = np.vstack(
            [np.hstack([Z, I]),
             np.hstack([la.solve(-self.M, self.K),
                        la.solve(-self.M, self.C)])])
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

        C = np.hstack((Cd - Ca @ la.solve(self.M, self.K),
                       Cv - Ca @ la.solve(self.M, self.C)))
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
        ic : array, optional
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
        >>> sys1 = VibeSystem(M, C, K) # create the system
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

    def freq_response(self, omega=None, modes=None):
        r"""Frequency response for a mdof system.

        This method returns the frequency response for a mdof system
        given a range of frequencies and the modes that will be used.

        Parameters
        ----------
        omega : array, optional
            Array with the desired range of frequencies (the default
             is 0 to 1.5 x highest damped natural frequency.
        modes : list, optional
            Modes that will be used to calculate the frequency response
            (all modes will be used if a list is not given).

        Returns
        ----------
        omega : array
            Array with the frequencies
        magdb : array
            Magnitude (dB) of the frequency response for each pair input/output.
            The order of the array is: [output, input, magnitude]
        phase : array
            Phase of the frequency response for each pair input/output.
            The order of the array is: [output, input, phase]

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
        >>> sys1 = VibeSystem(M, C, K) # create the system
        >>> omega, magdb, phase = sys1.freq_response()
        >>> magdb[0, 1, :4] # magnitude for output on 0 and input on 1.
        array([-69.54242509, -69.54234685, -69.54211212, -69.5417209 ])
        >>> np.around(phase[1, 1, :4],5) # phase for output on 1 and input on 1.
        array([-0.     , -0.00471, -0.00942, -0.01413])
        """

        rows = self.H.inputs  # inputs (mag and phase)
        cols = self.H.inputs  # outputs

        B = self.H.B
        C = self.H.C
        D = self.H.D

        evals = self.evalues
        psi = self.evectors
        psi_inv = la.inv(psi)  # TODO change to get psi_inv from la.eig

        # if omega is not given, define a range
        if omega is None:
            omega = np.linspace(0, max(evals.imag) * 1.5, 1000)

        # if modes are selected:
            if modes is not None:
                n = self.n  # n dof -> number of modes
                m = len(modes)  # -> number of desired modes
                # idx to get each evalue/evector and its conjugate
                idx = np.zeros((2 * m), int)
                idx[0:m] = modes  # modes
                idx[m:] = range(2 * n)[-m:]  # conjugates (see how evalues are ordered)

                evals_m = evals[np.ix_(idx)]
                psi_m = psi[np.ix_(range(2 * n), idx)]
                psi_inv_m = psi_inv[np.ix_(idx, range(2 * n))]

                magdb_m = np.empty((cols, rows, len(omega)))
                phase_m = np.empty((cols, rows, len(omega)))

                for wi, w in enumerate(omega):
                    diag = np.diag([1 / (1j * w - lam) for lam in evals_m])
                    H = C @ psi_m @ diag @ psi_inv_m @ B + D

                    magh = 20.0 * np.log10(abs(H))
                    angh = np.rad2deg((np.angle(H)))

                    magdb_m[:, :, wi] = magh
                    phase_m[:, :, wi] = angh

                return omega, magdb_m, phase_m

        magdb = np.empty((cols, rows, len(omega)))
        phase = np.empty((cols, rows, len(omega)))

        for wi, w in enumerate(omega):
            diag = np.diag([1 / (1j * w - lam) for lam in evals])
            H = C @ psi @ diag @ psi_inv @ B + D

            magh = 20.0 * np.log10(abs(H))
            angh = np.rad2deg((np.angle(H)))

            magdb[:, :, wi] = magh
            phase[:, :, wi] = angh

        return omega, magdb, phase
