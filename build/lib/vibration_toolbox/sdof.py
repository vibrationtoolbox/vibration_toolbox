import matplotlib.pyplot as plt
import numpy as np
from scipy import fft
import matplotlib as mpl
from scipy import integrate

try:
    from IPython.display import clear_output, display, HTML
    from ipywidgets.widgets.interaction import interact, interactive
except ImportError:
    print('Interactive iPython tools will not work without IPython.display \
          and ipywidgets installed.')


def in_ipynb():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':  # Jupyter notebook or qtconsole?
            return True
        elif shell == 'TerminalInteractiveShell':  # Terminal running IPython?
            return False
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['figure.figsize'] = (10, 6)


def free_response(m=10, c=1, k=100, x0=1, v0=-1, max_time=10):
    """Free response of a second order linear oscillator.

    Returns t, x, v, zeta, omega, omega_d and A resulting from the
    free response of a second order linear ordinary differential
    equation defined by
    :math:`m\ddot{x} + c \dot{x} + k x = 0`
    given initial conditions :math:`x_0` and :math:`\dot{x}_0 = v_0` for
    :math:`0 < t < t_{max}`

    Parameters
    ----------
    m, c, k :  floats, optional
        mass, damping coefficient, stiffness
    x0, v0:  floats, optional
        initial displacement, initial velocity
    max_time: float, optional
        end time for :math:`x(t)`

    Returns
    -------
    t, x, v : ndarrays
        time, displacement, and velocity
    zeta, omega, omega_d, A : floats
        damping ratio, undamped natural frequency, damped natural frequency,
        Amplitude

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import vibration_toolbox as vtb
    >>> vtb.free_response()[1][:5] # get the first five values of x
    array([[ 1.  ],
           [ 1.  ],
           [ 0.99],
           [ 0.99],
           [ 0.98]])

    >>> t, x, *_ = vtb.free_response() # *_ ignores all other returns
    >>> plt.plot(t,x)
    [<matplotlib.lines.Line2D object at ...>]
    >>> plt.xlabel('Time (sec)')
    <matplotlib.text.Text object at ...>
    >>> plt.ylabel('Displacement (m)')
    <matplotlib.text.Text object at ...>
    >>> plt.title('Displacement versus time')
    <matplotlib.text.Text object at ...>
    >>> plt.grid('on')
    """

    omega = np.sqrt(k / m)
    zeta = c / 2 / omega / m
    omega_d = omega * np.sqrt(1 - zeta ** 2)
    A = np.sqrt(x0 ** 2 + (v0 + omega * zeta * x0) ** 2 / omega_d ** 2)

    def sdofs_deriv(x_xd, t0, m=m, c=c, k=k):
        x, xd = x_xd
        return [xd, -c / m * xd - k / m * x]

    z0 = np.array([[x0, v0]])
    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250 * max_time))
    z_t = np.asarray([integrate.odeint(sdofs_deriv, z0i, t)
                      for z0i in z0])

    x, y = z_t[:, :].T
    return t, x, y, zeta, omega, omega_d, A


def phase_plot(m=10, c=1, k=100, x0=1, v0=-1, max_time=10):
    """Phase plot of free response of single degree of freedom system.

    For information on variables see `free_response`.

    Parameters
    ----------
    m, c, k :  floats, optional
        mass, damping coefficient, stiffness
    x0, v0:  floats, optional
        initial displacement, initial velocity
    max_time: float, optional
        end time for :math:`x(t)`

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import vibration_toolbox as vtb
    >>> vtb.phase_plot() # *_ ignores all other returns
    """
    t, x, v, zeta, omega, omega_d, A = free_response(
        m, c, k, x0, v0, max_time)
    fig = plt.figure()
    fig.suptitle('Velocity vs Displacement')
    ax = fig.add_subplot(111)
    ax.set_xlabel('Displacement')
    ax.set_ylabel('Velocity')
    ax.grid('on')
    ax.plot(x, v)


def phase_plot_i(max_time=(1.0, 200.0), v0=(-100, 100, 1.0),
                 m=(1.0, 100.0, 1.0),
                 c=(0.0, 1.0, 0.1), x0=(-100, 100, 1), k=(1.0, 100.0, 1.0)):
    """Interactive phase plot of free response of single degree of freedom system.
    ``phase_plot_i`` is only functional in a
    `Jupyter notebook <http://jupyter.org>`_.

    Parameters
    ----------
    m, c, k :  floats, optional
        mass, damping coefficient, stiffness
    x0, v0:  floats, optional
        initial displacement, initial velocity
    max_time: float, optional
        end time for :math:`x(t)`

    """
    if in_ipynb():
        w = interactive(phase_plot, max_time=max_time, v0=v0, m=m,
                        c=c, x0=x0, k=k)
        display(w)
    else:
        print('phase_plot_i can only be used in an iPython notebook.')


def time_plot(m=10, c=1, k=100, x0=1, v0=-1, max_time=100):
    t, x, v, zeta, omega, omega_d, A = free_response(
        m, c, k, x0, v0, max_time)
    fig = plt.figure()
    fig.suptitle('Displacement vs Time')
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time')
    ax.set_ylabel('Displacement')
    ax.grid('on')
    ax.plot(t, x)
    if zeta < 1:
        ax.plot(t, A * np.exp(-zeta * omega * t), '--g',
                linewidth=1)
        ax.plot(t, -A * np.exp(-zeta * omega * t), '--g',
                linewidth=1, label='$A e^{- \zeta \omega t}$')
        tmin, tmax, xmin, xmax = ax.axis()
        ax.text(.75 * tmax, .85 * (xmax - xmin) + xmin,
                '$\omega$ = %0.2f rad/sec' % (omega))
        ax.text(.75 * tmax, .80 * (xmax - xmin) +
                xmin, '$\zeta$ = %0.2f' % (zeta))
        ax.text(.75 * tmax, .75 * (xmax - xmin) + xmin,
                '$\omega_d$ = %0.2f rad/sec' % (omega_d))
    else:
        tmin, tmax, xmin, xmax = ax.axis()
        ax.text(.75 * tmax, .85 * (xmax - xmin) +
                xmin, '$\zeta$ = %0.2f' % (zeta))
        ax.text(.75 * tmax, .80 * (xmax - xmin) + xmin,
                '$\lambda_1$ = %0.2f' %
                (zeta * omega - omega * (zeta ** 2 - 1)))
        ax.text(.75 * tmax, .75 * (xmax - xmin) + xmin,
                '$\lambda_2$ = %0.2f' %
                (zeta * omega + omega * (zeta ** 2 - 1)))
    ax.legend()


def time_plot_i(max_time=(1.0, 100.0), x0=(-100, 100), v0=(-100, 100),
                m=(1.0, 100.0), c=(0.0, 100.0), k=(1.0, 100.0)):
    """Interactive single degree of freedom free reponse plot in iPython

    ``time_plot_i`` is only functional in a
    `Jupyter notebook <http://jupyter.org>`_.

    Parameters
    ----------
    m, c, k :  floats, optional
        mass, damping coefficient, stiffness
    x0, v0:  floats, optional
        initial displacement, initial velocity
    max_time: float, optional
        end time for :math:`x(t)`

    """
    if in_ipynb():
        w = interactive(time_plot, max_time=max_time, v0=v0, m=m,
                        c=c, x0=x0, k=k)
        display(w)
    else:
        print('time_plot_i can only be used in an iPython notebook.')
    '''
    I'd like to get the sliders to be side by side to take less vertical
    space
    # cont = widgets.HBox(children = w)
    print(help(w))
    '''


def analytical(m=1, c=0.1, k=1, x0=1, v0=0, n=8, dt=0.05):

    w = np.sqrt(k / m)
    zeta = c / (2 * w * m)  # (1.30)

    wd = w * np.sqrt(1 - zeta**2)  # (1.37)
    t = np.linspace(0, n * dt, n + 1)

    print('The natural frequency is ', w, 'rad/s.')
    print('The damping ratio is ', zeta)
    print('The damped natural frequency is ', wd)

    if zeta < 1:
        A = np.sqrt(((v0 + zeta * w * x0)**2 + (x0 * wd)**2) / wd**2)  # (1.38)
        phi = np.arctan2(x0 * wd, v0 + zeta * w * x0)  # (1.38)
        x = A * np.exp(-zeta * w * t) * np.sin(wd * t + phi)  # (1.36)
        print('A =', A)
        print('phi =', phi)

    elif zeta == 1:
        a1 = x0  # (1.46)
        a2 = v0 + w * x0  # (1.46)
        print('a1= ', a1)
        print('a2= ', a2)
        x = (a1 + a2 * t) * np.exp(-w * t)  # (1.45)

    else:
        a1 = (-v0 + (-zeta + np.sqrt(zeta**2 - 1)) * w * x0) / \
            (2 * w * np.sqrt(zeta**2 - 1))  # (1.42)
        a2 = (v0 + (zeta + np.sqrt(zeta**2 - 1)) * w * x0) / \
            (2 * w * np.sqrt(zeta**2 - 1))  # (1.43)
        print('a1= ', a1)
        print('a2= ', a2)
        x = (np.exp(-zeta * w * t) *
             (a1 * np.exp(-w * np.sqrt(zeta**2 - 1) * t) +
             a2 * np.exp(w * np.sqrt(zeta**2 - 1) * t)))  # (1.41)

    return x


def euler(m=1, c=.1, k=1, x0=1, v0=0, n=8, dt=0.05):
    """Euler method free response of a SDOF system.

    Free response using Euler's method to perform numerical integration.

    Parameters
    ----------
    m, c, k: float
        Mass, damping and stiffness.
    x0, v0: float
        Initial conditions
    n: int
        The number of steps
    dt: float
        The step size.

    Returns
    -------
    t, x, v: array
        Time, displacement, and velocity

    Examples
    --------
    >>> euler(m=1, c=.1, k=1, x0=1, v0=0, n=8, dt=0.05)
    (array([ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ]), array([[ 1.  ,  0.  ],
           [ 1.  , -0.05],
           [ 1.  , -0.1 ],
           [ 0.99, -0.15],
           [ 0.99, -0.2 ],
           [ 0.98, -0.25],
           [ 0.96, -0.29],
           [ 0.95, -0.34],
           [ 0.93, -0.39]]))
    """

    # creates the state space matrix
    A = np.array([[0, 1],
                  [-k / m, -c / m]])
    # creates the x array and set the first line according to the initial
    # conditions
    x = np.zeros((n + 1, 2))
    x[0] = x0, v0

    for i in range(0, n):
        x[i + 1] = x[i] + dt * A@x[i]

    t = np.linspace(0, n * dt, n + 1)

    return t, x


def rk4(m=1, c=.1, k=1, x0=1, v0=0, n=8, dt=0.05):
    """Runge-Kutta solution of underdamped system

    Parameters
    ----------
    m, c, k: float
        Mass, damping and stiffness.
    x0, v0: float
        Initial conditions
    n: int
        The number of steps
    dt: float
        The step size.

    Returns
    -------
    t, x, v: array
        Time, displacement, and velocity

    Examples
    --------
    >>> rk4(m=1, c=.1, k=1, x0=1, v0=0, n=8, dt=0.05)
    (array([ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ]), array([[ 1.  ,  0.  ],
           [ 1.  , -0.05],
           [ 1.  , -0.1 ],
           [ 0.99, -0.15],
           [ 0.98, -0.2 ],
           [ 0.97, -0.24],
           [ 0.96, -0.29],
           [ 0.94, -0.34],
           [ 0.92, -0.38]]))
    """

    t = np.linspace(0, n * dt, n + 1)
    x = np.zeros((n + 1, 2))
    x[0, :] = x0, v0
    A = np.array([[0, 1],
                  [-k / m, -c / m]])

    def f(x_): return A@x_

    for i in range(n):
        k1 = dt * f(x[i])
        k2 = dt * f(x[i] + k1 / 2)
        k3 = dt * f(x[i] + k2 / 2)
        k4 = dt * f(x[i] + k3)
        x[i + 1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

    return t, x


# def
# euler_beam_frf(xin=0.22,xout=0.22,fmin=0.0,fmax=1000.0,beamparams=sp.array((7.31e10,
# 1/12*0.03*.015**3, 2747, .015*0.03, 0.4)),


def frfplot(f, H):
    plt.subplot(211)
    plt.plot(f, 20 * np.log10(np.absolute(np.sum(H, axis=1))), '-')
    plt.grid('on')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('FRF (dB)')
    axlim = plt.axis()
    plt.axis(axlim + np.array([0, 0, -0.1 * (axlim[3] -
             axlim[2]), 0.1 * (axlim[3] - axlim[2])]))

    plt.subplot(212)
    plt.plot(f, np.unwrap(np.angle(np.sum(H, axis=1))) / np.pi * 180, '-')
    plt.grid('on')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase (deg)')
    axlim = plt.axis()
    plt.axis(axlim + np.array([0, 0, -0.1 * (axlim[3] -
             axlim[2]), 0.1 * (axlim[3] - axlim[2])]))


def response(xdd, f, t, x0, v0):
    """
    returns t, x, v
    :math:`\ddot{x} = g(x,v) + f(t)`
    given initial conditions :math:`x_0` and :math:`\dot{x}_0 = v_0` for the time `t`

*** Function hasn't been written yet... work in progress. Mimics vtb1_3
    Parameters

    m, c, k:           1) Floats. Mass, damping and stiffness.
    x0, v0:            2) Floats. Initial conditions
    max_time:          3) Float. end time or response to be returned

    Returns

    t, x, v: 1) Arrays. Time, displacement, and velocity

    :Example:
    >>> free_response()[1][:5] # get the first five values of x
    array([[ 1.  ],
           [ 1.  ],
           [ 0.99],
           [ 0.99],
           [ 0.98]])
    """

    omega = np.sqrt(k / m)
    zeta = c / 2 / omega / m
    omega_d = omega * np.sqrt(1 - zeta**2)
    A = np.sqrt(x0**2 + (v0 + omega * zeta * x0)**2 / omega_d**2)
#    print('The natural frequency is ', omega, 'rad/s.');
#    print('The damping ratio is ', zeta);
#    print('The damped natural frequency is ', omega_d);

    def sdofs_deriv(x_xd, t0, m=m, c=c, k=k):
        x, xd = x_xd
        return [xd, -c / m * xd - k / m * x]

    z0 = np.array([[x0, v0]])
    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250 * max_time))
    z_t = np.asarray([integrate.odeint(sdofs_deriv, z0i, t)
                      for z0i in z0])

    x, y = z_t[:, :].T
    return t, x, y, zeta, omega, omega_d, A


def forced_analytical(m=10, k=100, x0=1, v0=0,
                      wdr=0.5, F0=10, tf=100):
    """
    Returns the response of an undamped single degree of freedom system
    to a sinusoidal input with amplitude F0 and frequency wdr.

    Parameters
    ----------
    m, k: float
        Mass and stiffness
    x0, v0: float
        Initial conditions
    wdr:
        Force frequency
    F0: float
        Force magnitude
    tf: float
        End time

    Returns
    -------
    t, x: array
        Time and displacement

    Examples
    --------
    >>> forced_analytical(m=10, k=100, x0=1, v0=0, wdr=0.5, F0=10, tf=100)
    (array([   0.,    0.,    0., ...,  100.,  100.,  100.]), array([ 1.  ,  1.  ,  1.  , ..., -0.33, -0.33, -0.33]))
    """

    t = np.linspace(0, tf, tf/0.000125)

    f0 = F0 / m
    w = np.sqrt(k / m)
    x = (v0 / w * np.sin(w * t) +
         (x0 - f0 / (w**2 - wdr**2)) * np.cos(w * t) +
         f0 / (w**2 - wdr**2) * np.cos(wdr * t))   # (2.11)

    return t, x


def forced_response(m=10, c=0, k=100, x0=1, v0=0,
                    wdr=0.5, F0=10, max_time=100):
    """Harmonic response of SDOF system.

    Returns the the response of an underdamped single degree of
    freedom system to a sinusoidal input with amplitude F0 and
    frequency :math:`\omega_{dr}`.

    Parameters
    ----------
    m, c, k: float, optional
        Mass Damping, and stiffness
    x0, v0: float, optional
        Initial conditions
    wdr: float, optional
        Force frequency
    F0: float, optional
        Force magnitude
    max_time: float, optional
        End time

    Returns
    -------
    t, x, y: array
        Time, displacement and velocity

    Examples
    --------
    >>> f = forced_response(m=10, c=0, k=100, x0=1, v0=0, wdr=0.5, F0=10, max_time=100)
    >>> f[0][0]
    0.0
    """

    def sdofs_deriv(x_xd, t, m=m, c=c, k=k):
        x, xd = x_xd
        return [xd, (F0*np.cos(wdr*t) / m) - (c / m) * xd - (k / m) * x]

    z0 = np.array([x0, v0])
    # Solve for the trajectories
    t = np.linspace(0, max_time, int(250 * max_time))
    z_t = integrate.odeint(sdofs_deriv, z0, t)

    x, y = z_t.T
    return t, x, y


def steady_state_response(zs, rmin, rmax):
    """Plot steady state response SDOF damped system.

    Parameters
    ----------
    zs: array
        Array with the damping values
    rmin, rmax: floats
        Minimum and maximum frequency ratio

    Returns
    -------
    r: Array
        Array containing the values for the frequency ratio
    A: Array
        Array containing the values for anmplitude

        Plot with steady state magnitude and phase

    Examples
    --------
    >>> r, A = steady_state_response([0.1, 0.3, 0.8], 0, 2)
    >>> A[10]
    (0.98423159842039087-0.15988334018879749j)
    """

    if not isinstance(zs, list):
        zs = [zs]
    r = np.linspace(rmin, rmax, 100*(rmax-rmin))
    A0 = np.zeros((len(zs), len(r)), complex)
    for z in enumerate(zs):
        A0[z[0]] = (1/(1 - r**2 + 2*1j*r*z[1]))

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    plt.tight_layout()

    ax1.set_ylabel('Normalized Amplitude (dB)')
    ax1.set_title('Normalized Amplitude vs Frequency Ratio')

    ax2.set_xlabel('Frequency Ratio')
    ax2.set_ylabel('Phase Lag (deg)')
    ax2.set_title('Phase vs Frequency Ratio')

    for A in A0:
        ax1.plot(r, (np.absolute(A)))
        ax2.plot(r, -np.angle(A)/np.pi*180)

    ax1.legend((['$\zeta$ = ' + (str(s)) for s in zs]))

    return r, A


def steady_state_response_i(zs=(0, 1.0, 0.1), rmin=0, rmax=2.0):
    """Interactive phase plot of steady state response of
     single degree of freedom system.
    ``steady_state_response`` is only functional in a
    `Jupyter notebook <http://jupyter.org>`_.


    Parameters
    ----------
    zs: array
        Array with the damping values
    rmin, rmax: floats
        Minimum and maximum frequency ratio

    Returns
    -------
    r: Array
        Array containing the values for the frequency ratio
    A: Array
        Array containing the values for anmplitude

        Plot with steady state magnitude and phase
    """
    if in_ipynb():
        w = interactive(steady_state_response, zs=zs,
                        rmin=rmin, rmax=rmax)
        display(w)
    else:
        print('steady_state_response_i can only be used in an iPython notebook.')


def transmissibility(zs, rmin, rmax):
    """
    Returns a plot Displacement transmissibility ratio
    and force transmissibility ratio of a single degree
    of freedom damped system.

    Parameters
    ----------
    zs: array
        Array with the damping values
    rmin, rmax: float
        Minimum and maximum frequency ratio

    Returns
    -------
    r: Array
        Array containing the values for the frequency ratio
    D: Array
        Array containing the values for displacement
    F: Array
        Array containing the values for force

        Plot with Displacement transmissibility ratio
        and force transmissibility ratio

    Examples
    --------
    >>> r, D, F = transmissibility([0.01, 0.05, 0.1, 0.25, 0.5, 0.7], 0, 2)
    >>> D[10]
    1.0100027508815634
    """

    if not isinstance(zs, list):
        zs = [zs]
    r = np.linspace(rmin, rmax, 100*(rmax-rmin))
    DT = np.zeros((len(zs), len(r)))
    for z in enumerate(zs):
        # 2.71
        DT[z[0]] = ((1 + (2 * z[1] * r)**2) /
                    ((1 - r**2)**2 + (2 * z[1] * r)**2))**0.5

    FT = (r**2)*DT

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    plt.tight_layout()

    ax1.set_ylabel('Displacement Transmissibility Ratio (dB)')
    ax1.set_title('Displacement Transmissibility Ratio vs Frequency Ratio (X/Y)')
    ax1.set_ylim(0, 6)

    ax2.set_xlabel('Frequency Ratio')
    ax2.set_ylabel('Force Transmissibility Ratio (dB)')
    ax2.set_title('Force Transmissibility Ratio versus Frequency Ratio (F_T/kY)')
    ax2.set_yscale("log")
    ax2.set_ylim(0.01, 50)

    for D in DT:
        ax1.plot(r, D)
    for F in FT:
        ax2.plot(r, F)

    ax1.legend((['$\zeta$ = ' + (str(s)) for s in zs]))

    return r, D, F


def rotating_unbalance(m, m0, e, zs, rmin, rmax, normalized=True):
    """
    Returns a plot Displacement of a system with rotating
    unbalance.

    Parameters
    ----------
    m: float
        Mass of the system
    m0, e: float
        Mass and eccentricity of the unbalance.
    zs: array
        Array with the damping values
    rmin, rmax: float
        Minimum and maximum frequency ratio
    normalized: bool
        If true, the displacement is normalized (m*X/(m0*e))

    Returns
    -------
    r: Array
        Array containing the values for the frequency ratio
    Xn: Array
        Array containing the values for displacement

        Plot with Displacement displacement and phase
        for a system with rotating unbalance.

    Examples
    --------
    >>> r, Xn = rotating_unbalance(m=1, m0=0.5, e=0.1, zs=[0.1, 0.25, 0.707, 1], rmin=0, rmax=3.5, normalized=True)
    >>> Xn[1][10]
    (0.10104614704226758-0.0051182602098315527j)
    """

    if not isinstance(zs, list):
        zs = [zs]
    r = np.linspace(rmin, rmax, 100*(rmax-rmin))
    Xn = np.zeros((len(zs), len(r)), complex)
    for z in enumerate(zs):
        Xn[z[0]] = (r / (1 - r**2 + 2*1j*r*z[1]))

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    plt.tight_layout()

    if normalized==False:
        Xn = Xn * (m0 * e / m)
        ax1.set_ylabel('Displacement Magnitude')
        ax1.set_title('Displacement Magnitude vs Frequency Ratio')
    else:
        ax1.set_ylabel('Normalized Displacement Magnitude')
        ax1.set_title('Normalized Displacement Magnitude vs Frequency Ratio')

    ax2.set_xlabel('Frequency Ratio')
    ax2.set_ylabel('Phase')
    ax2.set_title('Phase vs Frequency Ratio')

    for X_z in Xn:
        ax1.plot(r, np.absolute(X_z))
        ax2.plot(r, -np.angle(X_z)/np.pi*180)

    ax1.legend((['$\zeta$ = ' + (str(s)) for s in zs]))

    return r, Xn


def impulse_response(m, c, k, Fo, max_time):
    """
    Returns a plot with the response of the system to an
    impulse of magnitude Fo (N.s).

    Parameters
    ----------
    m, c, k: float
        Mass, damping and stiffness.
    Fo: float
        Force applied over time (units N.s)
    max_time: float
        End time

    Returns
    -------
    t: Array
        Array containing the values for the time
    x: Array
        Array containing the values for displacement

        Plot with the response of the system to an
        impulse of magnitude Fo (N.s).

    Examples
    --------
    >>> t, x = impulse_response(m=100, c=20, k=2000, Fo=10, max_time=100)
    >>> x[10]
    0.0039629845398805623
    """

    t = np.linspace(0, max_time, int(250 * max_time))

    wn = np.sqrt(k / m)
    zeta = c / (2 * wn * m)
    wd = wn * np.sqrt(1 - zeta**2)
    fo = Fo / m

    x = fo / (wd * np.exp(zeta * wn * t)) * np.sin(wd * t)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Displacement')
    ax1.set_title('Displacement vs Time')
    ax1.plot(t, x)

    return t, x


def step_response(m, c, k, Fo, max_time):
    """
    Returns a plot with the response of the system to an
    step of magnitude Fo.

    Parameters
    ----------
    m, c, k: float
        Mass, damping and stiffness.
    Fo: float
        Force applied
    max_time: float
        End time

    Returns
    -------
    t: Array
        Array containing the values for the time
    x: Array
        Array containing the values for displacement

        Plot with the response of the system to an
        step of magnitude Fo.

    Examples
    --------
    >>> t, x = step_response(m=100, c=20, k=2000, Fo=10, max_time=100)
    >>> x[10]
    7.9581008173000833e-05
    """

    t = np.linspace(0, max_time, int(250 * max_time))

    wn = np.sqrt(k / m)
    zeta = c / (2 * wn * m)
    wd = wn * np.sqrt(1 - zeta**2)
    fo = Fo / m

    if zeta != 1:
        phi = np.arctan(zeta / np.sqrt(1 - zeta**2))

    if 0 < zeta < 1:
        x = fo / wn**2 * (1 - wn / wd*np.exp(-zeta * wn * t)*np.cos(wd * t - phi))
    elif zeta == 1:
        lam = -wn
        A1 = -fo / wn**2
        A2 = -A1 * lam
        x = fo / wn**2 + A1 * np.exp(lam * t) + A2 * t * np.exp(lam * t)
    elif zeta > 1:
        lam1 = -zeta * wn - wn * np.sqrt(zeta**2 - 1)
        lam2 = -zeta * wn + wn * np.sqrt(zeta**2 - 1)
        A2 = fo / (wn**2 * (lam2 / lam1 - 1))
        A1 = -lam2 / lam1 * A2
        x = fo / wn**2 + A1 * np.exp(lam1 * t) + A2 * np.exp(lam2 * t)
    else:
        raise ValueError('Zeta should be greater than zero')

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Displacement')
    ax1.set_title('Displacement vs Time')
    ax1.plot(t, x)

    return t, x


def fourier_series(dat, t, n):
    """
    Fourier series approximation to a function.
    returns Fourier coefficients of a function.
    The coefficients are numerical approximations of the true
    coefficients.

    Parameters
    ----------
    dat: array
        Array of data representing the function.
    t: array
        Corresponding time array.
    n: int
        The desired number of terms to use in the
        Fourier series.

    Returns
    -------
    a, b: tuple
        Tuple containing arrays with the Fourier coefficients.
        The function also produces a plot of the approximation.

    Examples
    --------
    >>> f = np.hstack((np.arange(-1, 1, .04), np.arange(1, -1, -.04)))
    >>> f += 1
    >>> t = np.arange(0, len(f))/len(f)
    >>> a, b = fourier_series(f, t, 5)
    >>> a[0]
    2.0
    """

    len_ = len(dat)/2
    fs = (fft(dat))/len_
    a0 = fs[0]
    a = np.real(np.hstack((a0, fs[1:len(fs/2)])))
    b = -np.imag(fs[1:len(fs/2)])
    len_ *= 2
    dt = 2*np.pi/len_
    tp = np.arange(0, 2*np.pi, dt)
    dataapprox = a[0]/2 + np.zeros_like(dat)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t, dat)

    for i in range(1, n):
        newdat = a[i]*np.cos(tp*i) + b[i]*np.sin(tp*i)
        dataapprox += newdat
        if i == n-1:
            ax1.plot(t, newdat)

    ax1.plot(t, dataapprox)

    return a, b


def response_spectrum(f):
    """
    Will display the response spectrum to a partial ramp
    input(see Figure 3.13) for the system with natural frequency
    f (in Hz) and no damping.

    Parameters
    ----------
    f: float
        Natural frequency.

    Returns
    -------
    t, rs: tuple
        Tuple with time and response arrays. It also returns
        a plot with the response spectrum.

    Examples
    --------
    >>> t, rs = response_spectrum(10)
    >>> rs[10]
    1.6285602401720802
    """

    t = np.linspace(.001 * 4 / f, 10 / f, 200)
    w = 2 * np.pi * f

    one = np.ones_like(t)

    rs1 = one / (w * t)
    rs2 = np.sqrt(2 * (1 - np.cos(w * t)))

    rs = one + rs1 * rs2

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Rise time (t_1)')
    ax1.set_ylabel('Dimensionless maximum response - (xk/Fo)max')
    ax1.set_title('Response spectrum of a SDOF system with f = %s Hz' % f)
    ax1.plot(t, rs)

    return t, rs


def fourier_approximation(a0, aodd, aeven, bodd, beven, N, T):
    """
    Plot the Fourier series defined by:
    N is the number of terms.

    Parameters
    ----------
    a0: float or str
        a0 Fourier coefficient.
    aodd: float or str
        an Fourier coefficient for n odd.
    aeven: float or str
        an Fourier coefficient for n even.
    bodd: float or str
        bn Fourier coefficient for n odd
    beven: float or str
        bn Fourier coefficient for n even

    Returns
    -------
    t, F: tuple
        Tuple with time and F(t). It also returns
        a plot with the Fourier approximation.

    Examples
    --------
    >>> # Square wave
    >>> t, F = fourier_approximation(-1, 0, 0, '-3*(-1+(-1)**n)/n/pi', '-3*(-1+(-1)**n)/n/pi', 20, 2)
    >>> F[10]
    1.2697210294282535
    >>> # Triangular wave
    >>> t, F = fourier_approximation(0,'-8/pi**2/n**2',0,0,0,20,10)
    >>> F[10]
    -0.90234928911935097
    """
    args = [str(arg) for arg in [a0, aodd, aeven, bodd, beven]] # change to str
    a0, aodd, aeven, bodd, beven = args

    dt = min(T/400, T/10*N)
    t = np.arange(0, T*3, dt)
    F = 0*t + eval(a0)/2
    pi = np.pi

    for n in range(1, N):
        if n % 2 == 0:
            a = aeven
            b = beven
        else:
            a = aodd
            b = bodd
        F = F + eval(a)*np.cos(n*2*np.pi*t/T) + eval(b)*np.sin(n*2*np.pi*t/T)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time, t')
    ax1.set_ylabel('F(t)')
    ax1.plot(t, F)

    return t, F


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE)
    #import vibration_toolbox as vtb


    # doctest.run_docstring_examples(frfest,globals(),optionflags=doctest.ELLIPSIS)
    # doctest.run_docstring_examples(asd,globals(),optionflags=doctest.ELLIPSIS)
    """ What this does.
    python (name of this file)  -v
    will test all of the examples in the help.

    Leaving off -v will run the tests without any output. Success will return nothing.

    See the doctest section of the python manual.
    https://docs.python.org/3.5/library/doctest.html
    """
