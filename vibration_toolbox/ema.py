import numpy as np
import matplotlib.pyplot as plt

__all__ = ['frf', 'plot_fft']


def frf(x, f, dt):
    """
    This function will return the frequency response
    function (H(iw)) of the sampled data.

    Parameters
    ----------
    x: array
        Array with the displacement data
    f: array
        Array with the force data
    dt: float
        Time step of the sampled data
    n: int
        Number of points in the fft

    Returns
    ----------
    freq: array
        Driving frequencies
    mag: array
        Magnitude of the frequency response function
    ang: array
        Phase of the frequency response function
    coh: array
        Coherence function

        Plot with the frf magnitude, phase and
        coherence.
    Examples
    --------
    >>> # First we need to load the sampled data which in a .mat file
    >>> import vibration_toolbox as vt
    >>> import scipy.io as sio
    >>> data = sio.loadmat(vt.__path__[0] + '/data/frf_data1.mat')
    >>> #print(data)
    >>> # Data is imported as arrays. We need to modify then to fit our function
    >>> x = data['x']
    >>> x = x.reshape(len(x))
    >>> f = data['f']
    >>> f = f.reshape(len(f))
    >>> dt = data['dt']
    >>> dt = float(dt)
    >>> # Now we are able to call the function
    >>> freq, mag, ang, coh = vt.frf(x, f, dt)
    >>> mag[10]
    1.018394853080...
    """

    w = np.sin(np.pi*np.arange(len(f))/len(f))**2  # window
    # apply window
    xw = x*w
    fw = f*w
    # take ffts
    FX = np.fft.fft(xw)
    FF = np.fft.fft(fw)
    # calculate the spectral densities
    SXF = FF*np.conj(FX)
    SXX = FX*np.conj(FX)
    SFF = FF*np.conj(FF)
    SFX = FX*np.conj(FF)
    # calculate the transfer functions
    TXF = SXX/SXF
    TXF2 = SFX/SFF

    lt = len(TXF)//2
    freq = np.arange(lt)/(2*lt*dt)

    TXF = TXF[:lt]
    mag = np.absolute(TXF)
    ang = np.angle(TXF)*180/np.pi

    coh = (np.absolute(SXF)**2)/(SXX*SFF)
    coh = np.real(coh)

    # plot H(w)
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax1)
    ax3 = fig.add_subplot(313, sharex=ax2)
    fig.tight_layout()

    ax1.set_title('$H(\omega)$ - Magnitude')
    ax2.set_title('$H(\omega)$ - Phase')
    ax3.set_title('$H(\omega)$ - Coherence')
    ax3.set_xlabel('Frequency (Hz)')
    ax3.set_ylim(0, 2)

    ax1.semilogy(freq, mag)
    ax2.plot(freq, ang)
    ax3.plot(freq, coh[:lt])

    _ = plt.show()

    return freq, mag, ang, coh


def plot_fft(t, time_response, ax=None):
    """
    This function will plot the fft given a time vector
    and the system time response.

    Parameters
    ----------
    t : array
        Time array.
    time_response : array
        Array with the system's time response.

    Returns
    ----------
    ax : array with matplotlib.axes, optional
        Matplotlib axes array created with plt.subplots.
        Plot has frequency in rad/s and magnitude in meters peak to peak.

    Examples
    --------
    >>> t = np.linspace(0, 10, 1000)
    >>> time_response = 2 * np.sin(40*t)
    >>> plot_fft(t, time_response)
    <matplotlib.axes...
    """
    if ax is None:
        fig, ax = plt.subplots()

    Ts = t[1] - t[0]  # sampling interval
    Fs = 1 / Ts  # sampling rate

    n = len(time_response)  # length of the signal
    k = np.arange(n)
    T = n / Fs
    freq_range = (k / T)[:(n // 2)] * 2 * np.pi  # one side frequency range
    y = np.fft.fft(time_response) * 4 / n  # * 4 / n to normalize to pk-pk
    y = y[:(n // 2)]
    ax.plot(freq_range, abs(y))
    ax.set_xlim(freq_range[0], freq_range[-1])
    ax.set_xlabel('Freq (rad/s)')
    ax.set_ylabel('Amplitude (m - pk-pk)')

    return ax



