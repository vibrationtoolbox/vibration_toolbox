import scipy as sp
import matplotlib.pyplot as plt

__all__= ["frf"]


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

    w = sp.sin(sp.pi*sp.arange(len(f))/len(f))**2 # window
    # apply window
    xw = x*w
    fw = f*w
    # take ffts
    FX = sp.fft(xw)
    FF = sp.fft(fw)
    # calculate the spectral densities
    SXF = FF*sp.conj(FX)
    SXX = FX*sp.conj(FX)
    SFF = FF*sp.conj(FF)
    SFX = FX*sp.conj(FF)
    # calculate the transfer functions
    TXF = SXX/SXF
    TXF2 = SFX/SFF

    lt = len(TXF)//2
    freq = sp.arange(lt)/(2*lt*dt)

    TXF = TXF[:lt]
    mag = sp.absolute(TXF)
    ang = sp.angle(TXF)*180/sp.pi

    coh = (sp.absolute(SXF)**2)/(SXX*SFF)
    coh = sp.real(coh)

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
    ax3.set_ylim(0,2)

    ax1.semilogy(freq, mag)
    ax2.plot(freq, ang)
    ax3.plot(freq, coh[:lt])

    _ = plt.show()

    return freq, mag, ang, coh


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS| doctest.NORMALIZE_WHITESPACE|doctest.REPORT_NDIFF)

    """ What this does.
    python (name of this file)  -v
    will test all of the examples in the help.

    Leaving off -v will run the tests without any output. Success will return nothing.

    See the doctest section of the python manual.
    https://docs.python.org/3.5/library/doctest.html
    """
