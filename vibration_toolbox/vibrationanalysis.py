# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:32:44 2017

@author: Loranger
"""

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

def gen_frf(x,f,dt,n=256):
    """Generates a frequency response function (H(iw))
    
    Parameters
    ----------
    x: first signal (x(t))
    f: second signal (f(t))
    dt: the time step of the sampled data.
    n: is the number of points in the fft.
        Recommendedvalue for n is 256 unless zooming.
        Use a higher n for zooming.
        
    Returns
    -------
    
    Examples
    --------
    
    """
    
    #input for beta testing
    
    if n<128:
        n=256
        
    lf = np.size(f)
    #lfcount = np.arange(lf-1)
    w = np.transpose(np.sin(np.pi * np.arange(lf-1)/(lf-1))**2)
    x = x.flatten()
    f = f.flatten()
    xw = x * w
    fw = f * w
    FX = np.fft(xw)
    FF = np.fft(fw)
    SXF = FF * np.conj(FX)
    SXX = FX * np.conj(FX)
    SFF = FF * np.conj(FF)
    #SXF only used to find unused var TXF2, remove?
    #SFX = FX * np.conj(FF)
    TXF = SXX / SXF
    #TXF2 appears unused, keeping it here as it is in MATLAB version
    #TXF2 = SFX / SFF
    lTXF2 = np.size(TXF)/2
    freq = np.transpose(np.array([0,np.arange(lTXF2-1)/lTXF2/2/dt]))
    TXF = TXF[0:lTXF2]
    mag = abs(TXF)
    ang = np.angle(TXF)*180/np.pi
    coh = (abs(SXF)**2)/SXX/SFF
    
    #plot everything
    plt.subplot(3,1,1)
    plt.semilogy(freq,mag)
    plt.title('Txf - transfer function magnitude')
    plt.xlabel('Frequency (Hz)')
    plt.grid('on')
    plt.subplot(3,1,2)
    plt.plot(freq,ang)
    plt.title('Txf - phase')
    plt.xlabel('Frequency (Hz)')
    plt.grid('on')
    plt.subplot(3,1,3)
    plt.plot(freq,coh)
    plt.title('Txf - coherance')
    plt.xlabel('Frequency (Hz)')
    plt.grid('on')
    plt.show()
    
def mdof_test(f,TF,Fmin,Fmax):
    

def mdof_cf(f,TF,Fmin,Fmax):
    """Curve fit to multiple degree of freedom FRF
    
    Parameters
    ----------
    f: the frequency vector in Hz. Does not have to start at 0 Hz.
    TF: the complex transfer function
    Fmin: the minimum frequency to be used for curve fitting in the FRF
    Fmax: the maximum frequency to be used for curve fitting in the FRF
        
    Returns
    -------
    z, nf: the damping ratio and natural frequency (Hz)
    a: the numerator of the identified transfer functions
    
    """
    
    inlow = 0
    inhigh = np.size(f)
    
    if len()
    
    if f[inlow] == 0:
        inlow = 1
    #f = f
    
    R = np.transpose(TF)
    U, s, v = np.linalg.svd(R)
    T = U[:,1]
    Hp = np.transpose(T)*R
    R = np.transpose(Hp)
    
    ll = length(f)
    w = f*2*np.pi()*sqrt(-1)
    w2 =w*0
    R3 = R*0
    for i in range(0, ll):
        R3[i] = np.conj(R[ll+1-i])
        w2[i] = np.conj(w[ll+1-i])
        TF2[i,:] = np.conj(TF[ll+1-i,:])
    w = np.vstack((w2,2))
    R = np.vstack((R3,R))
    
    N = 2
    garbage, y = np.meshgrid(r_[0:N],R)
    x, w2d = np.meshgrid(r_[0:N],w)
    c = -1*w^(N)*R
    
    aa = []
    b = np.linalg.pinv(aa)*c
    #rs = np.roots([1,b[]])
    garbage, irs = np.sort(np.abs(np.imag(rs)))
    rs = rs[irs]
    omega = np.abs(rs[2])
    nf = omega/2/np.pi()
    XoF1 = 1/((rs[0] - w)*(rs[2] - w))
    XoF2 = 1/(w^0)
    XoF3 = 1/(w^2)
    XoF = [XoF1,XoF2,XoF3]
    TF3 = np.vstack(TF2,TF)
    a = np.linalg.pinv(XoF)*TF3
    u = np.transpose(a[0,:])
    u = u/np.sqrt(np.abs(a[1,1]))
    
    #plot stuff?
    
    return z, nf, u

def sdof_cf(f,TF,Fmin,Fmax):
    """Curve fit to single degree of freedom FRF
    
    Parameters
    ----------
    f: the frequency vector in Hz. Does not have to start at 0 Hz.
    TF: the complex transfer function
    Fmin: the minimum frequency to be used for curve fitting in the FRF
    Fmax: the maximum frequency to be used for curve fitting in the FRF
        
    Returns
    -------
    z, nf: the damping ratio and natural frequency (Hz)
    a: the numerator of the identified transfer functions
    
    """
    
    inlow = 0
    inhigh = np.size(f)
    
    #if
    
    if f[inlow] == 0:
        inlow = 1
    #f = f
    
    R = TF
    y,cin = np.max(np.abs(TF))
    
    ll = np.size(f)
    w = f*2*np.pi()*np.sqrt(-1)
    w2 = w*0
    R3 = R*0
    for i in range(0, ll):
        R3[i] = np.conj(R[ll+1-i])
        w2[i] = np.conj(w[ll+1-i])
    w = np.vstack(w2,w)
    R = np.vstack(R3,R)
    aa = []
    
    n = 1
    #N = 2*n
    N = 2
    x, y = np.meshgrid(r_[0:N],R)
    x, w2d = np.meshgrid(r_[0:N],w)
    c = -1*w^N*R
    aa = []
    
    np.linalg.division()
    #b = aa
    rs = np.roots([1,np.transpose(b[N-1:-1:0]+1)])
    
    srs,irs = np.sort(np.abs(np.imag(rs)))  
    
    rs = rs[irs]
    omega = np.abs(rs[2])
    z = -1*np.real(rs[2])/np.abs(np.imag(rs))
    nf = omega/2/np.pi()
    
    XoF1 = [l/(w-rw[1]), 1/(w-rs[2])]
    XoF2 = 1/(w^0)
    XoF3 =1/w^2
    XoF = [XoF1, XoF2, XoF3]
    #a = XoF\R
    XoF = XoF[np.arange(ll+1,2*ll)-2*0,:]*a
    a = np.sqrt(-2*np.imag(a[0])*np.imag(rs[0])-2*np.real(a[1])*np.real(rs[0]))
    Fmin = np.min(f)
    Fmax = np.max(f)
    phase = np.unwrap(np.angle(TF))*180/np.pi()
    phase2 = np.unwrap(np.angle(XoF))*180/np.pi()
    
    #plot stuff
    
    
    
    return z, nf, a

def second_order_tf(M,D,K,numin,numout,freq):
    
    freqout
    recep
    mobil
    inert
    
    return freqout, recep, mobil, inert

def power_method(M,K,n=1):
    """Does a thing with the power method
    
    Parameters
    ----------
    x: first signal (x(t))
    f: second signal (f(t))
    dt: the time step of the sampled data.
    n: is the number of points in the fft.
        Recommendedvalue for n is 256 unless zooming.
        Use a higher n for zooming.
        
    Returns
    -------
    lbda: lambda, the first eigenvalue
    phi: the first mode shape
    
    Examples
    --------
    
    """
    shift = 1
    K = K+shift*M
    delta = 10
    X = np.random.normal(0,.1,(np.size(M,1),1))
    X = X/np.linalg.norm(X)
    A = K/M
    oldev=0
    
    for i in range(0,n-1):
        while delta>100*eps:
            Xnew = A*X
            newev = np.linalg.norm(Xnew)
            delta = abs(newev-oldev)
            X = Xnew/np.linalg.norm(Xnew)
            oldev = newev
        lbda(i) = 1/newev-shift
        phi()
        #A = matrdefl(A,newev,X)
        A = A-newev*X*np.transpose(X)
        delta = 10
    return lbda, phi