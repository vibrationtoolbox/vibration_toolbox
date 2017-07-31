# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 19:02:58 2017

@author: Loranger
"""

import vibration_toolbox as vt
import scipy.io as sio
data = sio.loadmat(vt.__path__[0] + '/data/case2.mat')
TF = data['Hf_chan_2']
f = data['Freq_domain']
z,nf,a = vt.ema.mdof_cf(f,TF,500,1000)
print(z)
print(nf)
print(a)