"""
The Vibration Toolbox, Python Edition

Joseph C. Slater and Raphael Timbo
"""


"""
With this, `from vibration_toolbox import *` will import
all toolbox functions in the name space

`import vibration_toolbox as vt` will keep them tucked behind `vt`

`import vibration_toolbox.sdof as sdof` will tuck the sdof functions in the
 `sdof` name space.

"""

__title__ = 'vibration_toolbox'
__version__ = '0.5b9'
__author__ = u'Joseph C. Slater and Raphael Timbó'
__license__ = 'MIT'
__copyright__ = 'Copyright 1991-2017 Joseph C. Slater'
__all__ = ['sdof', 'mdof', 'ema', 'vibsystem', 'continuous_systems',
           '__version__']

import scipy as sp
import matplotlib.pyplot as plt

from .sdof import *
from .mdof import *
from .ema import *
from .vibsystem import *
from .continuous_systems import *
