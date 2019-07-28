from .continuous_systems import *
from .vibesystem import *
from .ema import *
from .mdof import *
from .sdof import *
import matplotlib as mpl
import sys
"""The Vibration Toolbox, Python Edition.

Joseph C. Slater and Raphael Timbó

`import vibration_toolbox as vtb` will keep them inside the `vtb` namespace

`import vibration_toolbox.sdof as sdof` will keep the sdof functions in the
 `sdof` namespace.
"""

__title__ = 'vibration_toolbox'
# version may have no more then numerical digits after decimal point.
# 1.11 is actually a higher release than 1.2 (confusing)
__version__ = '0.6.7'
__author__ = u'Joseph C. Slater and Raphael Timbó'
__license__ = 'MIT'
__copyright__ = 'Copyright 1991-2019 Joseph C. Slater'
__all__ = ['sdof', 'mdof', 'ema', 'vibesystem', 'continuous_systems',
           '__version__']

"""
If the __all__ above is commented out, this code will then execute to
completion, as the default behaviour of import * is to import all symbols
that do not begin with an underscore, from the given namespace.

Reference:
https://docs.python.org/3.5/tutorial/modules.html#importing-from-a-package
"""


if 'pytest' in sys.argv[0]:
    # print('Setting backend to agg to run tests')
    mpl.use('agg')


# print options were change inside modules to produce better
# outputs at examples. Here we set the print options to the
# default values after importing the modules to avoid changing
# np default print options when importing the toolbox.
np.set_printoptions(edgeitems=3, infstr='inf', linewidth=75,
                    nanstr='nan', precision=8, suppress=False,
                    threshold=1000, formatter=None)
