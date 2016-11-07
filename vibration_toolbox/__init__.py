"""
The Vibration Toolbox, Python Edition

Joseph C. Slater and Raphael Timbo

__all__ = ['directory']
"""


"""
With this, `from vibration_toolbox import *` will import
all toolbox functions in the name space

`import vibration_toolbox as vt` will keep them tucked behind `vt`

`import vibration_toolbox.sdof as sdof` will tuck the sdof functions in the `sdof` name space.

"""
from .sdof import *
from .mdof import *
from .ema import *
