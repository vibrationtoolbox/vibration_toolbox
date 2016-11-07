"""
The Vibration Toolbox, Python Edition

Joseph C. Slater and Raphael Timbo

__all__ = ['directory']
"""


"""
With this, `from vtoolbox import *` will import
all toolbox functions in the name space

`import vtoolbox as vtb` will keep them tucked behind `vtb`

`import vtoolbox.sdof as sdof` will tuck the sdoc functions in the `sdof` box.

"""
from .sdof import *
from .mdof import *
from .ema import *
