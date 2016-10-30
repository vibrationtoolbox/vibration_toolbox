"""
I believe we need to keep this project simple for users.
Remember the intent is to a demonstrate/teach vibration.
This will not be a massively comprehensive set of tools,
so  `__all__` should contain any submodules should also be imported
when a user executes "from vtoolbox import *"

As there are no submodules, this is a placeholder in case such ever happens.

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
