#!/usr/bin/env python

import matplotlib
import nose
import numpy as np
matplotlib.use('agg')
np.set_printoptions(precision=2, suppress=True)

nose.main(argv=['fake', 'vibration_toolbox','--with-doctest','--doctest-tests','--doctest-options=+ELLIPSIS,+NORMALIZE_WHITESPACE'])

#  - nosetests --with-doctest --doctest-tests --doctest-options=+ELLIPSIS,+NORMALIZE_WHITESPACE
