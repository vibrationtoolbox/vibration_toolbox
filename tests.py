#!/usr/bin/env python

import matplotlib
import nose
matplotlib.use('agg')

nose.main(argv=' --with-doctest --doctest-tests --doctest-options=+ELLIPSIS,+NORMALIZE_WHITESPACE')

#  - nosetests --with-doctest --doctest-tests --doctest-options=+ELLIPSIS,+NORMALIZE_WHITESPACE
