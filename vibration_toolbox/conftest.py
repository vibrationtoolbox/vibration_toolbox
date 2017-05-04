import numpy as np
import matplotlib as mpl
import pytest

np.set_printoptions(precision=2, suppress=True)
print('running conftest')


@pytest.fixture(autouse=True, scope='session')
def set_backend():
    print('backend has been set to agg')
    mpl.use('agg')
