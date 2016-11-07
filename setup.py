#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup
import os
import sys

if sys.version_info < (3,5):
    sys.exit('Sorry, Python < 3.5 is not supported.')
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='vibration_toolbox',
      version='0.5b5',
      description='Educational code illustrating fundamentals of vibration for mechanical engineers.',
      author='Joseph C. Slater and Raphael Timbo',
      author_email='joseph.c.slater@gmail.com',
      url='https://github.com/vibrationtoolbox/vibration_toolbox',
      download_url='https://github.com/vibrationtoolbox/vibration_toolbox/archive/0.5b5.tar.gz',
      packages=['vibration_toolbox'],
      long_description = read('readme.rst'),
      keywords=['vibration','mechanical engineering'],
      install_requires=['numpy', 'scipy', 'matplotlib', 'vibration']
      )

# https://docs.python.org/3/distutils/setupscript.html#additional-meta-data
