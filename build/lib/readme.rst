=============================================
 The Engineering Vibration Toolbox for Python
=============================================

.. .. include:: <isonum.txt>
.. image:: https://badge.fury.io/py/vibration-toolbox.png/
    :target: http://badge.fury.io/py/vibration-toolbox

.. image:: https://travis-ci.org/vibrationtoolbox/vibration_toolbox.svg?branch=master
    :target: https://travis-ci.org/vibrationtoolbox/vibration_toolbox

.. image:: https://zenodo.org/badge/39572419.svg
    :target: https://zenodo.org/badge/latestdoi/39572419

.. image:: https://mybinder.org/badge.svg
    :target: https://mybinder.org/v2/gh/vibrationtoolbox/vibration_toolbox/binder
   
.. image:: https://img.shields.io/badge/Say%20Thanks-!-1EAEDB.svg 
    :target: https://saythanks.io/to/josephcslater

.. image:: https://img.shields.io/badge/patreon-donate-yellow.svg
   :target: https://www.patreon.com/josephcslater

.. .. image:: https://img.shields.io/pypi/v/vibration_toolbox.svg
    :target: https://img.shields.io/pypi/v/vibration_toolbox

.. #image:: https://coveralls.io/repos/vibrationtoolbox/vibration_toolbox/badge.png?branch=master
..  #:target: https://coveralls.io/r/vibrationtoolbox/vibration_toolbox

.. image:: http://pepy.tech/badge/vibration-toolbox
   :target: http://pepy.tech/project/vibration-toolbox
   :alt: PyPi Download stats


Joseph C. Slater and Raphael Timbó
----------------------------------

Welcome to the `Vibration Toolbox <http://vibrationtoolbox.github.io/vibration_toolbox/>`_.
This `Python <http://python.org>`_ version is a completely new design build for modern education. This is an *educational* set of codes intended primarily for
demonstration of vibration analysis and phenomenon. You may find them useful for application, but that isn't the intent of this toolbox. If you have professional-level needs please `contact the authors <mailto:joseph.c.slater@gmail.com>`_.

Full `documentation is available <http://vibrationtoolbox.github.io/vibration_toolbox/>`_, but please excuse that it is still under development. Such documentation has never existed for the other ports of the toolbox so this is taking some time. We don't need feedback at this time, but we will take assistance in improving documentation and code. *Please* clone the repository and support use by submitting pull requests fixing typos and clarifying documentation.


Try now!
--------

You won't get everything, but you can try parts of the toolbox immediately on  and `tutorial on mybinder.org <https://mybinder.org/v2/gh/vibrationtoolbox/vibration_toolbox/binder>`_, right in your browser window.


Installation
------------

If you aren't familiar at all with Python, please see  `Installing Python <https://github.com/vibrationtoolbox/vibration_toolbox/blob/master/docs/Installing_Python.rst>`_.

Installation is made easy with ``pip`` (or ``pip3``), with releases as we have time while we try
to create a full first release. Much of it works already, but we certainly need
issue reports (on `github <http://github.com/vibrationtoolbox/vibration_toolbox>`_).

To install type::

  pip install --user vibration_toolbox

at your command prompt **(not the python prompt)** where ``--user`` isn't necessary if you are using a locally installed version of Python such as `Anaconda <https://www.continuum.io/downloads>`_.

To run, I recommend you open a `Jupyter <https://jupyter.org>`_ notebook by using ``jupyter notebook`` at your command prompt/terminal prompt/Anaconda prompt and then type::

  import vibration_toolbox as vtb

For examples, see the `JupyterNotebooks folder <https://github.com/vibrationtoolbox/vibration_toolbox/tree/master/docs/tutorial>`_. Some of these have interactive capabilities that are only apparent when you run them yourself instead of just looking at them on GitHub. Unfortunately our organization of these still leaves a little to be desired. Help accepted!

Installation of current code/contributing
_________________________________________

The usage documentation is far behind the current code, while the reference is way ahead of the released code due to the `autodoc <http://www.sphinx-doc.org/en/stable/ext/autodoc.html>`_ capability of `Sphinx <http://www.sphinx-doc.org/en/stable/>`_. Especially as of 2017, the code is still in rapid development. So is the documentation. Releases to `pypi <https://pypi.python.org/pypi>`_.

If you wish to install the current version of the software, and especially contribute, please follow the instructions in `Contributing.rst <https://github.com/vibrationtoolbox/vibration_toolbox/blob/master/CONTRIBUTING.rst>`_

That should be it. Please note issues on the `issues tab <https://github.com/vibrationtoolbox/vibration_toolbox>`_ on GitHub.
