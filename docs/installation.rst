Installation
------------

Easy Installation
_________________

If you aren't familiar at all with Python, please see  `Installing Python <https://github.com/vibrationtoolbox/vibration_toolbox/blob/master/docs/Installing_Python.rst>`_.

Installation is made easy with ``pip`` (or ``pip3``), with releases as we have time while we try
to create a full first release. Much of it works already, but we certainly need
issue reports (on `github <http://github.com/vibrationtoolbox/vibration_toolbox>`_).

To install::

  pip install --user vibration_toolbox

where ``--user`` isn't necessary if you are using a locally installed version of Python such as `Anaconda <https://www.continuum.io/downloads>`_.

To run, I recommend you open a `Jupyter`_ notebook by using ``jupyter notebook`` and then type::

  import vibration_toolbox as vtb

For examples, see the `example ipynb notebooks <https://github.com/vibrationtoolbox/vibration_toolbox/tree/master/docs/tutorial>`_. Some of these have interactive capabilities that are only apparent when you load them with `Jupyter`_ instead of just looking at them on github.

Installation of current development version
___________________________________________

The usage documentation is far behind the current code, while the reference is way ahead of the released code due to the `autodoc <http://www.sphinx-doc.org/en/stable/ext/autodoc.html>`_ capability of `Sphinx <http://www.sphinx-doc.org/en/stable/>`_. Especially as of early 2017, the code is in rapid development. So is the documentation. Releases to `pypi <https://pypi.python.org/pypi>`_ are far behind current status as stopping to deploy would cost more time that it is worth.

If you wish to install the current version of the software please see `CONTRIBUTING.rst <https://github.com/vibrationtoolbox/vibration_toolbox/blob/master/CONTRIBUTING.rst>`_.

That should be it. Please note issues on the `issues tab <https://github.com/vibrationtoolbox/vibration_toolbox>`_ on github.

.. _Jupyter: jupyter.org
