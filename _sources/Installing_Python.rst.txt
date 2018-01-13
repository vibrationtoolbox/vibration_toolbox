.. _installing_python:

Installing Python
_________________

In order to be able to use the Vibration Toolbox you need a working `Scientific Python`_ installation.

The easiest path to this is to install Python via `Anaconda`_ or `Enthought Canopy`_. **You must install** Python 3.5 or later for the Vibration Toolbox to work. I prefer the `Anaconda`_ distribution myself, but many organizations prefer `Enthought Canopy`_. **Do not install Python 2.7.**  The Vibration Toolbox requires a Python 3.5 or later.

This proceeds as a normal install on your platform (Mac, Windows, Linux...).

Subsequently you must update some components of Anaconda by using the *conda* command from the *Anaconda Command Prompt*. On Windows, this runs as an actual program.

To ensure the conda package list is as up to date as possible::

  conda update conda

Then update everything else with::

  conda update --all

To use the `Jupyter`_ (the notebook), launch a terminal on Mac or Linux, or the Anaconda Terminal on Windows (or similar name for the `Enthought Canopy`_ distribution of Scientific Python) and type:

.. code-block:: bash

   jupyter notebook

A Matlab_-like experience is to use Spyder, which is also included. To run Spyder, at your terminal on Mac or Linux, or the Anaconda Terminal on Windows, type::

.. code-block:: bash

   spyder

It's all in the GUI from here. You just need to play around a bit.

.. _github: http://www.github.com
.. _Anaconda: http://continuum.io/downloads
.. _Jupyter: http://www.jupyter.org
.. _`Enthought Canopy`: https://store.enthought.com/downloads/
.. _`Scientific Python`: https://www.scipy.org
.. _`Matlab`: http://www.mathworks.com
