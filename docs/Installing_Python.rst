.. _installing_python:

Installing Python
_________________

In order to be able to use the Vibration Toolbox you need a working scientific python installation.

The easiest path to this is to install Python via `Anaconda`_. **You must install** Python 3.5 or later for the Vibration Toolbox to work. **Do note install Python 2.7.** 

This proceeds as a normal install on your platform (Mac, Windows, Linux...).

Subsequently you must update some components of Anaconda by using the *conda* command from the *Anaconda Command Prompt*. On Windows, this runs as an actual program.

To ensure the conda package list is as up to date as possible::

  conda update conda

Then update everything else with::

  conda update --all

Then,

.. code-block:: bash

  conda install jupyter

To use `Jupyter`_ (the notebook), launch a terminal on Mac or Linux, or the Anaconda Terminal on Windows (or similar name) and type:

.. code-block:: bash

   jupyter notebook

It's all in the GUI from here. You just need to play around a bit.

.. _github: http://www.github.com
.. _Anaconda: http://continuum.io/downloads
.. _Jupyter: http://www.jupyter.org
