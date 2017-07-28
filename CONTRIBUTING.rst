Reporting bugs
--------------
If you find a bug, please open an issue on the `Github issues tracker <https://github.com/vibrationtoolbox/vibration_toolbox/issues>`_.
Please provide some code that reproduces the error and versions of the packages installed.

Contributing code
-----------------
To contribute code we recommend you follow these steps:

#. Fork the repository on github

#. Create a new branch and add your code/make your modifications, committing to your branch.

#. Set up travis-ci for your branch. This is actually pretty quick and easy:

  #. Go to settings on github for your fork.

  #. Select ``Integration & Services``

  #. Click ``Add service`` and select ``Travis CI``.

4. Clone the repository to your favorite location on your drive where you want to work on it.

#. To work in `developer mode <https://packaging.python.org/distributing/#working-in-development-mode>`_, at the top level directory inside the ``vibration toolbox`` type::

    $ pip install -e .

   This will allow you to edit the code while having it pretend to be installed. Keep in mind, if you have actually installed the ``vibration toolbox`` you may have a conflict. You must uninstall it and install your development version with the command above.

#. If a new function is added
   please provide docstrings following the `Numpy standards for docstrings <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
   The docstrings should contain examples to be tested.

   Specifically note:

   1. Parameters should be listed similarly to:

    |    filename : str
    |    copy : bool
    |    dtype : data-type
    |    iterable : iterable object
    |    shape : int or tuple of int
    |    files : list of str
    |    time : array_like

   2. First line should be inline with the ``"""`` and brief enough to fit on one line.

   3. There must be a blank line after the first line.

   This is not exhaustive. It just highlights some consistent errors made.

#. Run the doctests regularly when you make edits.

   To run the doctests `<pytest https://docs.pytest.org/en/latest/>`_ is needed and can be installed with ``pip install -U pytest``.

   To run the tests from the shell you can `cd` to the project's root directory and type::

     $ pytest


    1. To run the tests from ``pycharm`` you can do: Run -> Edit Configurations -> Add -> python tests -> pytest Then just set the path to the project directory.

    2. To run the tests from ``spyder`` see `spyder-unittest <https://github.com/spyder-ide/spyder-unittest`_.

#. Commit and check `travis-ci <https://travis-ci.org/vibrationtoolbox/vibration_toolbox>`_ tests regularly. Having a great number of changes before a commit can make tracing errors very hard. Make sure you are looking at your branch when assessing whether it's working.

#. Update from the main repository before submitting a pull request. This assures that you can check how your code works with the current repository. If it doesn't work, the pull will (should) be denied.

#. If the tests are passing, make a git pull (in your GitHub app) to assure that your code is up to date with the master branch and that your code has no conflicts with the current base. Doing this regularly ensures that your accumulated edits won't be massively in conflict with the existing code base. After that, push your branch to github and then open a pull request.

Instructions bellow are directed to main developers
===================================================

To make distribution and release
--------------------------------

1) Edit the version number in ``vibration_toolbox/__init__.py``
2) Use the Makefile, ``make release``

The ``conf.py`` file for the documentation pulls the version from ``__init__.py``

To make a distribition (for testing or posting to github)
-----------------------------------------------------------

.. code-block:: bash

  >> make wheel

To test before release
----------------------

Check the Travis CI logs. They are more comprehensive.

To test distribution installabilty
-----------------------------------
Note: these are out of date and saved only for historical reasons.

python setup.py register -r pypitest
python setup.py sdist upload -r pypitest

look at https://testpypi.python.org/pypi

Other information sites
------------------------

`twine notes <https://packaging.python.org/distributing/#working-in-development-mode>`_

https://pypi.python.org/pypi/wheel
