Reporting bugs
--------------
If you find a bug, please open an issue on the `Github issues tracker <https://github.com/vibrationtoolbox/vibration_toolbox/issues>`_.
Please provide some code that reproduces the error and versions of the packages installed.

Contributing code
-----------------
To contribute code we recommend you follow these steps:

1. Clone the repository:

.. code-block:: bash

    >> git clone https://github.com/vibrationtoolbox/vibration_toolbox

2. Create a new branch and add your code. If a new function is added
please provide docstrings following the
`Numpy standards for docstrings <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
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

3. Run the doctests regularly when you make edits.

To run the doctests `<pytest https://docs.pytest.org/en/latest/>`_ is needed.
To run the tests from the shell you can access the project directory and type:

    $ pytest

To run the tests from pycharm you can do:
Run -> Edit Configurations -> Add -> python tests -> pytest
Then just set the path to the project directory.

4. Commit and check `travis-ci <https://travis-ci.org/vibrationtoolbox/vibration_toolbox>`_ tests regularly. Having a great number of changes before a commit can make tracing errors very hard. Make sure you are looking at your branch when assessing whether it's working. 

5. Update from the main repository before submitting a pull request. This assures that you can check how your code works with the current repository. If it doesn't work, the pull will (should) be denied. 

6. If the tests are passing, make a git pull to assure that your code is up to date with the master branch and that the code has no conflicts. After that, push your branch to github and then open a pull request.


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

.. code-block:: bash

  >> pip install --force-reinstall --upgrade --no-deps dist/vibration_toolbox-0.5b9-py3-none-any.whl

See `notes <https://packaging.python.org/distributing/#working-in-development-mode>`_ on working in development mode.

To test distribution installabilty
-----------------------------------
Note: these are out of date. 

python setup.py register -r pypitest
python setup.py sdist upload -r pypitest

look at https://testpypi.python.org/pypi

Other information sites
------------------------

`twine notes <https://packaging.python.org/distributing/#working-in-development-mode>`_ 

https://pypi.python.org/pypi/wheel
