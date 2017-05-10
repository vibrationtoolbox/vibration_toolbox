Running doctests
----------------
To run the doctests `<pytest https://docs.pytest.org/en/latest/>`_ is needed.
To run the tests from the shell you can access the project directory and type:
.. code-block:: shell
    $ pytest

To run the tests from pycharm you can do:
Run -> Edit Configurations -> Add -> python tests -> pytest
Then just set the path to the project directory.

To make distribution.
---------------------

1) Edit the version number in ``vibration_toolbox/__init__.py``
2) Use the Makefile, ``make release``

The ``conf.py`` file for the documentation pulls the version from ``__init__.py``

To make a wheel file to test before deployment::

  >>> make wheel

To test before release::

  >>> pip install --force-reinstall --upgrade --no-deps dist/vibration_toolbox-0.5b9-py3-none-any.whl

See ``create_distro.rst`` for explicit ``pypi`` commands that may not be necessary.

See `twine notes <https://packaging.python.org/distributing/#working-in-development-mode>`_ on modern pypi connectivity.
