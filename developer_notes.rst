To make distribution.

1) Edit the version number in ``__init__.py``
2) Use the Makefile, ``make release``


To make a wheel file to test before deployment::

  >>> python setup.py bdist_wheel
  >>> python setup.py sdist

To test before release::

  >>> pip install --force-reinstall --upgrade --no-deps dist/vibration_toolbox-0.5b9-py3-none-any.whl

See ``create_distro.rst`` for explicit pypi commands that may not be necessary.
