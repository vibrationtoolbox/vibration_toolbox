This is the folder for building the complete documentation for the Vibration Toolbox. The actual documentation is hosted via `githubpages <http://vibrationtoolbox.github.io/vibration_toolbox/>`_.

The remainder of this document contains notes for developers.

To make the doc, in the `docs` folder type::

  make docs

To serve the pages, move to the `docs/_build/html` folder and type::

  python -m http.server



===============
 Title
===============


Section Title
-------------

Subsection Title
________________

Subsubsection Title
~~~~~~~~~~~~~~~~~~~

Subsubsubsection Title
``````````````````````

Subsubsubsubection Title
''''''''''''''''''''''''

We shouldn't ever need this many.


All documentation will be built with via `sphinx <http://sphinx-doc.org>`_ using ``make html`` in the ``docs`` directory. Updating on ``github`` is documented
in the ``developers.rst`` file at the top of the repository.

Documentation of functions help (uses `autodoc <http://www.sphinx-doc.org/en/stable/ext/autodoc.html>`_):


Module and function help needs to then follow the `numpy convention.
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_


For information on how function docstrings are used, see `numpydoc <https://github.com/numpy/numpy/blob/master/doc/HOWTO_BUILD_DOCS.rst.txt>`_

For plots:
http://matplotlib.org/sampledoc/extensions.html
