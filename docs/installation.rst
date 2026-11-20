Installation
============

pylamda needs Python 3.11 or newer. The only required dependency is numpy.

From the repository:

.. code-block:: bash

   git clone https://github.com/chrab/pylamda.git
   cd pylamda
   pip install .

The dependencies are numpy and requests (requests resolves DOIs and bibcodes for the
BibTeX export).

Optional extras:

``plot``
   matplotlib, for plots of comparisons.
``dev``
   pytest, pytest-cov, ruff, Sphinx with the extensions of this documentation, and the
   ``plot`` extra.

.. code-block:: bash

   pip install ".[plot]"

For development, see :doc:`development`.
