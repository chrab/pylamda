Installation
============

pylamda needs Python 3.11 or newer. The only required dependency is numpy.

From the repository:

.. code-block:: bash

   git clone https://github.com/chrab/pylamda.git
   cd pylamda
   pip install .

Optional extras:

``plot``
   matplotlib, for plots of comparisons.
``bib``
   requests, to resolve DOIs and bibcodes for the BibTeX export.
``dev``
   pytest, ruff, Sphinx with the extensions of this documentation, and the other two
   extras.

.. code-block:: bash

   pip install ".[plot,bib]"

For development, see :doc:`development`.
