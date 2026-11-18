pylamda
=======

**A Python toolkit for files in the LAMDA format.**

The `Leiden Atomic and Molecular Database (LAMDA)
<https://home.strw.leidenuniv.nl/~moldata/>`_ defines a text format for energy levels,
Einstein A coefficients and collision rates. LAMDA, `EMAA <https://emaa.osug.fr>`_, and
codes such as RADEX and ProDiMo use it. ``pylamda`` reads and writes these files, and
will verify, compare, combine and fix them, and export their references as BibTeX.

pylamda follows the `LAMDA standard
<https://home.strw.leidenuniv.nl/~moldata/molformat.html>`_ strictly:

- collision partners are identified by their ID (1–7) only;
- deviations are reported as issues, never silently interpreted.

.. code-block:: python

   import pylamda

   data = pylamda.read("co.dat")
   for issue in data.issues:          # structure issues found by the reader
       print(issue)
   pylamda.write(data, "co_normalised.dat")

.. note::

   pylamda is at an early stage. Reading and writing work; ``verify``, ``compare``,
   ``fix`` and ``bib`` follow (see the :doc:`design`, section 14).

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   usage

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/index
   design

.. toctree::
   :maxdepth: 1
   :caption: Development

   development
