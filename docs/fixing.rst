Fixing
======

:func:`pylamda.fix`, :func:`pylamda.patch` and ``pylamda fix`` repair the issues that
:doc:`checks` reports, where a repair is possible.

**Rule: fix only what is unambiguous. Every ambiguous fix needs an explicit option.**

.. code-block:: bash

   pylamda fix in.dat -o out.dat --trust id        # normalised layout
   pylamda fix in.dat --in-place --trust name      # patch mode: only the changed lines
   pylamda fix in.dat --dry-run --inverted swap    # show the changes as a diff

.. code-block:: python

   import pylamda

   fixed, log = pylamda.fix("in.dat", trust="id")       # LamdaData and the fixes made
   text, log = pylamda.patch("in.dat", trust="name")    # the patched file as text

Two modes
---------

**Normalised** (``-o OUT``, :func:`pylamda.fix`): the fixed data are written in
pylamda's layout (:func:`pylamda.write`): ``-- Ref:`` partner lines, shortest lossless
numbers. Use it for new or converted files.

**Patch** (``--in-place``, or ``--patch`` with ``-o``; :func:`pylamda.patch`): only the
lines that change are edited, token by token. Column widths, number styles, comments and
line endings stay, so the git diff of a hand-maintained file (ProDiMo's ``data/``)
shows only the fix. ``--dry-run`` shows this diff and writes nothing.

Both modes give the same data. In patch mode, rates that are recomputed (``--inverted
chianti`` or ``balance``) keep the number of digits of the file (at least three
significant digits).

The input file is never overwritten without ``--in-place``.

What is fixed
-------------

In this order:

=========  =========================================================  ==============================
code       fix                                                        needs
=========  =========================================================  ==============================
S05        add ``!`` to a comment line without it                     (default)
S03        extra blocks after the last declared one: count them in    ``--extra-blocks include`` or
           (NPART raised) or comment them out                         ``comment``
S04        drop surplus values (what all readers do)                  (default)
S02        correct the running numbers of the level, line and rate    (default)
           rows (first column; no reader uses it)
C02, C12   the name in the partner text and the ID disagree:          ``--trust id`` or
           ``id`` corrects the name, ``name`` corrects the ID          ``--trust name``
L01        sort the levels by energy, renumber all transitions        (default; ``--no-sort``)
C08        rate rows with E(up) < E(low), see below                   ``--inverted swap``,
                                                                      ``chianti`` or ``balance``
R04, R05   compute missing frequencies and E_up from the energies     (default;
                                                                      ``--no-missing-columns``)
R04, R05   replace values that differ from the energies               ``--recompute``
F01-F03    append an identifier to a reference text                   ``--ref TARGET=ID``
=========  =========================================================  ==============================

There is no default for ``--trust``: wrong IDs and wrong names were both found in real
files. A corrected partner line has the form ``<id> <molecule> + <name> -- Ref:
<reference>``.

Inverted rate rows (C08):

``swap``
   Exchange up and low and keep the rate (the levels were only given in the wrong
   order).
``chianti``
   Exchange and multiply the rate by g(new low) / g(new up). For rates computed from
   a symmetric collision strength Υ that were divided by the wrong statistical weight.
``balance``
   The row holds an excitation rate: it is turned into a de-excitation rate by detailed
   balance, :math:`k_{ul} = k_{lu}\,(g_l/g_u)\,e^{\Delta E/kT}`.

Each fix is recorded: in the returned log (:class:`~pylamda.fix.FixRecord`) and as a
note line at the end of the file:

.. code-block:: text

   ! pylamda 0.1.0 fix: C02 collision block 1: name e -> H2 (trusted ID 1)

After the fixes, ``pylamda fix`` names the fixable issues that remain and the option
each one needs.

Adding identifiers to references
--------------------------------

``--ref TARGET=ID`` (repeatable) appends an identifier to a reference text, so that later
readers and :doc:`references` need no catalogue:

.. code-block:: bash

   pylamda fix co.dat --in-place --ref p-H2=2010ApJ...718.1062Y --ref levels=10.1051/0004-6361:20041729

.. code-block:: python

   data, log = pylamda.add_reference(data, "p-H2", bibcode="2010ApJ...718.1062Y")

- Targets: ``levels`` (the reference on the NLEV line; ``- Ref:`` is added if there is
  none), a partner name, or ``dataset`` (a dataset DOI on line 2; DOIs only).
- The existing text stays; the identifier is appended: a bibcode as SciX link
  (``https://scixplorer.org/abs/<bibcode>``), a DOI as ``doi:<doi>``, an arXiv ID as
  ``arXiv:<id>``.
- An identifier already in the text is not added again.
- Only free text changes, which all readers ignore.
- A bibcode has 19 characters; the last one is the initial of the first author
  (``2007A&A...472.1029F``).

After fixing, ``pylamda fix`` names the reference issues that remain with their
targets:

.. code-block:: text

   pylamda fix: F03 remains (levels, p-H2); it needs --ref levels=ID --ref p-H2=ID or --add-refs

``--add-refs`` does this for all free-text references (references without any
identifier) at once:

.. code-block:: bash

   pylamda fix co.dat --in-place --add-refs                 # from the reference catalogue
   pylamda fix co.dat --in-place --add-refs-from review.json  # from the decisions in a review file

- From the catalogue, only papers found at the same place or for the same molecule are
  written; a match by the text alone is not.
- From a review file (``pylamda bib --search``), only groups with a decision.
- Each paper is appended as SciX link; ``&`` in a bibcode is written as ``%26``.
- The note lines of pylamda (``! pylamda ... fix: ...``, and the same for
  ``add-partner``, ``drop-partner`` and ``scale-partner``) are not counted as references
  by ``check`` and ``bib``.
