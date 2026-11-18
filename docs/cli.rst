Command line
============

``pylamda`` (or ``python -m pylamda``) has one command so far, ``check``.

pylamda check
-------------

.. code-block:: text

   pylamda check [options] FILE [FILE ...]

Checks LAMDA files against the standard and prints the issues (see :doc:`checks`):

.. code-block:: text

   $ pylamda check data/CH+_Lambda.dat
   data/CH+_Lambda.dat:6: F03 info: the levels and lines have a reference only as free text (no DOI, bibcode or link)
   data/CH+_Lambda.dat:8: L04 info: 16 of 16 levels without quanta
   data/CH+_Lambda.dat:33: R04 warning: frequency differs from the level energies beyond rounding and rtol 0.0001 in 9 lines; ...
   data/CH+_Lambda.dat:45: C02 warning: partner line '1 CH+ + e ! Lim et al 1999 ...': ID 1 is H2, the name says e (read as H2; ProDiMo would read e)
   data/CH+_Lambda.dat:45: C04 warning: H2 (collision block 1) together with p-H2: some codes count H2 collisions twice
   data/CH+_Lambda.dat:46: S05 error: comment line without '!' before the number of collisional transitions of collision block 1: 'NUMBER OF COLL TRANS'
   ...
   Found 12 issues (2 errors, 4 warnings, 6 infos) in 1 file.

Options:

``--select CODES``
   Report only these codes: comma-separated prefixes (``C``, ``C0``, ``C02``); can be
   repeated.
``--ignore CODES``
   Do not report these codes.
``--error-on CODES``
   Report these codes as errors, e.g. ``--select C02 --error-on C02`` in a CI job for
   ProDiMo's files.
``--no-physics``
   Leave out the physics heuristics (the same as ``--ignore P``).
``--rtol X``
   Relative tolerance of R04 and R05, on top of the rounding (default 1e-4).
``--format text|json``
   Output format. JSON is a list with one entry per file: ``path``, ``readable``,
   ``error`` and ``issues``.
``--statistics``
   Show the number of issues and files per code instead of each issue.
``--list-codes``
   List all check codes with their severity and description.

Exit codes:

- 0: no issue of severity ``error``;
- 1: at least one error;
- 2: at least one file cannot be read (or opened).
