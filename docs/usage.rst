Reading and writing
===================

Reading
-------

:func:`pylamda.read` reads a file by position, as RADEX and ProDiMo do (Fortran
list-directed reads):

- Lines whose first non-blank character is ``!`` are comments. Blank lines are skipped.
- The declared counts decide. Header text is not used.
- Values that do not fit on one line (temperatures, rate rows) continue on the next
  line.
- Surplus values on a line are ignored.
- Everything after the last declared block is kept as ``notes``.

.. code-block:: python

   import pylamda

   data = pylamda.read("co.dat")
   data.molecule, data.n_levels, data.n_lines      # ('CO', 41, 40)
   data.level_energy[:3]                           # [cm-1]
   data.partners                                   # ['p-H2', 'o-H2']
   ph2 = data.partner("para-H2")                   # any known spelling
   ph2.temperatures, ph2.rates.shape               # [K], (ncol, ntemp)
   ph2.rate(100.0)                                 # [cm3 s-1], linear in T

The data model is :class:`~pylamda.model.LamdaData`, with one
:class:`~pylamda.model.CollisionSet` per collision block. Levels and lines are stored
as one numpy array per column, with 1-based level numbers as in the file.

Collision partners
------------------

The standard defines seven partners by ID: 1 = H2, 2 = p-H2, 3 = o-H2, 4 = e, 5 = H,
6 = He, 7 = H+. The rest of the partner line is free text.

pylamda identifies each partner **by its ID only**, as RADEX does.
:attr:`CollisionSet.partner <pylamda.model.CollisionSet.partner>` is the standard name
of the ID. The partner line as written is kept in
:attr:`~pylamda.model.CollisionSet.partner_text`, and the name in it (if the text
names one) in :attr:`~pylamda.model.CollisionSet.text_partner`:

.. code-block:: python

   c = data.collisions[0]          # partner line "1 CO + p-H2 ! Yang et al. 2010"
   c.partner_id, c.partner         # (1, 'H2'): the ID decides
   c.text_partner                  # 'p-H2': the name in the text disagrees
   c.reference                     # 'Yang et al. 2010'

Some codes (ProDiMo) read the name first. Such a file is then read differently by
different codes. The check ``C02`` of the coming ``verify`` reports it.

Issues
------

The reader collects every deviation from the standard as an
:class:`~pylamda.issues.Issue` with a code, a severity, the line number and a message:

.. code-block:: python

   for issue in pylamda.read("s-h.dat").issues:
       print(issue)
   # line 25: S04 warning: collision block 1: more than the 20 declared temperatures; the rest is ignored
   # line 27: S04 warning: collision block 1: 3 rate rows have more than 20 values; the rest is ignored

The reader finds the structure issues:

======  ========  ==================================================================
code    severity  issue
======  ========  ==================================================================
S01     error     file shorter than its counts require (the reader stops)
S02     error     level or transition numbers not 1..N in order
S03     warning   data after the last declared block (more blocks than NPART)
S04     warning   more values on a line than declared (NTEMP)
S05     error     text line without ``!`` between blocks (read as data by Fortran)
S06     error     number that Fortran cannot read
S07     error     level or radiative-transition row with fewer values than required
======  ========  ==================================================================

With ``strict=True``, every structure issue raises
:class:`~pylamda.issues.LamdaReadError`, for example in CI. Without it, only S01 and an
unreadable count stop the reader.

The :class:`~pylamda.model.SourceMap` in ``data.source_map`` holds the line number of
every level, line, partner line and rate row, and the raw lines of the file.

Writing
-------

:func:`pylamda.write` writes one normalised layout: the layout of LAMDA's format
description, with the references in EMAA's ``-- Ref:`` style, and every number as the
shortest text that reads back unchanged.

.. code-block:: python

   pylamda.write(data, "co_normalised.dat")

- The partner line becomes ``<id> <molecule> + <name> -- Ref: <reference>``. The ID is
  written as it is. The name is the one in the original text if there is one, else the
  standard name of the ID. The writer does not decide between a name and an ID that
  disagree.
- Notes and data after the last block are written as comment lines.
- Read → write → read gives the same data.

Checking
--------

:func:`pylamda.verify` runs all checks: the reader's structure issues and the checks of
levels, lines, collisions, references and physics. See :doc:`checks` for the codes and
:doc:`cli` for ``pylamda check``.

.. code-block:: python

   issues = pylamda.verify("co.dat")
   issues = pylamda.verify(data, physics=False, ignore=["F03"])
