Command line
============

``pylamda`` (or ``python -m pylamda``) has seven commands: ``check``, ``compare``, ``bib``, ``fix``, ``add-partner``, ``drop-partner`` and ``scale-partner``.

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
   data/CH+_Lambda.dat:45: C02 warning: partner "CH+ + e" (ID 1 = H2): the name says e (read as H2, ProDiMo reads e)
   data/CH+_Lambda.dat:45: C04 warning: partner "CH+ + e" (ID 1 = H2): H2 together with p-H2: some codes count H2 collisions twice
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
``--physics``
   Also run the physics heuristics P01-P04 (off by default; slow on large files).
   ``--select P`` (or single P codes) runs them too.
``--no-physics``
   Never run the physics heuristics, even with ``--select P``.
``--rtol X``
   Relative tolerance of R04 and R05, on top of the rounding (default 1e-4).
``--format text|json``
   Output format. JSON is a list with one entry per file: ``path``, ``readable``,
   ``error`` and ``issues``.
``--statistics``
   Show the number of issues and files per code instead of each issue.
``--list-codes``
   List all check codes with their severity and description.
``--no-catalogue``
   F03: do not use the reference catalogue.
``--skip-non-lamda``
   Leave out files that are not LAMDA files (see below); the summary gives their number.

Exit codes:

- 0: no issue of severity ``error``;
- 1: at least one error;
- 2: at least one file cannot be read (or opened), or is not a LAMDA file.

In text output, the issues of each file are printed as soon as the file is checked; the
summary line follows at the end. JSON and ``--statistics`` are printed after the last
file.

A file whose header does not read is reported with one line, not with its structure
issues, and the rest of it is not read (see :doc:`usage`, "Not a LAMDA file"):

.. code-block:: text

   $ pylamda check data/*.dat
   ...
   data/ATOM_Al.dat: not a LAMDA file: line 1: the molecule name '30' is a number
   data/DustChem.dat: not a LAMDA file: line 2: cannot read the molecular weight from '============'
   Found 633 issues (155 errors, 237 warnings, 241 infos) in 184 files; 70 not LAMDA files; 1 not readable.

With ``--format json``, such a file has ``"readable": false`` and ``"error": "not a LAMDA
file: <reason>"``.

``--skip-non-lamda`` leaves such files out: they give no output line and no JSON entry,
they are not counted as checked files, and they do not change the exit code. The summary
gives their number, so that nothing disappears without notice. Files that are LAMDA
files but cannot be read after the header ("not readable"), and files that cannot be
opened, are still reported (exit code 2):

.. code-block:: text

   $ pylamda check --skip-non-lamda data/*.dat
   ...
   Found 379 issues (60 errors, 153 warnings, 166 infos) in 105 files; 71 non-LAMDA files skipped.

pylamda compare
---------------

.. code-block:: text

   pylamda compare [options] A B

Compares two LAMDA files and prints what differs from A to B (see :doc:`comparing`), or
``no differences``.

Options:

``--match auto|quanta|energy|index``
   How levels are matched (default ``auto``).
``--rtol X``
   Relative change that counts in the exact comparison (default 1e-6).
``--atol X``
   Energy tolerance [cm⁻¹] for matching levels by energy (default 1e-3).
``--limit X``
   Relative change that counts on the union grid, per level and per E_up range
   (default 0.1).
``--format text|json``
   Output format.
``--max-items N``
   How many items of a list are named before "..." (default 5).

Exit codes, as ``diff``:

- 0: no differences;
- 1: differences;
- 2: a file cannot be read (or opened).

pylamda bib
-----------

.. code-block:: text

   pylamda bib [options] FILE [FILE ...]

Writes the references of LAMDA files as BibTeX (see :doc:`references`). Warnings and
the references without an identifier go to standard error.

Options:

``-o FILE``, ``--output FILE``
   The ``.bib`` file (default: standard output).
``--online``
   Resolve the identifiers with SciX (needs the API token) and doi.org.
``--scix-token-file FILE``
   File with the SciX API token (default: ``SCIX_API_TOKEN``, then
   ``~/.config/pylamda/scix_token``).
``--no-database-refs``
   Leave out the LAMDA and EMAA papers.
``--keys FILE``
   JSON file identifier → key; read and updated, so that keys stay the same.
``--no-cache``
   Do not use the cache of resolved references.
``--no-catalogue``
   Do not look up free-text references in the reference catalogue.
``--search``
   Search SciX for free-text references not in the catalogue (needs the token). The
   candidates go into the review file only.
``--review FILE``
   The review file of ``--search`` (default ``references_review.json``); decisions in it
   are kept.
``--no-hints``
   ``--search``: do not read the LAMDA/EMAA molecule pages.

Exit codes:

- 0: done;
- 1: online resolution failed (e.g. SciX rejected the token);
- 2: a file cannot be read (the others are still exported).

pylamda fix
-----------

.. code-block:: text

   pylamda fix FILE (-o OUT | --in-place | --dry-run) [options]

Fixes a LAMDA file (see :doc:`fixing`). The fixes made go to standard error, followed
by the fixable issues that remain and the option each one needs.

Output (one is required):

``-o OUT``
   Write the fixed file to OUT, in the normalised layout (``--patch``: patch mode).
``--in-place``
   Change FILE itself, in patch mode (``--normalise``: normalised layout).
``--dry-run``
   Show the changes of patch mode as a unified diff; write nothing.

Fixes:

``--trust id|name``
   C02/C12: correct the partner name (``id``) or the ID (``name``).
``--inverted swap|chianti|balance``
   C08: rate rows with E(up) < E(low).
``--extra-blocks include|comment``
   S03: count extra blocks in, or comment them out.
``--recompute``
   R04/R05: also replace frequencies and E_up values that differ from the energies.
``--ref TARGET=ID``
   Append a DOI, bibcode or arXiv ID to a reference (TARGET: ``levels``, ``dataset`` or
   a partner); repeatable.
``--add-refs``
   Append the SciX links of all free-text references from the reference catalogue.
``--add-refs-from REVIEW``
   The same, from the decisions in a review file (``pylamda bib --search``). Use only
   one of ``--add-refs`` and ``--add-refs-from``.
``--no-sort``
   L01: do not sort the levels.
``--no-missing-columns``
   R04/R05: do not compute missing frequencies and E_up values.

Exit codes:

- 0: done;
- 2: the file cannot be read, or a fix cannot be applied.

pylamda add-partner, pylamda drop-partner
-----------------------------------------

.. code-block:: text

   pylamda add-partner TARGET SOURCE PARTNER -o OUT [--match M] [--atol X] [--rtol-e X]
                       [--replace] [--block N] [--source-label TEXT]
   pylamda drop-partner FILE PARTNER -o OUT [--block N]

See :doc:`combining`. Both write the normalised layout to a new file (``OUT`` must not be
an input file) and print the note line. Exit codes: 0 done, 2 an error (file not
readable, partner missing or ambiguous, levels that do not fit).

pylamda scale-partner
---------------------

.. code-block:: text

   pylamda scale-partner FILE PARTNER FACTOR (-o OUT | --in-place) [--to NAME] [--keep]
                         [--block N] [--replace] [--reference TEXT] [--digits N] [--normalise]

Multiply the collision rates of one partner by ``FACTOR`` (a number or a fraction
``a/b``, e.g. ``1/1.37``). See :doc:`combining`. Keeps the layout of ``FILE``: only the
partner line and the rates of the block change, and a note line is added at the end.
``--normalise``, ``--keep`` and ``--replace`` write the normalised layout. Writes a new
file (``-o OUT``, not ``FILE``) or changes ``FILE`` (``--in-place``), prints the warnings
and the note line, and then the C issues of the result (``pylamda check --select C``). Exit codes: 0 done, 2 an error
(file not readable, partner missing or ambiguous, factor not > 0, partner exists
already).
