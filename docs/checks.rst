Checks
======

:func:`pylamda.verify` and ``pylamda check`` report every deviation from the `LAMDA
standard <https://home.strw.leidenuniv.nl/~moldata/molformat.html>`_ as an issue with a
code, a severity and a message, as ruff does for Python code.

.. code-block:: python

   import pylamda

   for issue in pylamda.verify("co.dat"):
       print(issue)
   # line 95: C04 warning: H2 (collision block 1) together with p-H2/o-H2: some codes count H2 collisions twice

Each issue has:

- ``code`` and ``severity`` (``error``, ``warning`` or ``info``);
- ``line``: the line in the file (for an issue about several rows, the first one);
- ``count``: the number of rows, levels or lines affected. A check of many rows gives
  one issue per block or section, with the first example, not one issue per row;
- ``fixable``: whether ``pylamda fix`` (step 5 of the design) can remove it.

Options:

- ``select`` / ``ignore`` take code prefixes: ``"C"`` (all C codes), ``"C0"``, ``"C02"``.
- ``physics``: the P checks are off by default; ``physics=True`` runs them,
  ``select=["P"]`` (or single P codes) too, and ``physics=False`` never runs them.
- ``error_on`` raises codes to ``error``, for example ``["C02"]`` for ProDiMo's files,
  which must have no C02 issue.

S: structure
------------

The reader (:func:`pylamda.read`) finds these while it reads the file.

=========  ========  =========================================================================
code       severity  issue
=========  ========  =========================================================================
``S01``    error     file shorter than its counts require (the reader stops)
``S02``    error     level or transition numbers not 1..N in order
``S03``    warning   data after the last declared block (more blocks than NPART)
``S04``    warning   more values on a line than declared (NTEMP)
``S05``    error     text line without ``!`` between blocks (read as data by Fortran)
``S06``    error     number that Fortran cannot read (``1,0e-10``), or a partner line without the ID
``S07``    error     level or radiative-transition row with fewer values than required
=========  ========  =========================================================================

L: levels
---------

=========  ===============  ==============================================================
code       severity         issue
=========  ===============  ==============================================================
``L01``    warning          energies not in increasing order
``L02``    error / info     weight g <= 0 or not finite (error); not an integer (info)
``L03``    warning          two levels with the same E, g and quanta
``L04``    info             quanta empty or not unique
``L05``    info             the lowest level has E ≠ 0
=========  ===============  ==============================================================

R: radiative transitions
------------------------

=========  ========  =========================================================================
code       severity  issue
=========  ========  =========================================================================
``R01``    error     upper or lower level outside 1..NLEV, or upper = lower
``R02``    error     E(up) <= E(low)
``R03``    error     Einstein A <= 0 or not finite
``R04``    warning   frequency missing, or not (E_up − E_low) · c
``R05``    warning   E_up [K] missing, or not E(up) · hc/k
``R06``    error     the same (up, low) twice
=========  ========  =========================================================================

R04 and R05 compare the columns with the values computed from the level energies. The
tolerance is ``rtol`` (default 1e-4, relative) plus the rounding of the numbers as
written in the file: half a unit of the last digit of the column value and of the
energies. A frequency of a hyperfine line computed from energies with six decimals is
therefore not reported. An energy of exactly 0 counts as exact.

C: collisions
-------------

=========  ========  =========================================================================
code       severity  issue
=========  ========  =========================================================================
``C01``    error     partner ID not in 1..7
``C02``    warning   the partner text names another standard partner than the ID
``C03``    error     the same partner ID twice
``C04``    warning   H2 together with p-H2 or o-H2 (counted twice by some codes)
``C05``    error     temperatures not strictly increasing, or <= 0
``C06``    error     rate < 0 or not finite
``C07``    info      rate row that is 0 at all temperatures
``C08``    warning   E(up) < E(low) in a rate row (wrong order, or an excitation rate)
``C09``    error     the same (up, low) twice in one block
``C10``    error     upper or lower level outside 1..NLEV, or upper = lower
``C11``    warning   the partner text names a partner outside the standard (read by its ID)
``C12``    warning   the partner text names He, the ID is H2, and the text says "scaled"
``C13``    info      the partner text says the rates are scaled
=========  ========  =========================================================================

pylamda reads every block by its partner ID, as RADEX does. ProDiMo reads the name in
the partner text first. If the two disagree (C02), the codes use different rates. EMAA's
partners p-H2O, o-H2O and CO are not in the standard and reuse the IDs 2, 3 and 5 (C11,
and C03 if the file also has p-H2 or o-H2).

C12 is a special case of C02. Some files give rates for H2 that were scaled from He
rates (by about 1.37, the ratio of the reduced masses), but keep the He name:
``1 SO2 + He ! rates (Green 1995 ApJS) scaled *1.39``. The ID (1 = H2) is probably
right and the name is wrong. C12 reports such a line instead of C02 if the name is He,
the ID is 1, 2 or 3, and the text says "scaled", "scaling" or gives a factor such as
``*1.39`` or ``x1.4``.

C13 is an info: the partner text says "scaled", "scaling" or "from H2" (also "from
p-H2", ...). The message quotes the matched text, with a factor that follows it, e.g.
``"scaled (*1.37)"``; check the reference for the details. :func:`pylamda.scale_partner`
can scale rates of one partner to another one, or back (:doc:`combining`).
C13 was the physics check P02 until 2026-11-23; it is a check of the partner text, so it
runs by default.

F: references
-------------

=========  ========  =========================================================================
code       severity  issue
=========  ========  =========================================================================
``F01``    warning   no reference for the levels and lines
``F02``    warning   partner without a reference
``F03``    info      reference without a DOI, bibcode or link (free text only)
=========  ========  =========================================================================

The reference of the levels and lines can be on the NLEV line (EMAA: ``Ref:``), on line
2, in a comment line or in the notes. F01 reports a file where none of these holds a
DOI, a link, a bibcode, an arXiv ID, "et al.", a year or a database name (CDMS, JPL,
NIST, …). F03 reports references with none of DOI, link, bibcode or arXiv ID that are
not in the reference catalogue (:doc:`references`; ``--no-catalogue`` leaves the
catalogue out).

A check that DOIs and bibcodes resolve (F04, ``--online``) follows.

P: physics heuristics
---------------------

The P checks are off by default (decided 2026-11-23): they can give false alarms, and
P01 fits every rate row, which is slow on large files (a 80 MB rovibrational file: about
6 s). Run them with ``physics=True`` (``pylamda check --physics``) or by selecting them
(``select=["P"]``, ``--select P``). Their severity is ``info``.

=========  ========  =========================================================================
code       severity  issue
=========  ========  =========================================================================
``P01``    info      rate row that behaves like an excitation rate
``P03``    info      rate outside the plausible range
``P04``    info      mass differs from the formula mass of the molecule
=========  ========  =========================================================================

``P01``
   An excitation rate falls like exp(−ΔE/kT) at T ≪ ΔE/k, so ln k against 1/T has the
   slope −ΔE/k. P01 fits this slope over the temperatures below ΔE/2k (at least two)
   and reports rows whose slope is between 2/3 and 3/2 of −ΔE/k. A row that is labelled
   with the wrong levels has another ΔE and may not be found.
``P03``
   A rate above 1e-8 cm³ s⁻¹ for a neutral partner, 1e-7 for H+, or 1e-4 for
   electrons (:data:`pylamda.verify.P03_MAX_RATE`).
``P04``
   The molecule name is read as a formula (``13CO``, ``C18O``, ``o-H3O+``, ``A-CH3OH``)
   and its mass compared with the file's mass (tolerance 0.5 Da). Names that do not
   read as a formula are not checked.
