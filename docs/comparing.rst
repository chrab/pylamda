Comparing
=========

:func:`pylamda.compare` reports what differs from one data set to another: header,
levels, radiative transitions and, per collision partner, the rates.

.. code-block:: python

   import pylamda

   diff = pylamda.compare("co_old.dat", "co_new.dat")
   if diff:                          # True if anything differs
       print("\n".join(diff.summary()))
   diff.as_dict()                    # the same as a JSON-compatible dict

.. code-block:: bash

   pylamda compare co_old.dat co_new.dat
   pylamda compare co_old.dat co_new.dat --format json

Matching levels
---------------

Levels are matched first (:func:`pylamda.match_levels`). Lines and rate rows are then
matched through the level map, so that an inserted level shows up as one added level,
not as a change of all later entries.

``auto`` (default)
   By the quantum-number text if it is unique in both data sets and at least half of
   the levels match that way. The text is compared in lower case, with blanks, ``_``,
   ``,`` and ``;`` taken as one separator. Levels left over are then matched by energy
   and g (method ``quanta+energy``). Otherwise all levels are matched by energy and g.
``quanta``
   By the quantum-number text only.
``energy``
   By energy, within ``atol + 1e-5 E`` (``atol`` default 1e-3 cm⁻¹), and equal g.
``index``
   By level number.

Collision partners are identified by their ID. If a file has two blocks with the same
ID (EMAA's H2O partners), the second one is ``"p-H2 #2"``.

Rates
-----

Rates are compared in two ways:

**On the common temperatures, exactly.** A rate counts as changed if it differs by more
than ``rtol`` (default 1e-6, i.e. real changes, not rounding), or if it goes to or from
0. The summary gives the ratio B/A (median and 5–95 %), the temperature whose median is
farthest from 1, and the number of changed rates and transitions with the largest change
(the full list is in :attr:`PartnerDiff.changed_transitions`). All ratios are B/A,
also the largest change ("max change"): the ratio of the rate that differs most from 1
(0.5 and 2 count the same). If B has all rates scaled by 0.72, every ratio is about 0.72.

**On the union grid, interpolated.** The temperatures of both tables in their common
range, with the rates interpolated linearly in log T and log k
(:func:`~pylamda.compare.interpolate_rates`). This shows trends also when the grids
differ and have no temperature in common.

Total rates and critical densities
----------------------------------

Per upper level u and partner, :func:`~pylamda.compare.total_rates` gives the total
de-excitation rate :math:`\sum_l k_{ul}(T)`. ``compare`` compares it for every matched
level on the union grid, from **all** rows of each data set, so added or removed
transitions count too.

:func:`pylamda.critical_density` gives

.. math::

   n_\mathrm{crit,u}(T) = \frac{\sum_l A_{ul}}{\sum_l k_{ul}(T)}.

``compare`` does not compare it: it is not a value in the file, and it changes only if
the rates or the A values change, which are compared themselves.

.. code-block:: python

   data = pylamda.read("co.dat")
   ncrit = pylamda.critical_density(data, "p-H2")          # (n_levels, ntemp) [cm-3]
   ncrit_100 = pylamda.critical_density(data, "p-H2", [100.0])

Summaries per E_up range
------------------------

A change of the high levels can hide behind many unchanged low levels. The summaries on
the union grid, of the total rates and the A values are
therefore also given per upper-level energy range: 0–100, 100–300, 300–1000,
1000–3000, 3000–10000 K and above (:data:`~pylamda.compare.EUP_BINS`). A transition
or level counts as changed in these summaries if its ratio differs from 1 by more than
``limit`` (default 0.1, as in TEST_LAMDA_FILES) at some temperature.

.. code-block:: text

   p-H2: 820 -> 3240 transitions, T 2-3000 K -> 2-3000 K
     ...
     by E_up of the upper level (transitions changed by > 10%): 0-100 K 15 of 15 (max change B/A 2.11), 100-300 K 27 of 30 (max change B/A 1.93), ...
     total de-excitation rate B/A: 23 of 40 levels changed by > 10%, median 0.95, max change 0.71 (level 1 at 2 K); ...

Only the exact comparison decides whether ``bool(diff)`` is True. The interpolated and
per-level summaries add information to a partner that differs.
