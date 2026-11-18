Combining
=========

:func:`pylamda.add_partner` adds the collision rates of one partner from another data
set of the same species, for example the electron rates of a CHIANTI conversion to a
LAMDA file that has none. :func:`pylamda.drop_partner` removes the collision data of one
partner.

.. code-block:: python

   import pylamda

   target = pylamda.read("oatom.dat")
   source = pylamda.read("o_1_chianti.dat")
   data, added = pylamda.add_partner(target, source, "e")
   pylamda.write(data, "oatom_with_e.dat")

   data, dropped = pylamda.drop_partner(data, "He")

.. code-block:: bash

   pylamda add-partner oatom.dat o_1_chianti.dat e -o oatom_with_e.dat
   pylamda drop-partner oatom_with_e.dat He -o oatom_without_he.dat

Adding a partner
----------------

1. The levels of the source are matched to those of the target
   (:func:`pylamda.match_levels`; by default by quanta, else by energy within 0.1 cm⁻¹
   + 1e-4 E and equal weight g). Different sources give the same level with slightly
   different energies (O I ³P₁: LAMDA 158.2687, CHIANTI 158.265 cm⁻¹).
2. The transitions of the partner between matched levels are added, with the level
   numbers of the target. Temperatures and rate coefficients are copied unchanged.
   Transitions to levels the target does not have are left out.
3. A note line records where the partner comes from:
   ``! Collision partner e added from CHIANTI 11.0.2, ion o_1: 3 of 10 transitions,
   levels matched by energy (max |dE| 0.00825 cm-1)``.

``add_partner`` refuses if a matched level has another weight g, or if the target orders
two levels the other way round (a transition would go from the lower to the higher
level). A partner the target has already is only replaced with ``replace_existing``
(``--replace``). The record :class:`~pylamda.combine.AddedPartner` holds the level map,
the numbers of transitions added and left out, and warnings (another molecule name, H2
next to p-H2/o-H2).

Partners are identified by their ID. If a file has several blocks with the same ID
(EMAA's H2O partners, checks C03 and C11), choose the block with ``block=``
(``--block``).

Both functions return a copy and never change their input. The command line writes the
normalised layout to a new file; it refuses to overwrite its input.
