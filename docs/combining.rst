Combining
=========

:func:`pylamda.add_partner` adds the collision rates of one partner from another data
set of the same species, for example the electron rates of a CHIANTI conversion to a
LAMDA file that has none. :func:`pylamda.drop_partner` removes the collision data of one
partner. :func:`pylamda.scale_partner` multiplies the rates of one partner by a factor,
for example He rates that a file gives as H2 rates.

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

Scaling a partner
-----------------

Rates of one partner are often scaled rates of another one. The common case: He rates
multiplied by about 1.37 (the square root of the ratio of the reduced masses) and given
as H2 rates, e.g. ``1 HCN + H2 ! scaled (*1.37) HCN-He from Green & Thaddeus (1974)``.
Checks P02 and C12 find such blocks. A code that scales He rates itself (ProDiMo) then
scales twice. :func:`pylamda.scale_partner` handles three cases:

.. code-block:: python

   import pylamda

   data = pylamda.read("hcn.dat")

   # change the partner: H2 rates that are He x 1.37 become He rates again
   he, scaled = pylamda.scale_partner(data, "H2", 1 / 1.37, to="He", reference="Green & Thaddeus (1974)")

   # scale in place (same ID): correct a factor 1.39 to 1.37
   fixed, scaled = pylamda.scale_partner(data, "H2", 1.37 / 1.39)

   # derive a new partner: keep He, add H2 = He x 1.37
   both, scaled = pylamda.scale_partner(he, "He", 1.37, to="H2", keep=True)

.. code-block:: bash

   pylamda scale-partner hcn.dat H2 1/1.37 --to He --reference "Green & Thaddeus (1974)" -o hcn_he.dat
   pylamda scale-partner hcn_he.dat He 1.37 --to H2 --keep -o hcn_he_h2.dat

1. The partner is found by its ID, as everywhere in pylamda. ``block=`` (``--block``)
   chooses among several blocks with one ID.
2. The rates are multiplied by the factor. The factor is always given explicitly; a
   factor written in the partner text is never used.
3. The scaled block gets the normalised partner line
   ``<id> <molecule> + <partner> -- Ref: <reference>``. The reference text stays as it
   was, unless ``reference=`` (``--reference``) gives a new one.
4. A note line records the operation:
   ``! Collision partner H2 (block 1, 435 transitions) scaled by 0.729927 and relabelled He``.

The scaled rates are rounded to as many significant digits as the rates of the source
have in the file (at least 3; 6 if the data were not read from a file). ``digits=``
(``--digits``) sets another number.

``scale_partner`` refuses a factor that is not finite and > 0, a partner outside the
standard, ``keep=True`` with the same partner (two blocks with one ID, check C03), and a
partner ``to`` that the data have already (unless ``replace_existing``, ``--replace``).
The record :class:`~pylamda.combine.ScaledPartner` holds warnings about:

- a partner text that names another partner than the ID (checks C02, C11, C12): the
  block was taken by its ID;
- a factor in the partner text (``*1.37``, ``x 1.4``, ``by 1.39``) that is neither the
  factor nor its inverse;
- an unusual direction: He to H2 with a factor < 1, or H2 to He with a factor > 1;
- a partner pair other than He and H2 (p-H2, o-H2);
- H2 next to p-H2 or o-H2 in the result (check C04);
- a kept reference text that mentions scaling: it may no longer describe the rates.
