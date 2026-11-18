References and BibTeX
=====================

:func:`pylamda.bibtex` and ``pylamda bib`` write the references of LAMDA files as one
BibTeX file.

.. code-block:: bash

   pylamda bib data/*.dat -o lamda_refs.bib            # offline
   pylamda bib data/*.dat -o lamda_refs.bib --online   # resolved with SciX and doi.org

.. code-block:: python

   import pylamda

   result = pylamda.bibtex(["co.dat", "hco+.dat"], online=True)
   result.write("lamda_refs.bib")
   result.free_text        # references without an identifier (not in the file)

What is collected
-----------------

From each file (:func:`pylamda.collect_references`):

- the dataset DOI (EMAA: ``-- EMAA_DATASET_DOI:`` on line 2);
- the reference of the levels (EMAA: ``Ref:`` on the NLEV line);
- the reference of each collision partner (the text after the partner name);
- identifiers in comment lines, on line 2 and in the notes.

In these texts pylamda looks for DOIs (``https://doi.org/10...``, ``doi:10...``,
``10....``), ADS/SciX bibcodes (``https://scixplorer.org/abs/<bibcode>``, ADS links, or
the bibcode alone) and arXiv IDs (``arXiv:2004.11230``). Every identifier gives one
entry; entries of the same paper (same bibcode or DOI) are merged. A comment above each
entry says what it is cited for (``% CO: collision rates with p-H2``).

References given only as free text ("Flower (2001)") are looked up in the reference
catalogue (below). Those not in it are not guessed: they are listed on standard error
(``pylamda bib``) or in :attr:`BibResult.free_text <pylamda.bib.BibResult.free_text>`.

The reference catalogue
-----------------------

pylamda ships a catalogue of reviewed free-text references
(``pylamda/data/references.json`` and ``references.bib``, :mod:`pylamda.catalogue`). It
was built from the free-text references of LAMDA, EMAA and ProDiMo files: each match of
a text to a SciX record was confirmed by a person. ``pylamda bib`` uses it without
asking and without network or token; such entries are marked "(pylamda reference
catalogue)" in the BibTeX file.

A reference text is split into the parts that cite one paper each ("Flower (2001) &
Wernli et al. (2006)" gives two), and each part is looked up:

1. at the same place: file name, molecule, levels or partner;
2. the same text for the same molecule in another file;
3. the same text anywhere in the catalogue, if it names one paper there. This can be
   wrong for another molecule ("Green 1986"), so ``pylamda bib`` warns about such hits.

``--no-catalogue`` (``catalogue=False``) switches the lookup off. ``pylamda check`` uses
the catalogue for F03 as well. The catalogue is maintained with ``tools/catalogue.py``
(see :doc:`development`).

Offline and online
------------------

Offline, each identifier gives a ``@misc`` entry with the DOI, the SciX link or the
arXiv ID, and a note.

With ``--online``, **SciX is the source of all bibliographic records**:

======================  ==========================================  ==========================================
identifier              with a SciX token                           without a token
======================  ==========================================  ==========================================
bibcode                 SciX export → BibTeX                        ``@misc`` with the SciX link; a warning
arXiv ID                SciX search → bibcode → export              ``@misc``; a warning
DOI                     SciX search → bibcode → export              doi.org (``Accept: application/x-bibtex``)
DOI not in SciX         doi.org                                     doi.org
======================  ==========================================  ==========================================

- All DOIs and arXiv IDs are searched in a few requests
  (``doi:"..." OR doi:"..."``), and all bibcodes are exported in one request.
- The ``adsurl`` of each SciX entry points to ``https://scixplorer.org/abs/<bibcode>``.
- SciX entries use the journal macros of AASTeX (``\apj``, ``\mnras``, ...). With other
  document classes, load them, e.g. from ``aas_macros.sty``.
- Results are cached in ``$XDG_CACHE_HOME/pylamda/bib_cache.json`` (default
  ``~/.cache/pylamda/``); ``--no-cache`` leaves the cache out.

The SciX API token
------------------

SciX needs an API token. Create it in your SciX account settings
(https://scixplorer.org/user/settings/token). pylamda looks for it in this order:

1. ``--scix-token-file FILE`` (``token=`` in Python);
2. the environment variable ``SCIX_API_TOKEN``;
3. the file ``$XDG_CONFIG_HOME/pylamda/scix_token`` (default
   ``~/.config/pylamda/scix_token``). pylamda warns if other users can read it
   (``chmod 600``).

No token is shipped. pylamda never writes the token to files, logs, the cache or
messages.

In a conda environment, the token can be set for that environment only:

.. code-block:: bash

   conda env config vars set -n PYLAMDA SCIX_API_TOKEN=...

Keys
----

Keys are ``<first author><year>`` from the resolved entry (``Yang2010``), or guessed from
the reference text just before the identifier, else made from the identifier. A second
key with the same name gets ``a``, ``b``, ... (``Yang2010a``). ``--keys FILE`` reads and
updates a JSON file identifier → key, so that keys stay the same across runs.

The database papers
-------------------

The main papers of the databases are always in the file, from
``pylamda/data/database_refs.bib`` (shipped with the package; no token or network
needed):

=================================  ==========================================================
when                               entries
=================================  ==========================================================
always                             Schöier et al. 2005 (``Schoier2005``) and van der Tak et al.
                                   2020, the 2020 LAMDA status review (``vanderTak2020``)
a file comes from EMAA             Faure et al. 2025 (``Faure2025``) and the EMAA dataset DOI
                                   10.17178/EMAA (``EMAA``)
=================================  ==========================================================

This follows the citation requests on the LAMDA and EMAA web pages (checked
2026-11-18); the LAMDA status review is added as well. LAMDA also asks to give the date
when the database was accessed; EMAA asks for an acknowledgement sentence. Both requests are written as comments at the top of the
BibTeX file. Both databases also ask to cite the original papers of each file: these
are the entries from the files. ``--no-database-refs`` leaves the database papers out.
