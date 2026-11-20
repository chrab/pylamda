Development
===========

Environment
-----------

The development environment holds the dependencies of the code, the tests and the
documentation. You can use pip or conda.

**With pip** (a virtual environment; CI uses the same ``dev`` extra):

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate
   pip install -e ".[dev]"

pip cannot install ``gfortran``. Install it with your system package manager (for
example ``apt install gfortran``) if you want to run the Fortran test.

**With conda** (``environment.yml``, includes ``gfortran``):

.. code-block:: bash

   conda env create -f environment.yml
   conda activate PYLAMDA
   pip install --no-deps -e .

``--no-deps`` makes pip install only the link to the source, so conda keeps all
dependencies. Do not use ``pip install -e ".[dev]"`` in the conda environment: pip then
installs a second copy of packages that conda already has.

Packaging uses Poetry. ``poetry.toml`` sets ``virtualenvs.create = false``, so Poetry
installs into the active environment (venv or conda).

Commands
--------

.. code-block:: bash

   pytest                                          # tests
   pytest --cov                                    # tests with a coverage report
   pytest --cov --cov-report=html                  # coverage as HTML in htmlcov/
   pytest --corpus /path/to/lamda_files            # also the round trip on real files
   ruff format . && ruff check .                   # formatting and linting
   sphinx-build -W -b html docs docs/_build/html   # documentation (warnings are errors)

``--corpus`` takes a directory (all its ``.dat``, ``.radex`` and ``.txt`` files), a
file, or a quoted glob pattern (``--corpus "data/*Lambda*.dat"``), and can be repeated.

``pytest --online`` also runs the tests that use SciX and doi.org; SciX needs
``SCIX_API_TOKEN`` (for example as a conda variable of ``PYLAMDA``).

Coverage uses pytest-cov. The settings are in ``[tool.coverage.*]`` of
``pyproject.toml``: branch coverage of the package ``pylamda``, with the missing lines
in the report. The ``tools/`` scripts are not measured.

Tests that need ``gfortran`` (a RADEX-style Fortran read of a written file) are skipped
if it is not on the ``PATH``.

Layout
------

.. code-block:: text

   src/pylamda/
     model.py        LamdaData, CollisionSet, SourceMap; partner IDs and names
     reader.py       read(): positional reader, collects issues, fills the SourceMap
     writer.py       write(): normalised layout
     verify.py       verify(): the checks L, R, C, F, P (plus the reader's S issues)
     compare.py      compare(), match_levels(), critical_density()
     bib.py          bibtex(): BibTeX export, SciX and doi.org, token, cache
     catalogue.py    the reference catalogue: text parsing, lookup
     refsearch.py    SciX search for free-text references (bib --search, the tool)
     fix.py          fix(), patch(), add_reference(): normalised and patch mode
     combine.py      add_partner(), drop_partner()
     data/           database_refs.bib (LAMDA and EMAA papers); references.json/.bib (catalogue)
   tools/
     catalogue.py    maintainer tool: builds the reference catalogue
     catalogue_review.json   the reviewed corpus
     issues.py       Issue, the check codes (CODES), LamdaReadError
     formula.py      formula mass of a molecule name (check P04)
     references.py   identifiers in reference texts (DOIs, bibcodes, arXiv IDs)
     cli.py          the command line (check, compare, bib, fix, add-partner, drop-partner)
   tests/
     data/           real LAMDA and EMAA files used as fixtures
     helpers.py      make_lamda() for small synthetic files (no issue with the defaults)
   docs/             Sphinx documentation; design.md is the design document

Conventions
-----------

- Type annotations everywhere; NumPy-style docstrings for all public functions, with
  units in square brackets.
- ruff for formatting and linting, line length 120.
- Tests for every new function and every issue code. A new code goes into
  ``issues.CODES`` and ``docs/checks.rst`` (a test checks this).
- The documentation must build with ``sphinx-build -W``, as on Read the Docs.

The reference catalogue (maintainers)
-------------------------------------

``tools/catalogue.py`` (not installed) builds the catalogue of free-text references
(design, section 10.1). It needs the SciX token (``SCIX_API_TOKEN``).

.. code-block:: bash

   python tools/catalogue.py search --database prodimo data/*Lambda*.dat -o review.json
   python tools/catalogue.py show review.json --open          # what still needs a decision
   python tools/catalogue.py show review.json --proposals     # the pre-filled proposals
   python tools/catalogue.py confirm review.json              # accept all proposals

``search`` groups equal texts of the same molecule, collects hints (the file's comments,
the LAMDA molecule page, the EMAA species page; pages cached in
``~/.cache/pylamda/catalogue``), searches SciX by first author and year, and pre-fills a
proposal where a hint names the paper or volume and page match. The reviewer sets
``decision`` per group: ``"ok"`` (the proposal), a bibcode, ``"none"`` or ``"skip"``.
Running ``search`` again keeps the decisions and adds new uses.

Continuous integration
----------------------

GitHub Actions (``.github/workflows/tests.yml``) runs ruff, the tests on Python 3.11 to
3.14 with a coverage report, and the documentation build. Read the Docs builds the
documentation from ``.readthedocs.yaml`` with the versions pinned in ``docs/requirements.txt``.
