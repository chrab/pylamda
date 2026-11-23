# pylamda: design

Status (2026-11-18): design with decisions taken (section 13). Moved from astrolinedata
(`docs/design_lamda.md`, commit `e0c038f`). Steps 1 to 6 of section 14 are implemented
(package skeleton, data model, reader, writer; `verify` and `pylamda check`; `compare`
and `pylamda compare`; `bibtex` and `pylamda bib`; `fix`, `patch` and `pylamda fix`;
the reference catalogue with `bib --search` and `fix --add-refs`) and step 7
(`add_partner`, `drop_partner`; astrolinedata and TEST_LAMDA_FILES use pylamda,
section 11). `scale_partner` was added on 2026-11-23 (section 9).
Step 8 and check F04 are left for later (section 15).

Section 1 shows what exists, section 2 the scope, sections 3–10 the design, section 11
the relation to astrolinedata, section 13 the decisions.

**Guiding rule: pylamda follows the LAMDA standard strictly**
(<https://home.strw.leidenuniv.nl/~moldata/molformat.html>):

- the partner ID decides;
- only the seven standard partners exist.

Deviations are reported, not interpreted.

## 1. Does something like this exist?

Checked on 2026-11-18 (web search, PyPI). A GitHub search was not possible (API rate
limit).

| tool | what it does with LAMDA files | what it does not do |
|---|---|---|
| [astroquery.lamda](https://astroquery.readthedocs.io/en/latest/lamda/lamda.html) | downloads files from the LAMDA website; `parse_lamda_datafile` reads a file into astropy tables; `write_lamda_datafile` writes them back | no checks, no fixes, no comparison, no references; needs astropy |
| [pythonradex](https://pythonradex.readthedocs.io/en/latest/notebooks/reading_LAMDA_files.html) | `LAMDA_file.read` reads levels, lines and rates for its own RADEX solver; interpolates rates in T | reader only, part of a solver |
| [ndradex](https://github.com/astropenguin/ndradex), [SpectralRadex](https://libraries.io/pypi/spectralradex) | run RADEX on LAMDA files | no file handling of their own |
| RADEX, ProDiMo, LIME, … | each has its own Fortran/C reader | each reads the format a bit differently (section 3.3) |

- No package had the name `pylamda` (PyPI: 404).
- `lamda` on PyPI is an unrelated Android tool.

**Conclusion:**

- Reading LAMDA files is solved several times.
- No tool **checks**, **fixes**, **compares** or **combines** LAMDA files, and none
  **collects their references**.
- The format is a community standard: LAMDA, EMAA, ProDiMo's own files, CHIANTI
  conversions. That gap therefore matters to more people than us.

## 2. Scope

**In scope:**

1. A well-defined Python data structure for the content of a file.
2. Read and write (lossless round trip of the data).
3. Verify: structure, counts, IDs, level order, references, physics checks.
4. Compare two files in a meaningful way (levels, lines, rates).
5. Combine files (`add_partner` and similar).
6. A BibTeX file of all references of one or more files.
7. Fix files: partner ID/name, level order, missing columns, references, …

**Out of scope:**

- Downloads from LAMDA/EMAA/CHIANTI: these stay in astrolinedata.
- Radiative transfer (RADEX, pythonradex).
- ProDiMo slab runs (`TEST_LAMDA_FILES/slab.py`): these stay in a ProDiMo tool and use
  pylamda.
- Conversions from other databases (CHIANTI → LAMDA): these stay in astrolinedata.

**Dependencies and packaging:**

- Python ≥ 3.11. Sphinx ≥ 8.2 needs 3.11, and code and documentation use one
  environment.
- Required: numpy and requests (requests resolves DOIs and bibcodes; the BibTeX export
  is an integral part of the package). No astropy, so that Fortran-centred groups can
  use the package without a large stack.
- Optional extras:
  - `plot`: matplotlib;
  - `dev`: pytest, pytest-cov, ruff, Sphinx with the extensions of the documentation, and
    the `plot` extra.
- Install for development: `pip install -e ".[dev]"`, or `poetry install --all-extras`.
- Packaging with **Poetry**:
  - `pyproject.toml` with a `[project]` table (PEP 621) and the build backend
    `poetry-core`;
  - `poetry.toml` with `virtualenvs.create = false`, so Poetry installs into the active
    conda environment;
  - the same setup as astrolinedata.
- Development environment: one conda environment, `PYLAMDA`.
- Documentation: Sphinx, on **Read the Docs**:
  - `.readthedocs.yaml`;
  - `docs/requirements.txt` with pinned versions;
  - `fail_on_warning: true`.

## 3. What we have already

### 3.1 Code

| where | content | state |
|---|---|---|
| astrolinedata `formats/lamda.py` (~1300 lines) | `LamdaData`, `CollisionSet`, `read_lamda` (by position, counts decide, notes after the data), `write_lamda` (lossless, `-- Ref:` style), `partner_name` (aliases), `find_dois`, `compare_lamda` / `LamdaDiff` / `PartnerDiff`, `add_partner` / `AddedPartner` | tested (203 tests in the package), used by 3 backends |
| astrolinedata `backends/nlte.py` | `nlte_metadata`, `data_references` (DOIs → reference dicts) | tested |
| astrolinedata `references.py` | `export_bibtex` (merges reference dicts by key); entries `Schoier2005`, `vanderTak2020`, `Faure2025`, `EMAA` | tested |
| `/home/rab/MODELS/TEST_LAMDA_FILES/lamda.py` | second reader that follows ProDiMo (`map_collider`: name first, then ID); `match_levels`; `compare_rates` (log-log interpolation onto the union T grid); `compare_spectroscopy`; `prodimo_corrections` | tested (pytest), standalone script |
| `TEST_LAMDA_FILES/report.py`, `slab.py` | text summary, PDF plots (ratio vs E_up), ProDiMo slab models | ProDiMo-specific |

There are two readers, two level matchers and two rate comparisons. pylamda merges each
pair into one.

### 3.2 Findings from the file work (2026-11)

These problems were found in real files. Each one is a check or fix in pylamda.

| finding | file(s) | becomes |
|---|---|---|
| Partner name and partner ID disagree | several ProDiMo files, `hcn@hfs.dat`, `no+.dat` | check C02, fix |
| EMAA's extra partners p-H2O, o-H2O, CO (not in the standard) reuse IDs 2, 3, 5 | some EMAA files (e.g. CO) | check C11; not supported (3.3) |
| More collision blocks than NPART | `oatom@lique.dat` (6th block O+e) | check S03 |
| More temperatures/values than NTEMP | `s-h.dat`, `si-h.dat` | check S04 |
| Comment line without `!` | ProDiMo CH+ file | check S05, fix |
| Partner line ending in ` !` | ProDiMo files | reader robustness |
| Excitation rates in the rows of de-excitation rates | `Oatom_Lambda_extra_NIST*` (O+H) | check P01 (heuristic) |
| Rates "scaled from H2" but labelled He; ProDiMo scales again | O2, SO2 | checks P02 (text hint), C12; `scale_partner` (section 9) |
| CHIANTI level order differs from energy order (inverted pairs) | CHIANTI conversions | check C08, fix (explicit) |
| References missing or only as free text | many ProDiMo files | checks F01–F03, `bib` report |

### 3.3 Collision partners: the ID decides

The standard defines seven partners by ID: 1 = H2, 2 = p-H2, 3 = o-H2, 4 = e, 5 = H,
6 = He, 7 = H+. The rest of the partner line is free text.

pylamda:

- identifies each partner **by its ID only**, as RADEX does;
- reads the name in the free text (`<mol> + <name>`) **only for check C02**.

Not every reader does it this way:

- ProDiMo (`MAP_COLLIDNAME`) uses the name first and the ID only as a fallback.
- If name and ID disagree, ProDiMo and RADEX therefore use different rates.
- Check C02 reports this. It would have found most of the ProDiMo partner problems at
  once.

**EMAA**:

- For the standard partners, EMAA's IDs are correct.
- Some EMAA files also have p-H2O, o-H2O or CO as partners. These are not in the
  standard and reuse IDs 2, 3 and 5.
- pylamda does not support these partners. Check C11 reports them.
- Read by ID, such a block would be taken as p-H2, o-H2 or H. The check makes this
  visible.

astrolinedata's `read_lamda(partners="auto"|"name")` is not carried over.

## 4. Data model

Based on astrolinedata's `LamdaData`: one numpy array per column, 1-based level numbers
as in the file. That design has worked well.

```text
LamdaData
├── molecule, mass, title, dataset_doi, source (the !SOURCE: text)
├── levels:     level_energy [cm-1], level_g, level_quanta (str), level_type, level_reference
├── lines:      line_up, line_low, line_a [s-1], line_freq [GHz] (NaN if absent),
│               line_eup [K] (NaN if absent)
├── collisions: list[CollisionSet]
│     partner_id, partner (standard name of the ID, a property), partner_text (line as written),
│     text_partner (the name in partner_text, a property),
│     temperatures [K], up, low, rates [cm3 s-1] (downward), reference
├── notes:      lines after the last declared block
└── source_map: SourceMap | None (line numbers of every item, raw lines, the reader's issues;
                only if read from a file)
```

The SourceMap attribute is `source_map`, not `source`: astrolinedata's `source` (the
`!SOURCE:` text) keeps its name. `data.issues` gives the reader's issues.

Changes from astrolinedata's `LamdaData`:

- **`SourceMap`**: the reader keeps the line number of each level, line, partner line
  and rate row, and the raw lines. Checks can then say "line 312" and not only "rate
  row 17 of p-H2". A fix can edit single lines (section 7).
- **Partner by ID only**:
  - `partner` is always the standard name of `partner_id` (3.3);
  - the raw partner text is kept, so that checks C02 and C11 can read the name in it.
- **Issues, not log messages**: the reader collects every deviation (surplus values,
  extra blocks, …) as an `Issue` (section 6). `verify` reuses the reader's issues.
  astrolinedata's reader only logs warnings.
- **Name**: the class stays `LamdaData` (decision 7), so that astrolinedata code does
  not change.

Methods stay small: `partner(name)`, `select_partners()`, `references()`, `rate(T)`,
`copy()`. All operations (`verify`, `compare`, `fix`, `add_partner`, `bibtex`) are
functions, not methods.

## 5. Read and write

**Reading** follows astrolinedata's `read_lamda` (rules in astrolinedata
`docs/design_nlte.md`, section 1.1):

- By position. The counts decide. Header text is not used.
- Lines whose first non-blank character is `!` are comments. Blank lines are skipped.
- Surplus values on a line are ignored and give an issue.
- Everything after the last declared block is kept as `notes`.
- `strict=True` turns every structure issue into an error, for CI.
- Not a LAMDA file (added 2026-11-23): the reader takes lines only as needed. If the
  header does not read (molecule name not a number, molecular weight a number > 0,
  NLEV >= 1, the level rows, NLIN), it raises `NotLamdaError` (a subclass of
  `LamdaReadError`) with a one-line `reason` and does not read the rest of the file.
  Binary data (NUL) in the header and lines longer than `MAX_LINE_LENGTH` (10^6
  characters, never read whole) count as well. `pylamda check` prints one line per such
  file (`<path>: not a LAMDA file: <reason>`), not its structure issues. On ProDiMo's
  `data/` (177 files): 70 files of other formats are recognised; before, they gave about 9200
  issues.

**Writing** has two modes:

| mode | use | how |
|---|---|---|
| `write(data, path)` | new or converted files | one normalised layout, as astrolinedata's `write_lamda`: `-- Ref:` style, shortest lossless numbers, A always in exponential notation |
| `patch(data, path)` (via `fix`) | existing, hand-maintained files (ProDiMo `data/`) | changes only the lines that changed; comments and layout stay; small git diffs |

The second mode is new. It is what we did by hand in the ProDiMo files: fix an ID, keep
everything else.

## 6. Verify

```python
issues = pylamda.verify("co.dat")  # list[Issue]
issues = pylamda.verify("co.dat", physics=False)  # without the P checks
```

```text
Issue(code="C02", severity="warning", line=312, message="partner line '1 CO + p-H2 ...': ID 1 is H2, the name says p-H2 (read as H2; ProDiMo would read p-H2)", fixable=True)
```

- Every check has a code, a default severity (`error`, `warning`, `info`) and a short
  description, as in ruff.
- Checks can be selected or ignored (`--select`, `--ignore`).
- Output: text, or JSON for machines.
- Exit code: 0 = no errors, 1 = errors, 2 = file not readable.

### 6.1 Checks (first set)

**S: structure** (mostly found by the reader)

| code | check | sev. |
|---|---|---|
| S01 | file shorter than its counts require | error |
| S02 | level/transition numbers not 1..N in order | error |
| S03 | data after the last declared block (more blocks than NPART) | warning |
| S04 | more values on a line than declared (NTEMP) | warning |
| S05 | text line without `!` between blocks (read as data by Fortran) | error |
| S06 | number that Fortran cannot read (e.g. `1.0D-10` is fine, `1,0e-10` is not); also a partner line without the ID | error |
| S07 | level or radiative-transition row with fewer values than required (i, E, g / i, up, low, A) | error |

**L: levels**

| code | check | sev. |
|---|---|---|
| L01 | energies not in increasing order | warning |
| L02 | g ≤ 0, or g not an integer | error / info |
| L03 | two levels with the same E, g and quanta | warning |
| L04 | quanta empty or not unique | info |
| L05 | lowest level E ≠ 0 | info |

**R: radiative transitions**

| code | check | sev. |
|---|---|---|
| R01 | up or low outside 1..NLEV, or up = low | error |
| R02 | E(up) ≤ E(low) | error |
| R03 | A ≤ 0 or not finite | error |
| R04 | frequency ≠ (E_up − E_low) · c (relative tolerance, default 1e-4) | warning |
| R05 | E_up [K] column ≠ E(up) · hc/k | warning |
| R06 | the same (up, low) twice | error |

**C: collisions**

| code | check | sev. |
|---|---|---|
| C01 | partner ID not in 1..7 | error |
| C02 | the name in the partner text is a different standard partner than the ID (3.3) | warning |
| C03 | the same partner twice | error |
| C04 | H2 together with p-H2 or o-H2 (counted twice by some codes) | warning |
| C05 | temperatures not strictly increasing, or ≤ 0 | error |
| C06 | rate < 0 or not finite | error |
| C07 | a rate row that is 0 at all temperatures | info |
| C08 | E(up) < E(low) in a rate row ("inverted": wrong order or an excitation rate) | warning |
| C09 | the same (up, low) twice in one block | error |
| C10 | up or low outside 1..NLEV, or up = low | error |
| C11 | the partner text names a partner outside the standard (EMAA: p-H2O, o-H2O, CO); not supported, the block is read by its ID | warning |
| C12 | the partner text names He, the ID is H2/p-H2/o-H2, and the text says the rates are scaled (`1 SO2 + He ! ... scaled *1.39`): He rates scaled to H2 with the name not changed; reported instead of C02 | warning |

**P: physics heuristics**

- On by default, as `info`.
- They can give false alarms.
- Switched off with `physics=False` / `--no-physics` (the same as `--ignore P`).

| code | check |
|---|---|
| P01 | a row that behaves like an **excitation** rate: at T ≪ ΔE/k, log k falls with slope ≈ −ΔE/k in 1/T. This would have found the O+H rows. |
| P02 | partner text says "scaled" / "scaling" / "from H2": warn that codes that scale He rates (ProDiMo `ncorr`) may scale twice (O2, SO2). |
| P03 | rates outside a plausible range (e.g. > 1e-8 cm3 s-1 for neutral partners) |
| P04 | mass differs from the formula mass of `molecule` (only if the name parses; atomic weights table from astrolinedata) |

**F: references**

| code | check | sev. |
|---|---|---|
| F01 | no reference for the levels/lines | warning |
| F02 | a partner without a reference | warning |
| F03 | a reference without a DOI, bibcode or link (free text only), and not in the catalogue (10.1) | info |
| F04 | a DOI or bibcode that does not resolve (only with `--online`) | warning |

More checks can come later, without API changes.

Implementation notes (step 2):

- A check of many rows (rates, transitions, levels) gives one issue per block or
  section, with the number of rows (`Issue.count`), the first example and its line.
- R04/R05: the tolerance is `rtol` plus the rounding of the numbers as written (from the
  SourceMap); an energy of exactly 0 is exact.
- P01: the slope is fitted over the temperatures below ΔE/2k; a row is reported if the
  slope is between 2/3 and 3/2 of −ΔE/k. Rows labelled with the wrong levels have
  another ΔE and may be missed (the O+H rows in `Oatom_Lambda_extra_NIST.dat`, whose
  grid starts at 50 K).
- P03 limits: 1e-8 cm3 s-1 (neutral partners), 1e-7 (H+), 1e-4 (electrons).
- F04 (`--online`) is left for later (section 15); the catalogue part of F03 is done
  (step 6).

## 7. Fix

```python
fixed, log = pylamda.fix(data, trust="id", sort_levels=True)
```

```bash
pylamda fix in.dat -o out.dat --trust id
pylamda fix in.dat --in-place --trust name   # patch mode: only changed lines
```

**Rule: fix only what is unambiguous. Every ambiguous fix needs an explicit option.**

| issue | fix | needs |
|---|---|---|
| C02 name ≠ ID | `--trust id`: correct the name; `--trust name`: correct the ID. The partner line becomes `<id> <mol> + <name> -- Ref: <ref>` | `--trust id\|name`: no default. The file is read by ID, but we found wrong IDs as well as wrong names |
| L01 levels not ordered | sort by energy, renumber all transitions (pure relabelling) | default on |
| C08 inverted rows | `swap`: exchange up/low and keep the rate (wrong order only); `chianti`: swap and multiply by g(new low)/g(new up) (Υ is symmetric, the rate was divided by the wrong g); `balance`: convert an excitation rate by detailed balance | `--inverted swap\|chianti\|balance`; no default |
| R04/R05 missing columns | compute frequency and E_up from the energies | default on for missing values; differing values need `--recompute` |
| S02 wrong running numbers of rows | correct the first column (no reader uses it); added 2026-11-18 | default on |
| S03 extra blocks | count them in, or turn them into comments | `--extra-blocks include\|comment` |
| S04 surplus values | drop them (what all readers do) | default on |
| S05 comment without `!` | add `!` | default on |
| partner text | normalise to the `-- Ref:` style | default on in normalised mode, off in patch mode |
| F01–F03 references | append a bibcode (as SciX link) or DOI to the reference text (10.2) | `--add-refs` (catalogue or confirmed review file) or `--ref TARGET=ID` |

- **Record**: each fix adds a note line, `! pylamda <version> fix: C02 partner 2: ID 1 -> 2 (trusted name)`,
  as astrolinedata's `add_partner` does. `log` holds the same as data.
- **Never overwrite** the input without `--in-place`.

Implementation notes (step 5):

- Two stages: S05, S03 (then the file is read again), S04, S02 on the raw lines; then the
  data fixes C02/C12, L01, C08, R04/R05 and references.
- Patch mode edits tokens in place (`replace_tokens`): right-aligned columns stay
  aligned, numbers keep their style (`format_like`; rates recomputed by `chianti` or
  `balance` keep the file's digits, at least 3 significant), line endings (CRLF) stay.
  Recomputed frequencies and E_up get their full precision.
- `patch()` and `fix()` give the same data (tested on the corpus).
- C12 is fixed by `--trust` like C02.
- CLI: `-o` = normalised, `--in-place` = patch, `--dry-run` = diff of the patch;
  `--patch`/`--normalise` switch. After fixing, the fixable issues that remain are
  named with the option they need.
- `--add-refs` (from the catalogue or a review file) follows with step 6.

## 8. Compare

Merge `compare_lamda` (astrolinedata) and `compare_rates` (TEST_LAMDA_FILES):

- **Level matching**: by quanta, by energy and g, or by index (`auto` as in
  astrolinedata). The match is shared by `compare`, `add_partner` and `fix`.
- **Lines and rates**: matched through the level map, so that inserted levels do not
  show up as changes of all later entries.
- **Rates on two grids**:
  - the common temperatures exactly (astrolinedata: detects real changes vs rounding);
  - the union grid in the common range with log-log interpolation (TEST_LAMDA_FILES:
    shows trends even when the grids differ).
- **Physically meaningful numbers**:
  - **critical densities** n_crit(T) = A_ul / Σ_l k_ul(T) per upper level and partner:
    one number per level that says whether a rate change matters for NLTE;
  - total de-excitation rate per level;
  - changes summarised per E_up range, so a change at 5000 K levels does not hide
    behind 1000 unchanged low levels.
- **Output**: `bool(diff)`, a text summary, JSON, and plots with the `plot` extra
  (ratio vs E_up per temperature, as `report.py`).

Implementation notes (step 3):

- `match_levels` returns a `LevelMatch` (1-based map, method). `auto`: by quanta
  (normalised as in TEST_LAMDA_FILES: lower case, blanks/`_`/`,`/`;` as one separator)
  if unique in both and at least half match; the rest by energy and g
  (`quanta+energy`); else by energy and g. Default tolerances as astrolinedata's
  `compare_lamda` (1e-3 cm-1 + 1e-5 E).
- `LamdaDiff`/`PartnerDiff` keep all fields of astrolinedata, so its code can switch.
  New: `a_ratio`, and per partner `interp_temperatures`, `rates_interpolated`,
  `critical_density`, `total_rate` (each a `RatioSummary` with per-T and per-E_up
  statistics).
- Only the exact comparison decides `bool(diff)`. The union-grid and per-level summaries
  count a change above `limit` (default 0.1, TEST_LAMDA_FILES' `LIMIT`).
- Partners are compared by ID; a repeated ID is `"<name> #2"`.
- CLI: `pylamda compare A B`, exit code 0/1/2 as `diff`.
- Plots follow with step 8.

## 9. Combine

- `add_partner(target, source, partner)`: moves unchanged from astrolinedata (its
  `docs/design_nlte.md`, decision 7).
- `select_partners(data, names)`, `drop_partner(data, name)`.
- `cut_levels(data, e_max | n_max)`: drop levels and all transitions to them.
- Later, if needed: replace the rates of a partner for some transitions only, merge two
  temperature grids. Not in the first version: these need physics decisions per case.

All combine functions return the new data and a record (as `AddedPartner`), and add a
note line.

Implemented (2026-11-18, `pylamda/combine.py`): `add_partner` from astrolinedata, with
partners by ID; a source or target with several blocks of the same ID (EMAA) needs
`block=`; warnings (other molecule, H2 next to p-H2/o-H2) go into
`AddedPartner.warnings`; the new partner line is normalised with the target's molecule;
level matching by `match_levels` (tolerances 0.1 cm-1 + 1e-4 E). `drop_partner` with a
`DroppedPartner` record. CLI `pylamda add-partner` and `drop-partner` (normalised
layout, never the input file). `select_partners` stays a method of `LamdaData`;
`cut_levels` is left for later (section 15).

Added on 2026-11-23: `scale_partner(data, partner, factor, to=, keep=)` with a
`ScaledPartner` record. It multiplies the rates of one partner by an explicit factor
(never taken from the partner text): in place, as another partner (`to`, e.g. H2 rates
that are He x 1.37 back to He), or as an added partner (`keep`, e.g. H2 = He x 1.37).
Errors: factor not finite and > 0, a partner outside the standard, `keep` with the same
partner (C03), `to` present already without `replace_existing`. Warnings: the partner
text names another partner than the ID (C02/C11/C12), a factor in the text that is
neither the factor nor its inverse, an unusual direction (He -> H2 with < 1, H2 -> He
with > 1), a pair other than He and H2/p-H2/o-H2, H2 next to p-H2/o-H2 (C04), a kept
reference text that mentions scaling. The rates are rounded to the significant digits
of the source as written (at least 3; 6 without a source map). Normalised layout only.
CLI `pylamda scale-partner FILE PARTNER FACTOR` (FACTOR also as `a/b`); it prints the C
issues of the result.

## 10. References and BibTeX

```bash
pylamda bib data/*.dat -o lamda_refs.bib            # offline: catalogue + @misc entries
export SCIX_API_TOKEN=...                           # the user's own SciX token
pylamda bib data/*.dat -o lamda_refs.bib --online
```

1. **Collect** from: dataset DOI, level reference, partner references, comment lines and
   notes.
2. **Extract** identifiers:
   - DOIs (`find_dois`);
   - bibcodes from SciX and ADS links (`scixplorer.org/abs/<bibcode>`,
     `ui.adsabs.harvard.edu/abs/<bibcode>`, as now in the ProDiMo files);
   - arXiv IDs.
3. **Resolve** (`--online` only). **SciX is the source for all bibliographic records.**

   | identifier | with a SciX token | without a token |
   |---|---|---|
   | bibcode | SciX export → BibTeX | not resolved; `@misc` with the SciX link; warning that a token is needed |
   | arXiv ID | SciX search (`identifier:`) → bibcode → export | not resolved, as above |
   | DOI | SciX search (`doi:`) → bibcode → export | doi.org (`Accept: application/x-bibtex`) |
   | DOI not in SciX (datasets: EMAA, Zenodo, …) | doi.org | doi.org |

   - All SciX requests in one batch (the export takes a list of bibcodes).
   - Each entry gets an `adsurl` pointing to `https://scixplorer.org/abs/<bibcode>`.
   - Results are cached (bibcode → BibTeX, DOI → bibcode).
4. **The SciX API token**:
   - The user must provide the token (SciX account settings).
   - pylamda looks for it in this order:
     1. the `token=` argument (`--scix-token-file` on the command line);
     2. the environment variable `SCIX_API_TOKEN`;
     3. the config file `~/.config/pylamda/scix_token` (`$XDG_CONFIG_HOME` respected).
        pylamda warns if other users can read the file.
   - No token is shipped, and there is no default.
   - The token is never written to files, logs, the cache or error messages.
   - SciX uses the same API as ADS. Confirmed on 2026-11-18 with a live token: base URL
     `https://api.scixplorer.org/v1` (also `https://api.adsabs.harvard.edu/v1`, same
     token); `GET /search/query` (`q`, `fl`, `rows`), header `Authorization: Bearer
     <token>`; `POST /export/bibtex` with `{"bibcode": [...]}` -> `{"export": "..."}`. In a
     search, `doi:"a" OR doi:"b"` works; `doi:("a" OR "b")` gives HTTP 400.
5. **Without resolution**: an `@misc` entry with `doi` (or the SciX link) and a `note`
   saying what the reference is for ("CO: collision rates with p-H2"), as astrolinedata's
   `data_references` does.
6. **Free text without an identifier** ("Flower (2001)"):
   - first looked up in the shipped reference catalogue (10.1);
   - otherwise not guessed: listed in a report (F03), with SciX candidates if
     `--search` is given (10.1);
   - identifiers can then be added to the file (10.2).
7. **Keys**: from the resolved entry (`Yang2010`), made unique (`Yang2010a`); a
   `--keys file` option keeps stable keys across runs.
8. **The main papers of the databases are always included.** This may give duplicates,
   but no database citation is missed.

   | when | papers |
   |---|---|
   | always (LAMDA format and database) | Schöier et al. 2005 (A&A 432, 369; DOI 10.1051/0004-6361:20041729); van der Tak et al. 2020 (Atoms 8, 15; DOI 10.3390/atoms8020015) |
   | any EMAA hint in a file (dataset DOI `10.17178/EMAA_…`, "EMAA" in a reference or the source) | Faure et al. 2025 (A&A 700; DOI 10.1051/0004-6361/202554403); the EMAA dataset DOI 10.17178/EMAA (EMAA's citation policy) |

   - The citation requests on the LAMDA and EMAA web pages are checked again during
     implementation. The list follows them.
   - Checked on 2026-11-18: the LAMDA home page asks to cite Schöier et al. 2005 only;
     it links van der Tak et al. 2020 (arXiv:2004.11230) as a review, not as a required
     citation. Decided on 2026-11-18 to keep both in every export (decision 10). EMAA
     (citation policy page) asks for Faure et al. 2025, the dataset DOI 10.17178/EMAA and
     an acknowledgement sentence. LAMDA also asks to give the access date. Both requests
     go into the BibTeX file as comments.
   - `--no-database-refs` leaves these papers out.
   - If an entry also comes from a file (same DOI), it appears once.
9. **Where the main papers come from**: they are **shipped with the package** as a
   plain BibTeX file (`pylamda/data/database_refs.bib`), with a short README entry.
   - Offline: no SciX token or network needed for the most important citations.
   - Easy to find and read, for people and for other tools and AI assistants that use
     pylamda.
   - astrolinedata already has these four entries (`references.py`: `Schoier2005`,
     `vanderTak2020`, `Faure2025`, `EMAA`); they move to pylamda.
   - The other option, to fetch them from SciX at run time, is not needed for four
     stable entries. SciX is still used to check them in the `--online` tests.

### 10.1 Free-text references: a curated catalogue

Many existing files give a reference only as free text ("Yang et al. (2010)",
"Flower 2001"). The molecule pages of LAMDA and EMAA often have more: full citations or
ADS links. Searching for these on every run would be slow, needs a token, and gives
different results over time. The search is therefore done **once, during development**,
and the results ship with the package.

**The catalogue** (`pylamda/data/references.json` + `references.bib`):

- One entry per reference: bibcode, DOI, BibTeX key, and the reviewed status.
- One **lookup row** per place where a free-text reference occurs, keyed by:
  - database, molecule, file name and partner;
  - the normalised reference text (lower case, whitespace and punctuation collapsed).
- The context is needed because the same text ("Yang et al. 2010") means different
  papers in different files.
- A lookup by text alone is a fallback. It is used only if the text is unique in the
  catalogue.
- Offline, no token needed. It is also readable by other tools and AI assistants
  (plain JSON and BibTeX).

**Building the catalogue** (a maintainer script in `tools/`, not part of the installed
package):

1. Collect all free-text references from the corpus: all LAMDA and EMAA files, and
   ProDiMo's files.
2. For each one, collect hints:
   - the LAMDA molecule page (astrolinedata's `parse_page` already reads these pages);
   - the EMAA species page;
   - comments in the file.
3. Search SciX with author, year, journal and title words from the text and the
   hints.
4. Write candidates with a score into a review file.
5. **A person confirms each match.** Only confirmed entries go into the catalogue.
   Unclear cases stay out and remain F03 issues.
6. Rerun when LAMDA or EMAA add files (astrolinedata sees new versions in LAMDA's
   `database.history`).

**State (2026-11-18):** `tools/catalogue.py` does steps 1-5 above and `build`:

- `search`: grouping (one decision per molecule and normalised text); hints from the
  file's comments, the LAMDA molecule page (reference links to ADS) and, for EMAA
  files, the EMAA species page; SciX queries by first author and year, by species, by
  volume/page/year, and by compound surname; scores (volume, page, journal, first
  author, co-authors, species and partner in the title with charge states and parent
  molecule, journal article, named by a hint); a pre-filled `proposal` if a hint names
  the paper or volume and page match, a weaker "best candidate" proposal otherwise.
- `confirm` (`--include-best`, `--include-first-pass`), `show` filters, `build`.
- The review file `tools/catalogue_review.json` is kept in the repository.
- First catalogue built from the existing corpus (75 ProDiMo files, the LAMDA and EMAA
  test files, `cn.dat`, one CHIANTI file): 126 groups, 167 uses; 114 confirmed (76
  tool proposals, 38 from a reviewed first pass by Claude), 9 skipped, 3 open (not in
  SciX: Faure 2014 C2H+e, S. Green HNCO, Faure 2012 OH+e). Result: 72 entries, 154
  lookup rows in `pylamda/data/references.json` and `.bib`.
- Run time: `pylamda/catalogue.py` (text splitting, normalisation and parsing moved
  from the tool, which imports them); lookup by place, then molecule and text, then
  text alone if unique (with a warning in `bib`); `pylamda bib` takes found papers
  from the shipped `references.bib` (source `catalogue`, offline); F03 does not report
  references that the catalogue knows. On the corpus: 70 catalogue entries in `bib`,
  F03 194 -> 61.
- `pylamda/refsearch.py`: the search part of the tool (hints, scoring, proposals,
  review file) moved into the package; `pylamda bib --search` writes candidates for
  references not in the catalogue to a review file (same format as the tool), never to
  the BibTeX file; `python -m pylamda.refsearch FILE` shows it.
- `pylamda fix --add-refs [REVIEW]`: SciX links of free-text references from the
  catalogue (only matches by place or molecule) or from the decisions of a review file.
  Links are written with `%26` for `&`; `find_bibcodes` also reads a raw `&`. Fix note
  lines are not reference texts for `check` and `bib`.
- `bibtex()` also takes the reviewed BibTeX of identifiers (DOI, bibcode) that the
  catalogue knows, offline.
- Step 6 is complete.

**At run time**:

- `pylamda bib` and `pylamda check` use the catalogue without asking.
- A reference resolved through the catalogue is marked as such in the report.
- `--search` (needs a SciX token) searches SciX for references that are not in the
  catalogue:
  - results are only **candidates**, shown in the report and written to a review file;
  - they are never added to the BibTeX output or a file without a confirmation (10.2).
- Users can send confirmed matches back as contributions to the catalogue.

### 10.2 Adding identifiers to files (reference fix)

A file can get its identifiers written in, so that later readers need no catalogue:

```bash
pylamda fix co.dat --in-place --add-refs                 # from the catalogue
pylamda fix co.dat --in-place --add-refs review.json     # confirmed candidates
pylamda fix co.dat --in-place --ref p-H2=2010ApJ...718.1062Y --ref levels=10.1051/0004-6361:20041729
```

```python
data, log = pylamda.add_reference(data, "p-H2", bibcode="2010ApJ...718.1062Y")
```

- Targets: `levels` (level/line reference), a partner name, or `dataset`.
- The identifier is appended to the existing text, which stays unchanged:
  - bibcodes as a SciX link, `https://scixplorer.org/abs/<bibcode>`, as now in the
    ProDiMo files;
  - DOIs as `doi:10.xxxx/...`.
- Example: `2 CO + p-H2 -- Ref: Yang et al. (2010) https://scixplorer.org/abs/2010ApJ...718.1062Y`.
- Only free text changes: on the partner line after the partner name, and on the NLEV
  line after `Ref:`. All readers ignore that text, so the data are read the same.
- Patch mode changes only these lines. The fix adds a note line, as all fixes do.
- An identifier already in the text is not added again.

## 11. Relation to astrolinedata and TEST_LAMDA_FILES

- astrolinedata (`/home/rab/git/astrolinedata`) depends on pylamda:
  - `formats/lamda.py` becomes a thin re-export for one release, then goes;
  - backends, `.json`/`.kv` metadata and `api.add_partner` (files plus metadata) stay
    in astrolinedata.
- `TEST_LAMDA_FILES/lamda.py` is replaced by pylamda:
  - `slab.py` and `report.py` stay ProDiMo-specific;
  - later perhaps they move to prodimopy as `compare_lamda_slab`.
- ProDiMo `data/`: a CI job `pylamda check data/*.dat` catches new ID, order or reference
  problems.
- ProDiMo reads partners by name first (3.3), so the ProDiMo files must have C02 at zero.
  The CI job can enforce this with `--select C02 --error-on C02` (raise C02 to an error).
- Neither pylamda nor astrolinedata supports partners outside the standard (water, CO,
  ...). pylamda reports them (C11, also C01 and C03); astrolinedata refuses such a file
  as not supported (decided 2026-11-18; before, it dropped those blocks).
- Git history: pylamda starts without the history of `formats/lamda.py` (decision 8).
  The commit that adds the code names the astrolinedata commit it comes from.

Done (step 7, 2026-11-18):

- astrolinedata (branch `use-pylamda`): `formats/lamda.py` and its tests are removed,
  with no re-export module (astrolinedata was not released yet); all modules import
  pylamda (dependency `pylamda @ git+https://github.com/chrab/pylamda.git`). The
  backends read by ID; `backends/nlte.check_supported` refuses a file with partners
  outside the format (checks C01, C03, C11) as not supported; `log_issues` logs pylamda's
  structure issues; `api.add_partner` calls `pylamda.add_partner` and logs its warnings.
  `ald.read_lamda`/`write_lamda`/`compare_lamda` are gone; users call pylamda.
- TEST_LAMDA_FILES: `lamda.py` keeps its name and interface but reads with
  `pylamda.read`, matches levels with `pylamda.match_levels` and interpolates with
  `pylamda.compare.interpolate_rates`. ProDiMo-specific parts stay there: the collider
  names as ProDiMo's MAP_COLLIDNAME reads them (name first), the second block of a name
  wins, `prodimo_corrections`, and the containers for the report. pylamda's structure
  issues appear in the report's messages.

## 12. Tests

- pytest. Tests run on Python 3.11, 3.12, 3.13 and 3.14 (GitHub Actions); locally in
  the conda environment `PYLAMDA`.
- Coverage with pytest-cov, opt-in (`pytest --cov`); branch coverage of `pylamda`,
  settings in `[tool.coverage.*]` of `pyproject.toml`. CI prints the report.
- One small synthetic file per issue code (as `make_lamda` in TEST_LAMDA_FILES): the
  check finds it; the fix removes it; the round trip stays lossless.
- Round trip: read → write → read gives equal data, for every file of the corpus.
- **Corpus** (not in the repository, opt-in with `--corpus PATH`): all LAMDA files,
  EMAA files, ProDiMo's 104 files (`/home/rab/git/prodimo/data`). Snapshot of the issue
  counts per file, so that changes in checks are visible.
- `--online` tests for BibTeX resolution, off by default (as astrolinedata's
  `--network`).

## 13. Decisions

Decided on 2026-11-18:

1. **Name**: `pylamda`.
2. **Python**: ≥ 3.11. Changed from 3.10: Sphinx ≥ 8.2 needs 3.11, and code and
   documentation share one environment.
3. **Standard**: strict LAMDA standard; partners by ID only; EMAA's extra partners
   (p-H2O, o-H2O, CO) are not supported (3.3).
4. **Patch mode**: part of the first version (sections 5, 7).
5. **P checks**: on by default as `info`; `physics=False` / `--no-physics` switches them
   off (6.1).
6. **BibTeX**: SciX for bibcodes (and for DOIs and arXiv IDs it knows); the user
   provides the SciX API token (10).
7. **Class name**: keep `LamdaData`.
8. **Git history**: not kept; the repository starts fresh.
9. **SciX token**: argument, environment variable `SCIX_API_TOKEN`, and config file
   `~/.config/pylamda/scix_token`.
10. **Database papers**: the main LAMDA and EMAA papers are always in the BibTeX output,
    from a `.bib` file shipped with the package (10, items 8 and 9).
11. **License**: MIT: permissive, OSI-approved, and the same as astrolinedata, so code
    can move between the two. (The repository started with CC0 1.0, which is not
    OSI-approved; changed on 2026-11-18.)
12. **Free-text references**:
    - resolved once, during development, into a shipped and reviewed catalogue (10.1);
    - SciX search at run time only on request, and only as candidates;
    - `fix --add-refs` / `--ref` write identifiers into files (10.2).
13. **Packaging and documentation**: Poetry; one conda environment `PYLAMDA`; a `dev`
    extra for `pip install -e ".[dev]"`; documentation on Read the Docs (2).

## 14. Steps

1. Package skeleton:
   - `pyproject.toml` (Poetry, extras), `poetry.toml`, ruff and pytest configuration;
   - Sphinx documentation skeleton, `.readthedocs.yaml`, `docs/requirements.txt`;
   - GitHub Actions for the tests.

   Then move the data model, reader and writer from astrolinedata. The reader identifies
   partners by ID, collects issues and keeps the SourceMap.
2. `verify` with all check groups (S, L, R, C, F, P). CLI `pylamda check`. Run it on the
   corpus.
3. `compare`: merge the two implementations; add critical densities.
4. `bib`: offline first, then `--online` with SciX and doi.org.
5. `fix`: normalised mode and patch mode, including `--add-refs` / `--ref`.
6. Reference catalogue: maintainer script, search and review of the corpus references,
   ship `references.json`/`.bib`; `--search` for candidates.
7. astrolinedata switches to pylamda; TEST_LAMDA_FILES switches to pylamda.
8. `cut_levels`, plots. Left for later (decided 2026-11-18, section 15).

## 15. Possible future features

Left for later (decided 2026-11-18). None of them changes the API of what exists.

- **Step 8: `cut_levels` and plots.**
  - `cut_levels(data, e_max | n_max)`: drop levels and all transitions to them, with a
    record and a note line like the other combine functions (section 9).
  - Plots with the `plot` extra (matplotlib): rate ratios b/a against E_up per
    temperature, as `TEST_LAMDA_FILES/report.py`; the same for critical densities and
    A values. The data are already in `LamdaDiff` (`rates_interpolated`,
    `critical_density`, `a_ratio`).
- **Check F04 (`pylamda check --online`)**: a DOI or bibcode in a file that does not
  resolve (warning). It can use the HTTP client and the cache of `pylamda bib`: SciX
  search for bibcodes, doi.org for DOIs. Needs the user's SciX token for bibcodes.
- From section 9, if needed: replace the rates of a partner for some transitions only;
  merge two temperature grids. These need physics decisions per case.
- Reference search: hints from the comment lines of other files of the same molecule
  (`CO_Lambda.dat` for `13CO_Lambda.dat`).
