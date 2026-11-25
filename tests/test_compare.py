"""compare(): levels, lines, rates (exact and on the union grid), critical densities."""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest
from helpers import EMAA_CO, LAMDA

from pylamda import CollisionSet, LamdaData, compare, critical_density, match_levels, read, write
from pylamda.cli import main
from pylamda.compare import EUP_BINS, interpolate_rates, normalise_quanta, summarise_ratios, total_rates


def _co() -> LamdaData:
    """EMAA's CO with only the real p-H2 and o-H2 blocks (the first two are H2O, read by ID)."""
    d = read(EMAA_CO)
    d.collisions = d.collisions[2:4]
    return d


def _with(d: LamdaData, c: CollisionSet, **kw) -> LamdaData:
    """Copy of ``d`` with collision set ``c`` changed."""
    new = dataclasses.replace(c, **kw)
    return dataclasses.replace(d, collisions=[new if x is c else x for x in d.collisions])


def _insert_level(d: LamdaData, pos: int, energy: float, g: float, quanta: str) -> LamdaData:
    """Copy with a new level at 0-based position ``pos``; later levels are renumbered."""

    def shift(x):
        return np.where(x > pos, x + 1, x)

    return dataclasses.replace(
        d,
        level_energy=np.insert(d.level_energy, pos, energy),
        level_g=np.insert(d.level_g, pos, g),
        level_quanta=[*d.level_quanta[:pos], quanta, *d.level_quanta[pos:]],
        line_up=shift(d.line_up),
        line_low=shift(d.line_low),
        collisions=[dataclasses.replace(c, up=shift(c.up), low=shift(c.low)) for c in d.collisions],
    )


# ------------------------------------------------------------------- level matching
def test_normalise_quanta():
    assert normalise_quanta(" J=1  F=2 ") == "j=1_f=2"
    assert normalise_quanta("1_1,0; x") == "1_1_0_x"


def test_match_levels_methods():
    d = _co()
    m = match_levels(d, d)
    assert m.method == "quanta" and m.mapping == {i: i for i in range(1, 15)} and m.max_energy_difference == 0.0
    e = _insert_level(d, 5, 50.0, 2.0, "X")
    m = match_levels(d, e)
    assert m.n_matched == 14 and m.mapping[6] == 7 and 6 not in m.mapping.values()
    plain = dataclasses.replace(d, level_quanta=[""] * d.n_levels)
    assert match_levels(plain, plain).method == "energy"
    with pytest.raises(ValueError, match="quanta are missing"):
        match_levels(plain, d, match="quanta")
    assert match_levels(d, e, match="index").mapping[6] == 6
    with pytest.raises(ValueError, match="match must be"):
        match_levels(d, d, match="guess")


def test_match_levels_quanta_then_energy():
    """auto: levels whose quanta are renamed are matched by energy and g (TEST_LAMDA_FILES)."""
    d = _co()
    renamed = dataclasses.replace(d, level_quanta=[*d.level_quanta[:-1], "J=13"])
    m = match_levels(d, renamed)
    assert m.method == "quanta+energy" and m.n_matched == 14 and m.mapping[14] == 14


# ------------------------------------------------------------------- astrolinedata's comparison
def test_compare_identical_and_header():
    d = _co()
    diff = compare(d, d.copy())
    assert not diff and diff.summary() == [] and diff.level_matching == "quanta"
    other = dataclasses.replace(d, mass=29.0, dataset_doi=None)
    assert compare(d, other).header == [
        "mass: 28.0 -> 29.0",
        "dataset DOI: '10.17178/EMAA_CO_Rotation_3d0e1a5c' -> None",
    ]


def test_compare_inserted_level():
    """An inserted level (renumbering) is one added level, not a change of everything."""
    d = _co()
    e = _insert_level(d, 5, 50.0, 2.0, "X")
    diff = compare(d, e)
    assert diff.levels_added == ["X"] and not diff.levels_removed and not diff.levels_changed
    assert not diff.lines_added and not diff.lines_removed and not diff.n_a_changed
    assert not diff.partners and diff.n_levels == (14, 15)
    assert "levels: 14 -> 15 (matched by quanta)" in diff.summary()
    assert len(compare(d, e, match="index").levels_changed) > 5


def test_compare_partners_rates_and_grid():
    d = _co()
    p = d.partner("p-H2")
    rates = p.rates.copy()
    rates[0, 2] *= 1.1
    t_new = np.append(p.temperatures, p.temperatures[-1] + 100.0)
    rates = np.hstack([rates, rates[:, -1:]])
    new_p = dataclasses.replace(p, temperatures=t_new, rates=rates)
    h = CollisionSet(5, np.array([10.0, 100.0]), np.array([2]), np.array([1]), np.array([[1e-10, 2e-10]]))
    e = dataclasses.replace(d, collisions=[new_p, h])  # o-H2 removed, H added
    diff = compare(d, e)
    assert diff.partners_added == ["H"] and diff.partners_removed == ["o-H2"]
    pd = diff.partners["p-H2"]
    assert pd.temperatures_added == [t_new[-1]] and not pd.temperatures_removed
    assert pd.n_rates_changed == 1 and pd.max_rate_change == pytest.approx(0.1)
    assert pd.max_rate_change_at == f"1 -> 0 at {p.temperatures[2]:g} K"
    text = "\n".join(diff.summary())
    assert "new collision partner H" in text and "collision partner o-H2 removed" in text
    t2 = f"{p.temperatures[2]:g}"
    assert (
        f"1 rate in 1 transition differs by more than 1e-06 (relative), max change B/A 1.1 (1 -> 0 at {t2} K)" in text
    )
    assert "in 1 transition:" not in text  # no list of transitions, only the largest change
    assert pd.changed_transitions == [
        {
            "transition": "1 -> 0",
            "ratio": pytest.approx(1.1),
            "temperature": p.temperatures[2],
            "rate_a": p.rates[0, 2],
            "rate_b": pytest.approx(1.1 * p.rates[0, 2]),
            "n_changed": 1,
        }
    ]
    assert pd.ratio["median"] == 1.0 and len(pd.ratio_by_temperature) == len(p.temperatures)
    assert diff.as_dict()["partners"]["p-H2"]["n_rates_changed"] == 1
    # the union grid is the common range: the added temperature lies outside a's grid
    assert pd.interp_temperatures == [float(t) for t in p.temperatures]


def test_compare_lines_and_energy_matching():
    d = _co()
    a = d.line_a.copy()
    a[3] *= 2
    diff = compare(d, dataclasses.replace(d, line_a=a))
    assert diff.n_a_changed == 1 and diff.max_a_change == pytest.approx(1.0)
    assert diff.a_ratio.n_changed == 1 and diff.a_ratio.max_at == "4 -> 3"
    assert any(
        line.startswith("  A by E_up of the upper level (changed by > 10%): 0-100 K 1 of ") for line in diff.summary()
    )
    plain = dataclasses.replace(d, level_quanta=[""] * d.n_levels)
    shifted = dataclasses.replace(plain, level_energy=plain.level_energy * (1 + 5e-6))
    diff = compare(plain, shifted)  # within the matching tolerance, beyond rtol
    assert diff.level_matching == "energy" and not diff.levels_added and not diff.lines_removed
    assert len(diff.levels_changed) == d.n_levels - 1  # all but E = 0


def test_compare_statistics_and_decrease():
    """Ratio statistics per temperature; a decrease is a factor a/b; zero rates are skipped."""
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    rates = p.rates.copy()
    rates[:, 0] *= 0.5
    rates[1, 1] = 0.0
    pd = compare(d, _with(d, p, rates=rates), rtol=0.1).partners["p-H2"]
    first, second = pd.ratio_by_temperature[:2]
    assert first["median"] == pytest.approx(0.5) and first["max_factor"] == pytest.approx(2.0)
    assert second["median"] == 1.0
    assert pd.n_rates_changed == p.n_transitions + 1 and pd.max_rate_factor == pytest.approx(2.0)
    assert all(c["ratio"] == pytest.approx(0.5) for c in pd.changed_transitions if c["n_changed"] == 1)
    # the summary shows only the temperature whose median is farthest from 1
    lines = compare(d, _with(d, p, rates=rates), rtol=0.1).summary()
    t0 = f"{p.temperatures[0]:g}"
    assert f"  largest deviation by temperature at {t0} K: median 0.5, 90% within 0.5-0.5" in lines
    assert not any("by temperature (median" in line for line in lines)


def test_compare_zero_rates_summary():
    """Old rates all 0: no ratios exist; the summary must still work."""
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    diff = compare(_with(d, p, rates=p.rates * 0.0), d)
    pd = diff.partners["p-H2"]
    assert pd.n_rates_changed == p.rates.size and pd.ratio == {} and pd.ratio_by_temperature == []
    assert pd.rates_interpolated is None and pd.total_rate is None
    text = "\n".join(diff.summary())
    assert "no pair with both rates > 0" in text and "all to or from a zero rate" in text
    assert f"differ: all to or from a zero rate (1 -> 0: 0 -> {p.rates[0, 0]:.3g} at {p.temperatures[0]:g} K)" in text


def test_compare_mixed_zero_and_ratio_changes():
    """One rate appears (0 -> value), one is tripled, one quartered: one consistent maximum."""
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    a_rates, b_rates = p.rates.copy(), p.rates.copy()
    a_rates[5, 3] = 0.0
    b_rates[0, 0] *= 3.0
    b_rates[2, 1] *= 0.25
    diff = compare(_with(d, p, rates=a_rates), _with(d, p, rates=b_rates), rtol=0.1)
    pd = diff.partners["p-H2"]
    assert (pd.n_rates_changed, pd.n_rates_zero_changed) == (3, 1)
    assert pd.max_rate_factor == pytest.approx(4.0) and pd.max_rate_change == pytest.approx(0.75)
    assert pd.max_rate_ratio == pytest.approx(0.25)  # the same rate, as B/A
    assert [c["ratio"] for c in pd.changed_transitions] == [None, pytest.approx(0.25), pytest.approx(3.0)]
    text = "\n".join(diff.summary())
    assert "3 rates in 3 transitions differ by more than 0.1 (relative), max change B/A 0.25" in text
    assert "; 1 to or from a zero rate (" in text  # with one example


def test_compare_reference_changed():
    d = _co()
    p = d.partner("p-H2")
    diff = compare(d, _with(d, p, reference="another"))
    assert diff.partners["p-H2"].reference_changed and "  reference changed" in diff.summary()


# ------------------------------------------------------------------- the merged parts
def test_interpolate_rates():
    t = np.array([10.0, 100.0])
    r = np.array([[1e-11, 1e-10], [0.0, 1e-10]])
    out = interpolate_rates(t, r, np.array([10.0, np.sqrt(1000.0), 100.0, 1000.0]))
    np.testing.assert_allclose(out[0], [1e-11, np.sqrt(1e-21), 1e-10, 1e-10])  # log-log; clamped outside
    assert out[1, 0] == 0.0 and out[1, 2] == pytest.approx(1e-10)
    assert interpolate_rates(np.array([50.0]), np.array([[2e-11]]), np.array([10.0, 99.0]))[0].tolist() == [
        2e-11,
        2e-11,
    ]


def test_union_grid_when_grids_differ():
    """Different grids: no common temperature, but the interpolated comparison works."""
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    t = p.temperatures
    mid = np.sqrt(t[:-1] * t[1:])  # geometric midpoints: none in common with t
    b_rates = interpolate_rates(t, p.rates, mid) * 2.0
    pd = compare(d, _with(d, p, temperatures=mid, rates=b_rates)).partners["p-H2"]
    assert pd.n_temperatures_compared == 0 and pd.n_rates_compared == 0
    assert pd.interp_temperatures[0] == mid[0] and pd.interp_temperatures[-1] == mid[-1]
    s = pd.rates_interpolated
    assert s.median == pytest.approx(2.0, rel=0.05) and s.n_changed == p.n_transitions
    text = compare(d, _with(d, p, temperatures=mid, rates=b_rates)).summary()
    assert any("rates B/A on the union grid" in line for line in text)


def test_critical_density_and_total_rate():
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    k = total_rates(d, p)
    assert k.shape == (d.n_levels, len(p.temperatures)) and np.all(k[0] == 0)
    np.testing.assert_allclose(k[1], p.rates[(p.up == 2)].sum(axis=0))
    n = critical_density(d, "p-H2")
    a2 = d.line_a[d.line_up == 2].sum()
    np.testing.assert_allclose(n[1], a2 / k[1])
    assert np.all(np.isnan(n[0]))  # the ground level has no radiative decay
    n100 = critical_density(d, p, np.array([100.0]))
    assert n100.shape == (d.n_levels, 1)


def test_compare_added_transition_changes_total_rate():
    """A new transition out of a level changes its total rate, though no common rate changed."""
    d = _co().select_partners(["p-H2"])
    p = d.partner("p-H2")
    rows = np.flatnonzero(p.up == 3)
    keep = np.setdiff1d(np.arange(p.n_transitions), rows[:1])  # a lacks one row out of level 3
    a = _with(d, p, up=p.up[keep], low=p.low[keep], rates=p.rates[keep])
    pd = compare(a, d).partners["p-H2"]
    assert pd.transitions_added == ["2 -> 0"] and pd.n_rates_changed == 0
    k = pd.total_rate
    assert k.n_changed >= 1 and k.max_at.startswith("level 2 at") and k.max_ratio > 1
    assert k.by_eup[0]["eup_min"] == 0.0 and k.by_eup[0]["n_changed"] >= 1
    lines = compare(a, d).summary()
    assert any("total de-excitation rate B/A" in line for line in lines)
    assert not any("critical density" in line for line in lines)
    assert not hasattr(pd, "critical_density")


def test_summarise_ratios_by_eup():
    q = np.array([[1.0, 1.0], [1.5, 1.0], [np.nan, 0.5], [1.0, 1.0]])
    s = summarise_ratios(q, ["a", "b", "c", "d"], np.array([50.0, 150.0, 5000.0, 20000.0]), [10.0, 20.0], 0.1)
    assert (s.n, s.n_items, s.n_changed) == (7, 4, 2)
    assert s.max_factor == pytest.approx(2.0) and s.max_ratio == 0.5 and s.max_at == "c at 20 K"
    assert [(x["eup_min"], x["n_items"], x["n_changed"]) for x in s.by_eup] == [
        (0.0, 1, 0),
        (100.0, 1, 1),
        (3000.0, 1, 1),
        (10000.0, 1, 0),
    ]
    assert (
        s.eup_text()
        == "0-100 K 0 of 1, 100-300 K 1 of 1 (max change B/A 1.5), 3000-10000 K 1 of 1 (max change B/A 0.5), "
        "10000+ K 0 of 1"
    )
    assert summarise_ratios(np.full((2, 1), np.nan), ["a", "b"], np.zeros(2), None, 0.1) is None
    assert EUP_BINS[0] == 0.0


def test_repeated_partner_ids():
    """EMAA's file has two blocks with ID 2 and two with ID 3: compared in order."""
    d = read(EMAA_CO)
    p2 = d.collisions[2]
    diff = compare(d, _with(d, p2, rates=p2.rates * 2))
    assert list(diff.partners) == ["p-H2 #2"]
    assert not compare(d, d.copy())


def test_files_and_level_errors(tmp_path):
    write(_co(), tmp_path / "a.dat")
    assert not compare(tmp_path / "a.dat", str(tmp_path / "a.dat"))
    d = _co()
    d.line_up[0] = 99
    with pytest.raises(ValueError, match="b: radiative transitions: 1 transition"):
        compare(_co(), d)


def test_as_dict_is_json():
    d = read(LAMDA / "catom.dat")
    e = d.copy()
    e.collisions[0].rates *= 3
    out = compare(d, e).as_dict()
    text = json.dumps(out, allow_nan=False)  # no inf/nan
    assert '"eup_max": null' in text or '"eup_max"' in text
    assert out["partners"]["H"]["total_rate"]["n_changed"] >= 1 and "critical_density" not in out["partners"]["H"]


# ------------------------------------------------------------------- command line
def test_compare_command(tmp_path, capsys):
    d = _co()
    write(d, tmp_path / "a.dat")
    d.partner("p-H2").rates[0, 0] *= 2
    write(d, tmp_path / "b.dat")
    assert main(["compare", str(tmp_path / "a.dat"), str(tmp_path / "a.dat")]) == 0
    assert capsys.readouterr().out == "no differences\n"
    assert main(["compare", str(tmp_path / "a.dat"), str(tmp_path / "b.dat"), "--rtol", "0.5"]) == 1
    out = capsys.readouterr().out
    assert "1 rate in 1 transition differs by more than 0.5 (relative), max change B/A 2" in out
    assert "  by E_up of the upper level (transitions changed by > 10%): 0-100 K 1 of " in out
    assert main(["compare", str(tmp_path / "a.dat"), str(tmp_path / "b.dat"), "--format", "json"]) == 1
    data = json.loads(capsys.readouterr().out)
    assert data["partners"]["p-H2"]["n_rates_changed"] == 1 and data["rtol"] == 1e-6
    assert main(["compare", str(tmp_path / "a.dat"), str(tmp_path / "missing.dat")]) == 2
    (tmp_path / "bad.dat").write_text("!MOLECULE\nCO\n")
    assert main(["compare", str(tmp_path / "a.dat"), str(tmp_path / "bad.dat")]) == 2
    assert "pylamda compare: error:" in capsys.readouterr().err
