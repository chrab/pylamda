"""The command line: pylamda check."""

from __future__ import annotations

import json
import subprocess
import sys

import pytest
from helpers import make_lamda

from pylamda import CODES
from pylamda.cli import main


@pytest.fixture
def files(tmp_path):
    clean = tmp_path / "clean.dat"
    clean.write_text(make_lamda())
    warning = tmp_path / "warning.dat"
    warning.write_text(make_lamda(partner="1 X + p-H2 ! https://doi.org/10.1/x"))  # C02
    error = tmp_path / "error.dat"
    error.write_text(make_lamda(temps="20.0 10.0"))  # C05
    broken = tmp_path / "broken.dat"
    broken.write_text(make_lamda(ncol_comment="-"))  # a LAMDA header, then a count that does not read
    other = tmp_path / "other.txt"
    other.write_text("!MOLECULE\nX\n")  # ends before the molecular weight: not a LAMDA file
    return {"clean": clean, "warning": warning, "error": error, "broken": broken, "other": other}


def test_exit_codes(files, capsys):
    assert main(["check", str(files["clean"])]) == 0
    assert capsys.readouterr().out == "Found 0 issues in 1 file.\n"
    assert main(["check", str(files["clean"]), str(files["warning"])]) == 0
    assert main(["check", str(files["warning"]), str(files["error"])]) == 1
    assert main(["check", str(files["error"]), str(files["broken"])]) == 2
    assert main(["check", str(files["clean"]), str(files["other"])]) == 2
    assert main(["check", str(files["clean"]), str(files["broken"].with_name("missing.dat"))]) == 2
    out = capsys.readouterr().out
    assert "missing.dat: cannot open" in out


def test_text_output(files, capsys):
    main(["check", str(files["warning"]), str(files["broken"]), str(files["other"])])
    out = capsys.readouterr().out.splitlines()
    assert out[0].startswith(f"{files['warning']}:19: C02 warning: partner line '1 X + p-H2")
    assert out[1].startswith(f"{files['broken']}:20: S06 error: cannot read the number of collisional transitions")
    assert out[2] == f"{files['broken']}: not readable"
    # one statement for a file that is not LAMDA, without its structure issues
    assert out[3] == f"{files['other']}: not a LAMDA file: line 2: the file ends before the molecular weight"
    assert out[4] == "Found 2 issues (1 error, 1 warning) in 3 files; 1 not a LAMDA file; 1 not readable."


def test_options(files, capsys, tmp_path):
    assert main(["check", "--error-on", "C02", str(files["warning"])]) == 1
    assert "C02 error:" in capsys.readouterr().out
    assert main(["check", "--ignore", "C", str(files["warning"]), str(files["error"])]) == 0
    assert main(["check", "--select", "F,P", "--select", "C05", str(files["error"])]) == 1
    capsys.readouterr()
    scaled = tmp_path / "scaled.dat"
    scaled.write_text(make_lamda(partner="1 X + H2 ! scaled from He, https://doi.org/10.1/x"))  # P02
    main(["check", str(scaled)])
    assert " P02 " not in capsys.readouterr().out  # the physics checks are off by default
    main(["check", "--physics", str(scaled)])
    assert " P02 " in capsys.readouterr().out
    main(["check", "--select", "P", str(scaled)])
    assert " P02 " in capsys.readouterr().out
    main(["check", "--select", "P", "--no-physics", str(scaled)])
    assert " P02 " not in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["check", "--physics", "--no-physics", str(scaled)])


def test_text_output_per_file(files, capsys, monkeypatch):
    """Text output: the issues of each file are printed before the next file is checked."""
    import pylamda.cli as cli

    seen = []
    real_verify = cli.verify

    def spy(path, **kw):
        seen.append(capsys.readouterr().out)  # what was printed before this file is checked
        return real_verify(path, **kw)

    monkeypatch.setattr(cli, "verify", spy)
    main(["check", str(files["warning"]), str(files["other"]), str(files["clean"])])
    assert seen[0] == ""
    assert "C02 warning" in seen[1]
    assert "not a LAMDA file" in seen[2]
    assert capsys.readouterr().out.startswith("Found 1 issue (1 warning) in 3 files")


def test_json_output(files, capsys):
    args = ["check", "--format", "json", str(files["warning"]), str(files["broken"]), str(files["other"])]
    assert main(args) == 2
    out = json.loads(capsys.readouterr().out)
    assert [(f["path"], f["readable"], f["error"]) for f in out] == [
        (str(files["warning"]), True, None),
        (str(files["broken"]), False, "not readable"),
        (str(files["other"]), False, "not a LAMDA file: line 2: the file ends before the molecular weight"),
    ]
    assert out[0]["issues"][0]["code"] == "C02" and out[0]["issues"][0]["line"] == 19
    assert out[1]["issues"][0]["code"] == "S06" and out[2]["issues"] == []


def test_statistics(files, capsys):
    main(["check", "--statistics", str(files["warning"]), str(files["error"]), str(files["clean"])])
    out = capsys.readouterr().out.splitlines()
    assert out[0].split()[:4] == ["1", "C02", "warning", "in"]
    assert out[1].split()[:4] == ["1", "C05", "error", "in"]
    assert out[-1] == "Found 2 issues (1 error, 1 warning) in 3 files."


def test_list_codes(capsys):
    assert main(["check", "--list-codes"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == len(CODES) and lines[0].startswith("S01  error")


def test_no_files(capsys):
    assert main(["check"]) == 2
    assert "no files given" in capsys.readouterr().err


def test_python_m(files):
    result = subprocess.run(
        [sys.executable, "-m", "pylamda", "check", str(files["error"])], capture_output=True, text=True
    )
    assert result.returncode == 1 and "C05 error" in result.stdout
