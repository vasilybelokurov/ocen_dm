"""The gitignored .env supplies machine-specific locations without overriding the shell."""
from __future__ import annotations

import os

from ocen_dm.paths import load_local_env


def test_parses_sets_and_never_overrides(tmp_path, monkeypatch):
    for name in ("OCEN_T_A", "OCEN_T_B", "OCEN_T_C", "OCEN_T_D", "OCEN_T_KEEP"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("OCEN_T_KEEP", "from-shell")
    monkeypatch.setenv("OCEN_T_BASE", "/base")
    f = tmp_path/".env"
    home = "~"  # tilde expansion, built so no home path appears in the repository
    f.write_text(f"# comment\n\nOCEN_T_A={home}/data/x\nexport OCEN_T_B = \"/quoted path\"\nOCEN_T_C=$OCEN_T_BASE/sub\n"
                 "OCEN_T_KEEP=from-file\nnot a line\n1BAD=x\nOCEN_T_D='single'\n")
    applied = load_local_env(f)
    assert os.environ["OCEN_T_A"] == os.path.join(os.path.expanduser(home), "data", "x")
    assert os.environ["OCEN_T_B"] == "/quoted path"
    assert os.environ["OCEN_T_C"] == "/base/sub"
    assert os.environ["OCEN_T_D"] == "single"
    assert os.environ["OCEN_T_KEEP"] == "from-shell" and "OCEN_T_KEEP" not in applied
    assert set(applied) == {"OCEN_T_A", "OCEN_T_B", "OCEN_T_C", "OCEN_T_D"}


def test_missing_file_is_a_no_op(tmp_path):
    assert load_local_env(tmp_path/"absent.env") == {}
