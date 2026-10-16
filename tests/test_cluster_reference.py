"""The cluster's reference values live only in ``ocen_dm.cluster``; nobody types them in.

Guards two defects of the same kind: a hard-coded centre with a fabricated provenance, 10.7
arcsec off (2026-09-19), and a systemic proper motion labelled "Vasiliev & Baumgardt 2021"
that matches no source, copied into five modules (2026-10-16). Every numeric literal in
``src/`` and ``bin/`` (comments and strings excluded) is compared with every version of omega
Cen's position, motion and distance seen in the repository or its sources. Fitted values
(stream best-fit PMs, distance grids) are different numbers and are not affected.
"""
from __future__ import annotations

import io
import tokenize
from pathlib import Path

import pytest

from ocen_dm import cluster

ROOT = Path(__file__).resolve().parents[1]

#: every version of omega Cen's reference values met so far (absolute values; signs vary by convention)
BANNED = {
    "RA (deg)": {201.69683333, 201.69683, 201.696833, 201.696838, 201.696991, 201.69699, 201.697, 201.6970},
    "Dec (deg)": {47.47956944, 47.47958, 47.479569, 47.479339, 47.479472, 47.479473, 47.47947, 47.48, 47.476583},
    "l, b (deg)": {309.102, 14.968},
    "systemic PM (mas/yr)": {3.250, 3.2499, 3.257, 3.236, 6.746, 6.7461, 6.730, 6.731},
    "systemic v_los (km/s)": {232.7, 232.78, 232.6, 232.99},
    "distance (kpc)": {5.43, 5.494},
}


def _offenders(path):
    found = []
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(path.read_text()).readline))
    except (tokenize.TokenError, SyntaxError, UnicodeDecodeError):
        return found
    for tok in tokens:
        if tok.type != tokenize.NUMBER:
            continue
        try:
            value = abs(float(tok.string))
        except ValueError:
            continue
        for kind, values in BANNED.items():
            if any(value == v for v in values):
                found.append(f"{path.relative_to(ROOT)}:{tok.start[0]}: {tok.string} ({kind})")
    return found


def test_no_cluster_reference_value_is_typed_in():
    files = [p for d in ("src", "bin") for p in sorted((ROOT/d).rglob("*.py"))
             if p.name != "cluster.py" and "__pycache__" not in p.parts]
    offenders = [o for p in files for o in _offenders(p)]
    assert not offenders, "import these from ocen_dm.cluster instead:\n" + "\n".join(offenders)


def test_the_banned_set_covers_the_adopted_values():
    adopted = [cluster.OCEN_RA, cluster.OCEN_DEC, *cluster.MU_SYS, cluster.OCEN_VSYS_KMS,
               cluster.OCEN_DISTANCE_KPC, cluster.OCEN_DISTANCE_OMEGACAT_KPC]
    every = set().union(*BANNED.values())
    assert all(abs(v) in every for v in adopted)


def test_galactic_coordinates_of_the_centre():
    l, b = cluster.galactic()
    assert l == pytest.approx(309.102, abs=2e-3) and b == pytest.approx(14.968, abs=2e-3)
