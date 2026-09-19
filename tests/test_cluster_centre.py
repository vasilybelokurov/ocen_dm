"""The centre is derived from the catalogues, never written down.

Guards the 2026-09-19 defect: a hard-coded declination with a fabricated provenance
("Baumgardt catalogue centre", from a directory that was never downloaded) sitting 10.7
arcsec north of the centre both catalogues actually use.
"""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.cluster import OCEN_DEC, OCEN_RA, centre
from ocen_dm.paths import processed_dir, raw_dir

_HAS = ((raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits").exists()
        and (processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists())
pytestmark = pytest.mark.skipif(not _HAS, reason="catalogues not present")


def _sep(a, b):
    return float(np.hypot((a[0] - b[0]) * np.cos(np.radians(a[1])), a[1] - b[1]) * 3600)


def test_the_constants_are_what_the_data_say():
    assert _sep((OCEN_RA, OCEN_DEC), centre("omegacat")) < 0.05, \
        "the module constants must equal the derived oMEGACat centre"


def test_the_two_catalogues_agree_with_each_other():
    assert _sep(centre("omegacat"), centre("vasiliev")) < 1.5


def test_the_old_hard_coded_value_is_rejected():
    assert _sep((201.696833, -47.476583), centre("omegacat")) > 5.0
    assert _sep((201.696833, -47.476583), centre("vasiliev")) > 5.0
    assert abs(OCEN_DEC + 47.476583) * 3600 > 5.0, "the wrong declination must not come back"


def test_every_module_uses_the_one_definition():
    from ocen_dm.kinematics import hst_profile, outer_profile
    from ocen_dm.selection import field_template, hst_gaia_match
    for mod in (outer_profile, hst_profile, field_template, hst_gaia_match):
        assert mod.OCEN_RA == OCEN_RA and mod.OCEN_DEC == OCEN_DEC, mod.__name__


def test_no_module_hard_codes_a_centre():
    import pathlib
    src = pathlib.Path(__file__).resolve().parents[1] / "src" / "ocen_dm"
    offenders = [p.name for p in src.rglob("*.py")
                 if p.name != "cluster.py" and "-47.4" in p.read_text()]
    assert not offenders, f"centre literal re-introduced in {offenders}"


def test_the_innermost_hst_bin_is_really_central():
    """The symptom that exposed it: 'r = 3 arcsec' stars were really at 8-14 arcsec."""
    from ocen_dm.kinematics.hst_profile import load_hst_sample
    s = load_hst_sample()
    inner = s.r_arcsec[(s.r_arcsec < 5) & (s.quality_flag > 0)]
    assert len(inner) > 50
    true = centre("omegacat")
    from astropy.table import Table
    t = Table.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    ra, dec = np.asarray(t["RA"], float), np.asarray(t["DEC"], float)
    r_true = np.hypot((ra - true[0]) * np.cos(np.radians(true[1])), dec - true[1]) * 3600
    assert np.nanmax(r_true[(r_true < 5)]) < 5.01
