"""Literature presets: assumptions map onto the family correctly."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics.presets import PRESETS, build_preset, derived_quantities, v_band_luminosity
from ocen_dm.light_model import MGEFit


def test_v_band_luminosity_of_omega_cen():
    # Harris V_t = 3.68, E(B-V) = 0.12 -> A_V = 0.372; at 4.59 kpc: m - M = 13.31; M_V = -10.00
    L = v_band_luminosity(4.59)
    assert 8e5 < L < 9.5e5                                  # 10**(-0.4*(-10.00-4.83)) = 8.6e5


def test_watkins_preset_has_constant_beta_no_remnants_no_bh_fixed_distance():
    p, fam, datasets = build_preset("watkins2013")
    assert datasets == ["hst_pm_radial", "hst_pm_tangential"]
    assert fam.names == ("M_star", "beta_0")
    theta = fam.to_dict(fam.transform(np.array([0.5, 0.9])))
    full = fam.complete(theta)
    assert full["M_rem"] == 1e4 and full["M_bh"] == 1e2 and full["beta_inf"] == full["beta_0"]
    assert fam.distance(full) == 4.59
    assert 0.3 < theta["beta_0"] < 0.5                      # beta0_max raised to 0.5 for this preset
    d = derived_quantities(p, theta, fam)
    assert d["ML_V"] == pytest.approx(theta["M_star"] / v_band_luminosity(4.59))


def test_all_presets_build():
    for name in PRESETS:
        p, fam, datasets = build_preset(name)
        assert len(fam.names) >= 2 and all(isinstance(d, str) for d in datasets)
        assert p.url.startswith("https://ui.adsabs.harvard.edu/")
    with pytest.raises(KeyError):
        build_preset("king1966")
