"""bin/sph_profile_limit.upper_limit and the NFW M(<50 pc) conversion."""
import sys
from pathlib import Path
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"bin")); sys.path.insert(0, str(ROOT/"src"))
from scipy.integrate import quad


def test_m50_matches_density_integral():
    from sph_profile_limit import m50
    rs = 300.0; rho = lambda r: 1/((r/rs)*(1+r/rs)**2)
    M = lambda R: quad(lambda r: 4*np.pi*r*r*rho(r), 0, R)[0]
    assert m50(1e6, rs) == pytest.approx(1e6*M(50)/M(100), rel=1e-8)


def test_upper_limit_parabola():
    from sph_profile_limit import upper_limit, LOGM
    # q = ((lm - 5)/0.5)^2: q = 2.71 at lm = 5 + 0.5*sqrt(2.71)
    rows = [dict(logM100=float(l), M50=10**l, lnpost=-0.5*((l-5)/0.5)**2) for l in LOGM]
    m100, m50v = upper_limit(rows, 0.0, 2.71)
    assert np.log10(m100) == pytest.approx(5 + 0.5*np.sqrt(2.71), abs=0.03)
