"""Sensitivity prior options of bin/sph_local.Local (plan step 6): M_rem/M_star half-normal + cap, M_star floor."""
import json, sys
from pathlib import Path
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"bin")); sys.path.insert(0, str(ROOT/"src"))
SRC = ROOT/"results/limits/sph_local_K1_inner.json"
pytestmark = pytest.mark.skipif(not SRC.exists(), reason="needs the local K1 inner-node fit")


def test_mrem_halfnormal_and_floor():
    from sph_local import Local
    th = json.loads(SRC.read_text())["best"]["theta"]
    L0 = Local("K1"); L = Local("K1", mrem_hn=0.1, prior_overrides={"M_star": (1.3e6, 1e7)})
    x = np.array([th[n] for n in L.names]); x[L.names.index("M_rem")] = 0.2*x[L.names.index("M_star")]
    l0, p0 = L0.terms(x); l1, p1 = L.terms(x)
    assert l0 == pytest.approx(l1) and p1 - p0 == pytest.approx(-0.5*(0.2/0.1)**2)
    x[L.names.index("M_rem")] = 0.31*x[L.names.index("M_star")]
    assert L.terms(x)[1] == -np.inf
    i = L.names.index("M_star"); assert L.lo[i] == pytest.approx(np.log10(1.3e6))
