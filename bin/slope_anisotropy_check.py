#!/usr/bin/env python3
"""Necessary DF-consistency conditions for spherical models with a given beta(r) (2026-11-19, plan step 1).

(i)  Global density slope-anisotropy inequality gamma(r) >= 2 beta(r), gamma = -dln nu/dln r of the tracer
     (Ciotti & Morganti 2010, https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.1070C -- proven for a broad class of
     separable augmented densities; a necessary condition, not sufficient).
(ii) Central point-mass limit (An & Evans 2006, https://ui.adsabs.harvard.edu/abs/2006ApJ...642..752A):
     gamma_0 >= beta_0 + 1/2 as r -> 0. Reported also as a local diagnostic gamma - beta - 1/2 where the BH holds
     > 50% of the enclosed mass (not a theorem away from r -> 0).
Usage: PYTHONPATH=src:bin python bin/slope_anisotropy_check.py   (spherical local optima, results/limits/sph_local_*_inner.json)
Module: check(model, D) -> dict for a SphericalJeans model.
"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))


def check(model, D, bh_mass, beta=None):
    pc = D*1e3/206264.806
    r = np.geomspace(0.02*pc, 3000*pc, 600)
    nu = np.asarray(model.tracer.density(r), float)
    gam = -np.gradient(np.log(nu), np.log(r))
    b = np.asarray((beta or model.anisotropy.beta)(r), float)
    menc = np.asarray(model.mass.enclosed_mass(r), float)
    bhdom = bh_mass/menc > 0.5
    m1 = gam - 2*b
    m2 = gam - b - 0.5
    i1 = int(np.argmin(m1))
    out = dict(min_gamma_minus_2beta=float(m1[i1]), at_arcsec=float(r[i1]/pc), passes_CM=bool(m1.min() >= -1e-3),
               central_gamma0=float(gam[0]), central_beta0=float(b[0]), passes_AE_central=bool(m2[0] >= -1e-2),
               bh_dominated_to_arcsec=float(r[bhdom][-1]/pc) if bhdom.any() else 0.0,
               min_local_AE_in_bh_region=float(m2[bhdom].min()) if bhdom.any() else None,
               at_arcsec_local=float(r[bhdom][np.argmin(m2[bhdom])]/pc) if bhdom.any() else None)
    return out


if __name__ == "__main__":
    from sph_local import Local
    res = {}
    for kind in ("K1", "NFW"):
        for mode in ("inner", "legacy"):
            f = ROOT/f"results/limits/sph_local_{kind}_{mode}.json"
            if not f.exists():
                continue
            J = json.loads(f.read_text()); th = J["best"]["theta"]
            L = Local(kind, inner_node=2.0 if mode == "inner" else None)
            model, D, _ = L.fam.build(th)
            res[f"{kind}_{mode}"] = c = check(model, D, th["M_bh"])
            print(kind, mode, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in c.items()}))
    (ROOT/"results/limits/slope_anisotropy_spherical.json").write_text(json.dumps(res, indent=1))
