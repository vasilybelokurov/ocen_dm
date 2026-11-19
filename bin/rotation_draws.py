#!/usr/bin/env python3
"""Plan step 5 (docs/PLAN_2026-11-19_limits.md): streaming terms for N draws of the rotation-field chain, so that the
likelihood can marginalise over the rotation field (coherent across bins and datasets) instead of adding an
independent per-bin error (tune_jam_to_data.Setup._corrected_data).

For every profile: st2_0 (MAP field, as in the data), st2_m for each draw m (same estimator: <vbar>^2 for Gaia rings,
<vbar^2> for HST/MUSE star samples), and the measurement-only errors e_lo, e_hi (before the rotation term is added).
The model prediction at draw m is exactly sqrt(pred^2 + st2_0 - st2_m), pred = the standard prediction.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/rotation_draws.py --ndraw 200
Output: results/limits/rotation_draws.npz (+ a check that the draw s.d. reproduces the production rotation error)
"""
import argparse, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
from ocen_dm.kinematics.likelihood import streaming_error_on_sigma
import jam_mock_rotation_test as J
from tune_jam_to_data import Setup


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ndraw", type=int, default=200); a = ap.parse_args()
    S = Setup(per_bin=None)
    import fit_rotation_field as F
    from ocen_dm.kinematics.hst_profile import load_hst_product
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    chain = np.load(ROOT/"results/rotation/rotation_field_cylindrical_chain_20261113.npy")
    draws = chain[np.random.default_rng(11).choice(len(chain), a.ndraw, replace=False)]
    mf0, mfs = J.MeanField(1.0), [J.MeanField(1.0, params=d) for d in draws]
    jeans, D, _ = S.fam.build(dict(zip(S.fam.names, S.x_map)))
    lum_ref = S.models(jeans, D, F.QPROJ, F.QPROJ)[0]
    prod, g = load_hst_product(), load_edr3_profile()
    out = {}
    for p, pd in zip(S.pb0.data.profiles, S.data.profiles):
        xp, yp, _, nd = S.pts[p.name]
        w = J.weights(lum_ref, xp, yp, nd)
        gaia = p.instrument == "GaiaEDR3"
        st = lambda t: t[0]**2 if gaia else t[1]
        st0 = st(S._ring_terms(p, mf0, w))
        stm = np.array([st(S._ring_terms(p, m, w)) for m in mfs])
        if gaia:
            col = "sigma_pmr" if p.kind == "pmr" else "sigma_pmt"
            e_lo = e_hi = np.asarray(g[col+"_err"], float)
        elif p.instrument == "HST" and p.kind == "pmt":
            e_lo = e_hi = np.asarray(prod["sigma_pmt_err"], float)
        else:
            e_lo, e_hi = p.err_lo, p.err_hi
        assert np.allclose(st0, pd.streaming2), p.name
        # check: draw spread vs the production per-bin rotation error
        dv = np.std(np.sqrt(stm), axis=0)
        de = streaming_error_on_sigma(pd.value, np.sqrt(st0), dv)
        prod_de = np.sqrt(np.maximum(pd.err_lo**2 - np.asarray(e_lo)**2, 0))
        print(f"{p.name}: max |st0| {np.max(st0):.3g}; rot err (draws) / production: median "
              f"{np.median(de/np.maximum(prod_de, 1e-12)) if np.any(prod_de > 0) else float('nan'):.2f}; "
              f"median rot/meas {np.median(de/np.asarray(e_lo)):.2f}")
        out[p.name+"__st0"], out[p.name+"__stm"] = st0, stm
        out[p.name+"__elo"], out[p.name+"__ehi"] = np.asarray(e_lo, float), np.asarray(e_hi, float)
        out[p.name+"__value"] = pd.value
    np.savez(ROOT/"results/limits/rotation_draws.npz", **out)
    print("saved results/limits/rotation_draws.npz")


if __name__ == "__main__":
    main()
