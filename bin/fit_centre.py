#!/usr/bin/env python3
"""Joint centre of omega Cen on the oMEGACat VI catalogue (2026-10-20; src/ocen_dm/centre_fit.py).

Stages (offsets from the adopted AvdM10 centre, arcsec; x = dRA* east +, y = dDec):
  A  symmetry centres (fixed truncation), F625W 19-22, annuli 5-55", 60-120", 120-200";
  B  counts only, F625W 19-22, R < 60, 120, 200": with the depth-tied completeness term; variants without it, with cells of
     |D| > 0.05 masked, and with a lopsided (m=1, amplitude ~ R) term (R < 200");
  C  dispersions only (hq astrometry + membership), R < 200": all magnitudes; F625W < 20.5 (small PM errors);
  D  joint (B + C, R < 200", shared centre).
Errors: numerical Hessian of -lnL over all free parameters (statistical only; calibrated by bin/centre_mocks.py).
Output: results/plot_data/centre_fit.json, plots/centre_fit.png
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src python bin/fit_centre.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from ocen_dm.centre_fit import CountModel, DepthMap, DispersionModel, fit, load_omegacat, symmetry_centre   # noqa: E402

REFS = {"Goldsbury+2010": (0.41, 0.35), "oMEGACat VI kinematic": (-1.31, 1.44)}


def row(lab, r, out):
    out[lab] = dict(centre=r.centre, err=r.centre_err, loglike=r.loglike, converged=r.converged, message=r.message, params=r.params)
    print(f"{lab:38s} ({r.centre[0]:+.2f} +- {r.centre_err[0]:.2f}, {r.centre[1]:+.2f} +- {r.centre_err[1]:.2f})  "
          f"converged {r.converged}  [{time.strftime('%H:%M:%S')}]", flush=True)


def main():
    s = load_omegacat(); depth = DepthMap(s); out = dict(symmetry={}, fits={})
    mid = (s.mag >= 19) & (s.mag < 22); sm = s.subset(mid)
    print(f"stars with F625W in R < 210\": {len(s.x)}; F625W 19-22: {mid.sum()}; hq: {s.hq.sum()}; depth D range "
          f"{np.percentile(depth.D[depth.D != 0], [5, 95]).round(3)}", flush=True)
    for a, b in ((5., 55.), (60., 120.), (120., 200.)):
        c, edge = symmetry_centre(sm.x, sm.y, a, b)
        out["symmetry"][f"{a:g}-{b:g}"] = dict(centre=c.tolist(), at_grid_edge=edge)
        print(f"symmetry {a:g}-{b:g}\": ({c[0]:+.2f}, {c[1]:+.2f}) edge {edge}", flush=True)
    for rmax in (60., 120., 200.):
        row(f"counts 19-22 R<{rmax:g} depth", fit(counts=CountModel(sm, depth, rmax, grid=0.5 if rmax <= 60 else 1.0)), out["fits"])
    flat = lambda x, y: np.zeros(np.shape(x))
    row("counts 19-22 R<200 no depth term", fit(counts=CountModel(sm, flat, 200.)), out["fits"])
    row("counts 19-22 R<200 |D|>0.05 masked", fit(counts=CountModel(sm, depth, 200., mask_depth=0.05)), out["fits"])
    row("counts 19-22 R<200 depth + lopsided", fit(counts=CountModel(sm, depth, 200., lopsided=True)), out["fits"])
    row("disp all R<200", fit(disp=DispersionModel(s, depth, 200.)), out["fits"])
    row("disp F625W<20.5 R<200", fit(disp=DispersionModel(s, depth, 200., mag=(-99, 20.5))), out["fits"])
    row("JOINT counts 19-22 + disp all, R<200", fit(counts=CountModel(sm, depth, 200.), disp=DispersionModel(s, depth, 200.)), out["fits"])
    (ROOT/"results/plot_data/centre_fit.json").write_text(json.dumps(out, indent=1, default=float))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 7.5), constrained_layout=True)
    for i, (lab, r) in enumerate(out["fits"].items()):
        ax.errorbar(*r["centre"], xerr=r["err"][0], yerr=r["err"][1], fmt="*" if lab.startswith("JOINT") else "o", ms=12 if lab.startswith("JOINT") else 5,
                    capsize=2, color=f"C{i % 10}", label=lab)
    for lab, r in out["symmetry"].items():
        ax.plot(*r["centre"], "^", color="0.4", ms=6, label=f"symmetry {lab}\"")
    ax.plot(0, 0, "k+", ms=16, mew=2, label="AvdM10 (adopted)")
    for j, (lab, o) in enumerate(REFS.items()):
        ax.plot(*o, "x", ms=11, mew=2, color=["m", "r"][j], label=lab)
    ax.set(xlim=(4, -4), ylim=(-4, 4), aspect="equal", xlabel="dRA* [arcsec] (east left)", ylabel="dDec [arcsec]",
           title="omega Cen centre: counts (F625W 19-22, depth-tied completeness), PM dispersions, joint (Hessian errors)")
    ax.axhline(0, color="0.85", lw=0.5); ax.axvline(0, color="0.85", lw=0.5); ax.legend(fontsize=6.5, loc="lower left")
    fig.savefig(ROOT/"plots/centre_fit.png", dpi=110); print(ROOT/"plots/centre_fit.png")


if __name__ == "__main__":
    main()
