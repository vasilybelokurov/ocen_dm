#!/usr/bin/env python3
"""Density centre of omega Cen from the oMEGACat star positions, by magnitude (2026-10-20).

Question: is the bright-star (F625W < 19) central dip about the adopted Anderson & van der Marel (2010) centre a centring
effect? Do bright and faint stars define the same centre?
Data: oMEGACat VI catalogue, all stars with measured F625W (the C1 count selection, no quality flag); offsets
x = dRA* (east positive) = -0.04"(x_pix - 15000), y = dDec = 0.04"(y_pix - 15000) from the adopted centre (pixel 15000,
15000 = AvdM10). The HST coverage is complete inside these apertures (C1 coverage = 1 inside 30", >= 0.956 inside 250").
Slices: F625W < 19, 19-21, 21-23, and all < 23.
Methods (the iterative centroid was dropped: in the flat core it has no restoring force and wanders by several arcsec):
  1. symmetry (AvdM10 "pie slices"): the trial centre that minimises the chi^2 of the count differences between opposite
     wedges (8 wedge pairs of 22.5 deg, annulus 5-55"), on a 0.5" grid +-4" then 0.1" +-0.6", refined by a parabola;
  2. model fit: unbinned Poisson likelihood of an elliptical cored power law Sigma ~ (1 + m^2/a^2)^(-b), m^2 = x'^2 + y'^2/q^2,
     within R < 60, 120, 200" of the adopted centre (normalisation on a grid of 0.5" (60") or 1" (120, 200")), free centre,
     a, b, q, PA. Inside 60" the density changes by only ~15%, so the centre information comes from the gradient further
     out. HST coverage is complete to 250" (every 3" cell occupied); completeness of the faint slices beyond 60" untested.
Errors: bootstrap over stars (30 draws for the symmetry centre, 20 per model fit).
Reference centres (offsets from AvdM10, arcsec): Goldsbury+2010 (+0.41, +0.35); oMEGACat VI kinematic (-1.31+-0.72,
+1.44+-0.66); Haeberle+2024 IMBH position 0.77" NE (direction taken as 45 deg: approximate).
Output: plots/hst_density_centre.png, results/plot_data/hst_density_centre.json
Usage: PYTHONPATH=src python bin/hst_density_centre.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
SLICES = {"F625W<19": (-99, 19), "19-21": (19, 21), "21-23": (21, 23), "all<23": (-99, 23)}
REFS = {"Goldsbury+2010": ((0.41, 0.35), None), "oMEGACat VI kinematic": ((-1.31, 1.44), (0.72, 0.66)),
        "IMBH (Haeberle+2024, 0.77\" NE)": ((0.544, 0.544), None)}
RMAX = 205.


def symmetry(x, y, r_in=5., r_out=55., nw=16):
    """Coarse grid (0.5" over +-4") then fine grid (0.1" over +-0.6" around the coarse minimum)."""
    k = np.hypot(x, y) < r_out+5; x, y = x[k], y[k]
    c = _sym_grid(x, y, np.zeros(2), 4., 0.5, r_in, r_out, nw)
    return _sym_grid(x, y, c, 0.6, 0.1, r_in, r_out, nw)


def _sym_grid(x, y, c0, half, step, r_in, r_out, nw):
    g0 = np.arange(-half, half+step/2, step); gx, gy = c0[0]+g0, c0[1]+g0; chi = np.empty((len(g0), len(g0)))
    for i, cx in enumerate(gx):
        dx = x-cx
        for j, cy in enumerate(gy):
            dy = y-cy; R = np.hypot(dx, dy); k = (R >= r_in) & (R < r_out)
            w = np.floor((np.arctan2(dy[k], dx[k])+np.pi)/(2*np.pi)*nw).astype(int) % nw
            n = np.bincount(w, minlength=nw); a, b = n[:nw//2], n[nw//2:]
            chi[i, j] = np.sum((a-b)**2/np.maximum(a+b, 1))
    i, j = np.unravel_index(np.argmin(chi), chi.shape)
    def para(v, k, g):                                                 # parabolic refinement along one axis
        if 0 < k < len(v)-1:
            d = v[k-1]-2*v[k]+v[k+1]
            return g[k]+(0.5*(v[k-1]-v[k+1])/d*step if d > 0 else 0.)
        return g[k]
    return np.array([para(chi[:, j], i, gx), para(chi[i, :], j, gy)])


def model_fit(x, y, rmax, start=(0., 0.), p0=None):
    k = np.hypot(x, y) < rmax; x, y = x[k], y[k]; h = 0.5 if rmax <= 60 else 1.0
    gx = np.arange(-rmax+h/2, rmax, h); X, Y = np.meshgrid(gx, gx); inside = np.hypot(X, Y) < rmax
    Xi, Yi = X[inside], Y[inside]; dA = h**2

    def prof(px, py, p):
        cx, cy, la, b, lq, th = p
        c, s = np.cos(th), np.sin(th); u, v = (px-cx)*c+(py-cy)*s, -(px-cx)*s+(py-cy)*c
        return (1+(u*u+v*v/np.exp(2*lq))/np.exp(2*la))**(-b)

    def nll(p):
        if not (-0.6 < p[4] < 0.6 and 0 < p[3] < 5):
            return 1e30
        return -(np.sum(np.log(prof(x, y, p)))-len(x)*np.log(np.sum(prof(Xi, Yi, p))*dA))

    best = None
    for z in ([p0] if p0 is not None else [[start[0], start[1], np.log(100.), 1., 0., th0] for th0 in (0.3, 1.2, 2.1)]):
        z = np.asarray(z, float); simplex = np.vstack([z]+[z+d for d in np.diag([1.0, 1.0, 0.3, 0.2, 0.05, 0.3])])
        # explicit simplex: scipy's default step for a zero entry is 0.00025, which froze the centre at its start
        r = minimize(nll, z, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-3, fatol=1e-3, initial_simplex=simplex))
        if best is None or r.fun < best.fun:
            best = r
    return best.x


def main():
    from astropy.table import Table
    from ocen_dm.paths import raw_dir
    t = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    m = np.asarray(t["f625w"], float); xp, yp = np.asarray(t["x"], float), np.asarray(t["y"], float)
    ok = np.isfinite(m) & np.isfinite(xp) & np.isfinite(yp)
    x = -0.04*(xp[ok]-15000.); y = 0.04*(yp[ok]-15000.); m = m[ok]
    k = np.hypot(x, y) < RMAX+5; x, y, m = x[k], y[k], m[k]
    rng = np.random.default_rng(7); out = dict(slices={}, refs={k_: dict(offset=v[0], err=v[1]) for k_, v in REFS.items()})
    for lab, (a, b) in SLICES.items():
        s = (m >= a) & (m < b); xs, ys = x[s], y[s]; res = {}
        c = symmetry(xs, ys)
        bs = np.array([symmetry(xs[i], ys[i]) for i in (rng.integers(0, len(xs), len(xs)) for _ in range(30))])
        res["symmetry"] = dict(c=c.tolist(), err=bs.std(0).tolist())
        for rmax in (60., 120., 200.):
            p = model_fit(xs, ys, rmax)
            bs = np.array([model_fit(xs[i], ys[i], rmax, p0=p) for i in (rng.integers(0, len(xs), len(xs)) for _ in range(20))])
            res[f"model_{rmax:g}"] = dict(c=p[:2].tolist(), err=bs[:, :2].std(0).tolist(), a_arcsec=float(np.exp(p[2])), b=float(p[3]),
                                          q=float(np.exp(p[4])), pa_deg=float(np.degrees(p[5]) % 180), n=int((np.hypot(xs, ys) < rmax).sum()))
        out["slices"][lab] = dict(n_stars=int(s.sum()), n_within_200=int((np.hypot(xs, ys) < 200).sum()), methods=res)
        print(f"{lab:9s} N(<200\") {out['slices'][lab]['n_within_200']:7d}: " + " | ".join(
            f"{mk} ({v['c'][0]:+.2f}+-{v['err'][0]:.2f}, {v['c'][1]:+.2f}+-{v['err'][1]:.2f})" for mk, v in res.items()), flush=True)
        for r_ in (60, 120, 200):
            v = res[f"model_{r_}"]; print(f"          model_{r_}: a {v['a_arcsec']:.0f}\" b {v['b']:.2f} q {v['q']:.3f} PA {v['pa_deg']:.0f}")
    (ROOT/"results/plot_data/hst_density_centre.json").write_text(json.dumps(out, indent=1, default=float))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 6.5), constrained_layout=True)
    mk = {"symmetry": "^", "model_60": "o", "model_120": "s", "model_200": "*"}
    for zoom, x_ in zip((3.5, 1.2), ax):
        for i, (lab, r) in enumerate(out["slices"].items()):
            for mt, v in r["methods"].items():
                x_.errorbar(v["c"][0], v["c"][1], xerr=v["err"][0], yerr=v["err"][1], fmt=mk[mt], color=f"C{i}", ms=8 if mt == "model_200" else 5,
                            capsize=2, label=f"{lab}: {mt}" if zoom == 3.5 else None)
        x_.plot(0, 0, "k+", ms=16, mew=2, label="AvdM10 (adopted)" if zoom == 3.5 else None)
        for j, (lab, (o, e)) in enumerate(REFS.items()):
            x_.errorbar(o[0], o[1], xerr=None if e is None else e[0], yerr=None if e is None else e[1], fmt="x", color=f"C{6+j}", ms=10, mew=2,
                        capsize=3, label=lab if zoom == 3.5 else None)
        x_.set(xlim=(zoom, -zoom), ylim=(-zoom, zoom), aspect="equal", xlabel="dRA* [arcsec] (east left)", ylabel="dDec [arcsec]",
               title=f"density centres, +-{zoom}\"")
        x_.axhline(0, color="0.8", lw=0.5); x_.axvline(0, color="0.8", lw=0.5)
    ax[0].legend(fontsize=6, ncol=2, loc="lower left")
    fig.suptitle("omega Cen density centre from oMEGACat star positions, by magnitude (star bootstraps)", fontsize=10)
    fig.savefig(ROOT/"plots/hst_density_centre.png", dpi=110); print(ROOT/"plots/hst_density_centre.png")


if __name__ == "__main__":
    main()
