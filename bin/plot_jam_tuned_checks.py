#!/usr/bin/env python3
"""Diagnostic plots behind the adversarial review of the tuned flattened JAM models (2026-11-13).

Panels: (a) HST tangential rotation term streaming2, old (snapshot) vs corrected (measured field), with sigma_T^2 data;
(b) anisotropy: JAM's meridional beta vs spherical-equivalent beta_sph = 1-(s_th^2+<vphi^2>)/(2 s_r^2), shell-averaged
with tracer weight, from JamPy intrinsic moments; (c) chi2 vs small shifts of beta_inf around the q'0.95 and q'0.90
optima (stability); (d) sigma_T/sigma_R vs radius, data and tuned models; (e) enclosed mass M(<r) of each tuned
model's (round-parent) mass model; (f) M_star, M_rem and total per shape.
Usage: PYTHONPATH=src:bin python bin/plot_jam_tuned_checks.py  -> plots/jam_tuned_checks_20261113.png
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from jampy.axi.jam_axi_intr import intrinsic_moments
import jam_mock_rotation_test as J
from tune_jam_to_data import Setup

TAGS = [("q1.00_r1.00_i90", "round", "k"), ("q0.95_r0.95_i44", "q'=0.95", "C0"),
        ("q0.90_r0.90_i44", "q'=0.90", "C3"), ("q0.85_r0.85_i44", "q'=0.85", "C2")]
S = Setup()
res = {t: json.loads((ROOT/f"results/rotation/jam_tuned_{t}.json").read_text()) for t, *_ in TAGS}
theta = {t: dict(zip(r["names"], r["x"])) for t, r in res.items()}
new = {p.name: p for p in S.data.profiles}; old = {p.name: p for p in S.pb0.data.profiles}
fig, ax = plt.subplots(2, 3, figsize=(15, 9)); ax = ax.ravel()

# (a) rotation term
p, p0 = new["hst_pm_tangential_ours"], old["hst_pm_tangential_ours"]
a = ax[0]
a.plot(p.r, p.value**2, "ko", ms=3, label=r"data $\sigma_T^2$")
a.plot(p.r, p0.streaming2, "s-", color="C1", ms=3, label="rotation term, old (VB21 curve)")
a.plot(p.r, p.streaming2, "^-", color="C4", ms=3, label="rotation term, corrected (measured field)")
a.set_xscale("log"); a.set_xlabel("R [arcsec]"); a.set_ylabel(r"[(mas/yr)$^2$]")
a.set_title("(a) HST tangential: rotation term added back"); a.legend(fontsize=8)

# (b) anisotropy
rr = np.geomspace(0.3, 60, 30); th = np.linspace(0.02, np.pi/2-0.02, 40)
a = ax[1]
for t, lab, col in TAGS:
    d = res[t]; jeans, D, _ = S.fam.build(theta[t])
    lum, pot, mbh = J.jam_models(jeans, D, d["qobs"], d["qrem"])
    inc = np.radians(d["inc"]); pc = D*1e-3*np.pi/0.648
    def deproj(m):
        s, sg, q = map(np.asarray, m); qi = np.sqrt(q**2-np.cos(inc)**2)/np.sin(inc); sgp = sg*pc
        return (s*q/(sgp*qi*np.sqrt(2*np.pi)), sgp, qi)
    L, P = deproj(lum), deproj(pot)
    R = (rr[:, None]*np.sin(th)[None]).ravel(); z = (rr[:, None]*np.cos(th)[None]).ravel()
    s2r, s2t, _, v2p, nu = [x.reshape(len(rr), len(th)) for x in intrinsic_moments(
        R, z, L, P, mbh, lambda r, t_: (jeans.anisotropy.beta(r), 0.), 0., None, 31)]
    w = nu*np.sin(th)[None]; av = lambda x: np.sum(w*x, 1)/np.sum(w, 1)
    a.plot(rr, 1-(av(s2t)+av(v2p))/(2*av(s2r)), "-", color=col, label=f"{lab}: spherical-equivalent")
    a.plot(rr, jeans.anisotropy.beta(rr), "--", color=col, lw=0.9)
a.plot([], [], "k--", lw=0.9, label="dashed: JAM's own beta (2 of 3 directions)")
a.axhline(0, color="0.6", lw=0.5); a.set_xscale("log"); a.set_ylim(-1.4, 0.8)
a.set_xlabel("r [pc]"); a.set_ylabel(r"$\beta$"); a.set_title("(b) orbit anisotropy: two definitions"); a.legend(fontsize=7)

# (c) stability
a = ax[2]
db = np.linspace(-2e-3, 2e-3, 17)
for t, lab, col in TAGS[1:3]:
    c2 = []
    for dd in db:
        x = dict(theta[t]); x["beta_inf"] += dd
        try:
            pred, _ = S.predict(x, res[t]["qobs"], res[t]["qrem"], res[t]["inc"])
            c2.append(sum(v["chi2"] for v in S.stats(pred).values()))
        except ValueError:
            c2.append(np.nan)
    c2 = np.array(c2)
    a.plot(db*1e3, c2, "o-", color=col, ms=3, label=lab)
    for dd in db[~np.isfinite(c2)]:
        a.axvline(dd*1e3, color=col, alpha=0.25, lw=4)
a.set_xlabel(r"shift of $\beta_\infty$ from best value [$10^{-3}$]"); a.set_ylabel(r"total $\chi^2$ (89 bins)")
a.set_title("(c) stability (shaded: model breaks, negative moments)"); a.legend(fontsize=8)

# (d) sigma_T / sigma_R
a = ax[3]
for inst, mk in (("hst_pm_{}_ours", "o"), ("gaia_edr3_ours_{}", "s")):
    pr, pt = new[inst.format("radial")], new[inst.format("tangential")]
    assert np.allclose(pr.r, pt.r)
    rat = pt.value/pr.value
    e = rat*np.hypot(0.5*(pt.err_lo+pt.err_hi)/pt.value, 0.5*(pr.err_lo+pr.err_hi)/pr.value)
    a.errorbar(pr.r, rat, e, fmt=mk, color="0.4", ms=3, lw=0.7, mfc="none" if mk == "s" else "0.4")
    for t, lab, col in TAGS:
        a.plot(pr.r, np.array(res[t]["pred"][inst.format("tangential")])/np.array(res[t]["pred"][inst.format("radial")]),
               "-", color=col, lw=1.1, label=lab if mk == "o" else None)
a.axhline(1, color="0.6", lw=0.5); a.set_xscale("log"); a.set_xlabel("R [arcsec]")
a.set_ylabel(r"$\sigma_T/\sigma_R$ (HST ●, Gaia □)"); a.set_title(r"(d) $\sigma_T/\sigma_R$: data vs tuned models"); a.legend(fontsize=8)

# (e) M(<r), (f) masses
r = np.geomspace(0.1, 100, 200)
for t, lab, col in TAGS:
    jeans, D, _ = S.fam.build(theta[t])
    ax[4].plot(r, jeans.mass.enclosed_mass(r), color=col, label=lab)
    x = theta[t]
    for k, (v, hatch) in enumerate(((x["M_star"], ""), (x["M_rem"], "//"), (x["M_star"]+x["M_rem"], ".."))):
        i = TAGS.index((t, lab, col))
        ax[5].bar(k*5+i, v/1e6, color=col, hatch=hatch, alpha=0.8)
ax[4].set_xscale("log"); ax[4].set_yscale("log"); ax[4].set_xlabel("r [pc]"); ax[4].set_ylabel(r"M(<r) [M$_\odot$]")
ax[4].set_title("(e) enclosed mass of each tuned model"); ax[4].legend(fontsize=8)
ax[5].set_xticks([1.5, 6.5, 11.5]); ax[5].set_xticklabels(["stars", "remnants", "stars + remnants"])
ax[5].set_ylabel(r"[10$^6$ M$_\odot$]"); ax[5].set_title("(f) masses per shape (bars: round, .95, .90, .85)")
fig.tight_layout()
out = ROOT/"plots/jam_tuned_checks_20261113.png"; fig.savefig(out, dpi=120); print(out)
