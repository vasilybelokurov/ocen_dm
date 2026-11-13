#!/usr/bin/env python3
"""Data vs model for the joint rotation-field fit (bin/fit_rotation_field.py).

Panels: (a) LOS rotation amplitude per annulus: MUSE per-star velocities binned in rings (ML fit of v_sys + A cos(phi')
with phi' measured from the fitted PA, intrinsic sigma from the published profile, per-star errors) vs the model's
amplitude in the same rings (posterior draws); (b) ring-averaged tangential PM rotation (oMEGACat II Table 10 and our
Gaia EDR3) vs the model; (c) intrinsic v_phi(R) with 68% band; panel (b) also shows the VB21 PM curve (median, 16-84%); inset text: i, PA, v_sys.
Usage: PYTHONPATH=src:bin python bin/plot_rotation_field.py [--law cylindrical|spherical]
Output: plots/rotation_field_<law>_20261113.png
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_rotation_field as F

EDGES = np.array([0, 10, 20, 35, 55, 80, 110, 150, 200, 260, 345.])


def ring_amplitude(d, pa, vsys, sel):
    """ML amplitude A of v - vsys = A cos(phi') in one ring (variance fixed): linear least squares."""
    a = np.radians(pa)
    xp = d["x"][sel]*np.sin(a)+d["y"][sel]*np.cos(a)
    c = xp/np.hypot(d["x"][sel], d["y"][sel])
    w = 1/d["var"][sel]
    A = np.sum(w*c*(d["v"][sel]-vsys))/np.sum(w*c*c)
    return A, 1/np.sqrt(np.sum(w*c*c))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--law", default="cylindrical")
    law = ap.parse_args().law
    res = json.loads((ROOT/f"results/rotation/rotation_field_{law}_20261113.json").read_text())
    chain = np.load(ROOT/f"results/rotation/rotation_field_{law}_chain_20261113.npy")
    d = F.load_data(); pb = F.Problem(d, F.Field(law))
    pmap = np.array(res["map"]); nn = len(F.NODES)
    R = np.hypot(d["x"], d["y"])
    rmid, Ad, Ae = [], [], []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        s = (R >= lo) & (R < hi)
        A, e = ring_amplitude(d, pmap[-2], pmap[-1], s); rmid.append(np.median(R[s])); Ad.append(A); Ae.append(e)
    rmid = np.array(rmid)
    draws = chain[np.random.default_rng(0).choice(len(chain), 150, replace=False)]
    Am, ringm, vint = [], [], []
    rr = np.geomspace(5, 2500, 200)
    for p in draws:
        mlos, vt, _ = pb.f.project(p[:nn], p[-3])
        A = 2*np.mean(mlos*np.cos(F.PHI)[None, :], axis=1)
        Am.append(np.interp(np.log(rmid), np.log(F.RP), A))
        ringm.append(np.abs(np.interp(np.log(rr), np.log(F.RP), np.mean(vt, axis=1))))
        vint.append(F.vcurve(p[:nn], rr))
    qA, qr, qv = (np.percentile(np.array(a), [16, 50, 84], axis=0) for a in (Am, ringm, vint))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.3))
    ax[0].errorbar(rmid, Ad, Ae, fmt="o", color="C3", ms=4, label="MUSE stars, ring fit (fixed PA)")
    ax[0].plot(rmid, qA[1], color="k"); ax[0].fill_between(rmid, qA[0], qA[2], color="k", alpha=.2, lw=0, label="model")
    ax[0].set_xscale("log"); ax[0].set_xlabel("R [arcsec]"); ax[0].set_ylabel("LOS rotation amplitude [km/s]")
    ax[0].axhline(0, color="k", lw=.5); ax[0].legend(fontsize=7.5, frameon=False)
    src = np.array(res["pm_r"]); pv, pe = np.array(res["pm_v"]), np.array(res["pm_e"])
    t10 = np.arange(len(src)) < 12
    ax[1].errorbar(src[t10], pv[t10], pe[t10], fmt="s", color="k", ms=4, label="oMEGACat II Table 10 (Gaia FPR − HST)")
    ax[1].errorbar(src[~t10], pv[~t10], pe[~t10], fmt="o", color="C2", ms=4, label="our Gaia EDR3 mixture fit")
    ax[1].plot(rr, qr[1], color="k"); ax[1].fill_between(rr, qr[0], qr[2], color="k", alpha=.2, lw=0, label="model")
    from astropy.table import Table
    from ocen_dm.paths import processed_dir
    vb = Table.read(processed_dir()/"kinematics"/"vasiliev2021_ocen_pm_profiles.ecsv")
    k = np.asarray(vb["r"]) > 0
    ax[1].fill_between(np.asarray(vb["r"])[k], np.asarray(vb["vrot_pm_p16"])[k]*F.KMS, np.asarray(vb["vrot_pm_p84"])[k]*F.KMS,
                       color="C0", alpha=.25, lw=0)
    ax[1].plot(np.asarray(vb["r"])[k], np.asarray(vb["vrot_pm"])[k]*F.KMS, color="C0", lw=1.5, ls="--",
               label="Vasiliev & Baumgardt 2021 (current HST/Gaia term)")
    ax[1].set_xscale("log"); ax[1].set_xlabel("R [arcsec]"); ax[1].set_ylabel("ring-mean tangential PM rotation [km/s]")
    ax[1].legend(fontsize=7.5, frameon=False)
    ax[2].plot(rr, qv[1], color="C0"); ax[2].fill_between(rr, qv[0], qv[2], color="C0", alpha=.25, lw=0)
    ax[2].set_xscale("log"); ax[2].set_xlabel(f"{'cylindrical R' if law == 'cylindrical' else 'spherical r'} [arcsec]")
    ax[2].set_ylabel(r"intrinsic $v_\phi$ [km/s]")
    p16, p50, p84 = (np.array(res[k]) for k in ("p16", "p50", "p84"))
    txt = "\n".join(f"{n} = {b:.2f} (+{c-b:.2f}/−{b-a:.2f})" for n, a, b, c in zip(("i [deg]", "PA [deg]", "v_sys [km/s]"),
                                                                                p16[-3:], p50[-3:], p84[-3:]))
    ax[2].text(0.03, 0.97, f"{law} rotation, q'={res['q_proj']}\n"+txt+f"\nχ²(PM rings) = {res['chi2_pm']:.1f}/{res['n_pm']}",
               transform=ax[2].transAxes, va="top", fontsize=8)
    for a in ax:
        a.set_xlim(4, 2600)
    fig.tight_layout()
    out = ROOT/f"plots/rotation_field_{law}_20261113.png"
    fig.savefig(out, dpi=130)
    print(out)
    print("LOS ring amplitude data:", np.round(Ad, 2).tolist(), "model:", np.round(qA[1], 2).tolist())


if __name__ == "__main__":
    main()
