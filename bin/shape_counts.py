#!/usr/bin/env python3
"""Projected shape of omega Cen from star counts: data loaders and the ring table of azimuthal harmonics (2026-11-15).

Coordinates: arcsec on the sky, x = east, y = north, about the adopted centre (ocen_dm.cluster, = the oMEGACat
pixel centre 15000, 15000); position angle phi = atan2(x, y), east of north (the JAM/kinematic convention, PA 104.2).
Gaia positions are gnomonic (tangent-plane) projections about the centre.

Samples (docs/JAM_FIT_PLAN.md step 4):
  hst   oMEGACat VI, all stars with measured F625W < mag (default 19; 18 = bright cut), 60-250" (full coverage).
  gaia  Vasiliev & Baumgardt 2021 field catalogue, G < mag (default 17; 16 = bright cut), 300-2400";
        members=True: position-independent PM + parallax cut, chi2 = |pm - MU_SYS|^2 / (s_int^2 + err^2) < 9 per
        component pair (s_int = 0.5 mas/yr: dispersion + rotation), |parallax - 1/D| < 3 parallax_error.
  kuzma Kuzma et al. 2025 periphery members (Pristine + Gaia), beyond 2400": PA / tail direction only.

Ring table: complex harmonics Z_m = <exp(i m phi)> per ring, amplitude |Z_m| with the noise-bias-corrected value
sqrt(max(|Z|^2 - 1/N, 0)), phase PA_m = arg(Z_m)/m. For m = 2 the major-axis PA is PA_2 (density peaks at phi = PA,
PA + 180). Small-eps conversion for a power-law ring (Sigma ~ R^gamma, gamma < 0): eps = -4 |Z_2| / gamma, diluted by
the field fraction f (eps = -4 |Z_2| / (gamma (1-f))). gamma from the counts; f from a uniform field (Gaia, all stars).
The ring table is a diagnostic; the measurement is the unbinned fit (bin/fit_shape_counts.py).

Usage: PYTHONPATH=src python bin/shape_counts.py      -> results/shape/ring_table_20261115.json, plots/shape_rings_20261115.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import numpy as np
from astropy.io import fits
from astropy.table import Table

from ocen_dm.cluster import MU_SYS, OCEN_DEC, OCEN_DISTANCE_KPC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

PA_KIN = 104.2
S_INT = 0.5          # mas/yr, intrinsic PM width allowed in the member cut


def gnomonic(ra, dec, ra0=OCEN_RA, dec0=OCEN_DEC):
    """Tangent-plane (x east, y north) in arcsec."""
    a, d, a0, d0 = (np.radians(np.asarray(v, float)) for v in (ra, dec, ra0, dec0))
    c = np.sin(d0)*np.sin(d)+np.cos(d0)*np.cos(d)*np.cos(a-a0)
    x = np.cos(d)*np.sin(a-a0)/c
    y = (np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(a-a0))/c
    return np.degrees(x)*3600., np.degrees(y)*3600.


def load_hst(mag=19., rmin=60., rmax=250.):
    t = fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.); y = 0.04*(np.asarray(t["y"], float)-15000.)
    f = np.asarray(t["f625w"], float)
    R = np.hypot(x, y)
    k = np.isfinite(x) & np.isfinite(y) & np.isfinite(f) & (f < mag) & (R >= rmin) & (R < rmax)
    return dict(x=x[k], y=y[k], rmin=rmin, rmax=rmax, name=f"HST F625W<{mag:g}")


def member_cut(t, s_int=S_INT):
    dx = np.asarray(t["pmra"], float)-MU_SYS[0]; dy = np.asarray(t["pmdec"], float)-MU_SYS[1]
    ex, ey = np.asarray(t["pmra_error"], float), np.asarray(t["pmdec_error"], float)
    chi2 = dx**2/(s_int**2+ex**2)+dy**2/(s_int**2+ey**2)
    plx, eplx = np.asarray(t["parallax"], float), np.asarray(t["parallax_error"], float)
    return (chi2 < 9.) & (np.abs(plx-1/OCEN_DISTANCE_KPC) < 3*eplx)


def load_gaia(mag=17., rmin=300., rmax=2400., members=False, s_int=S_INT):
    t = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    x, y = gnomonic(t["ra"], t["dec"])
    g = np.asarray(t["g_mag"], float)
    R = np.hypot(x, y)
    k = np.isfinite(g) & (g < mag) & (R >= rmin) & (R < rmax)
    if members:
        k &= member_cut(t, s_int)
    return dict(x=x[k], y=y[k], rmin=rmin, rmax=rmax,
                name=f"Gaia G<{mag:g}"+(f" PM members (s_int {s_int:g})" if members else " all"))


def load_kuzma(pmin=0.5, rmin=2400.):
    t = Table.read(processed_dir()/"tails"/"kuzma2025_periphery.ecsv")
    pcol = "membership_prob"
    x, y = gnomonic(t["ra"], t["dec"])
    R = np.hypot(x, y)
    k = (np.asarray(t[pcol], float) > pmin) & (R >= rmin)
    return dict(x=x[k], y=y[k], rmin=rmin, rmax=float(R[k].max()), name=f"Kuzma25 {pcol}>{pmin:g}")


def harmonics(x, y, m):
    phi = np.arctan2(x, y)
    z = np.mean(np.exp(1j*m*phi)); n = len(x)
    amp = np.sqrt(max(abs(z)**2-1/n, 0.)); err = np.sqrt(1/(2*n))*np.sqrt(2)   # s.d. of |Z| for a uniform ring ~ 1/sqrt(N)
    pa = np.degrees(np.angle(z))/m % (360./m)
    pa_err = np.degrees(err/max(abs(z), 1e-9))/m
    return dict(amp=float(amp), raw=float(abs(z)), err=float(1/np.sqrt(n)), pa=float(pa), pa_err=float(min(pa_err, 90./m)))


def ring_table(s, edges, field=0.0):
    """field: uniform field density (stars/arcsec^2) for the dilution correction and gamma."""
    R = np.hypot(s["x"], s["y"]); out = []
    n_all = np.array([np.sum((R >= lo) & (R < hi)) for lo, hi in zip(edges[:-1], edges[1:])])
    area = np.pi*(edges[1:]**2-edges[:-1]**2)
    clus = np.maximum(n_all/area-field, 1e-12)
    lr = np.log(np.sqrt(edges[1:]*edges[:-1]))
    gam = np.gradient(np.log(clus), lr)
    for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        sel = (R >= lo) & (R < hi)
        h = {m: harmonics(s["x"][sel], s["y"][sel], m) for m in (1, 2, 3, 4)}
        f = min(field*area[k]/max(n_all[k], 1), 0.99)
        eps = -4*h[2]["amp"]/gam[k]/(1-f) if gam[k] < 0 else np.nan
        eps_err = 4*h[2]["err"]/abs(gam[k])/(1-f)
        out.append(dict(r_lower=float(lo), r_upper=float(hi), n=int(n_all[k]), field_frac=float(f), gamma=float(gam[k]),
                        eps=float(eps), eps_err=float(eps_err), **{f"m{m}": h[m] for m in h}))
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    F = json.loads((processed_dir()/"kinematics"/"tracer_mge_counts_c.json").read_text())["fit"]["gaia_field_per_arcmin2"]/3600.
    he = np.array([60, 75, 90, 105, 120, 140, 160, 180, 200, 225, 250.])
    ge = np.array([300, 350, 400, 450, 500, 600, 700, 850, 1000, 1150, 1300, 1500, 1700, 1950, 2150, 2400.])
    samples = {}
    for mag in (19., 18.):
        s = load_hst(mag); samples[s["name"]] = (s, ring_table(s, he))
    for mag in (17., 16.):
        s, m = load_gaia(mag), load_gaia(mag, members=True)
        # uniform field of the all-star sample = density of NON-members at 1500-2400" (cluster share there is small
        # after the PM cut); members: residual field neglected in the table (dilution not corrected)
        Rs, Rm = np.hypot(s["x"], s["y"]), np.hypot(m["x"], m["y"])
        a = np.pi*(2400**2-1500**2)
        fld = (np.sum(Rs >= 1500)-np.sum(Rm >= 1500))/a
        print(f"G<{mag:g}: non-member density at 1500-2400\" {fld*3600:.3f}/arcmin^2 (count-MGE fit, G<17: {F*3600:.3f})")
        samples[s["name"]] = (s, ring_table(s, ge, field=fld))
        samples[m["name"]] = (m, ring_table(m, ge, field=0.0))
    try:
        s = load_kuzma(); samples[s["name"]] = (s, ring_table(s, np.geomspace(2400, max(s["rmax"], 2500), 7)))
    except Exception as e:                       # noqa: BLE001
        print("kuzma skipped:", e)
    out = ROOT/"results/shape"; out.mkdir(parents=True, exist_ok=True)
    (out/"ring_table_20261115.json").write_text(json.dumps({k: v[1] for k, v in samples.items()}, indent=1))
    fig, ax = plt.subplots(4, 1, figsize=(9, 13), sharex=True)
    for j, (name, (s, tab)) in enumerate(samples.items()):
        r = np.array([np.sqrt(t["r_lower"]*t["r_upper"]) for t in tab])*(1+0.015*j)
        e = np.array([t["eps"] for t in tab]); ee = np.array([t["eps_err"] for t in tab])
        ax[0].errorbar(r, e, ee, fmt="o-", ms=3, lw=.8, label=name)
        pa = np.array([t["m2"]["pa"] for t in tab]); pe = np.array([t["m2"]["pa_err"] for t in tab])
        ax[1].errorbar(r, pa, pe, fmt="o", ms=3, lw=.8)
        for i, m in enumerate((1, 3)):
            ax[2+i].errorbar(r, [t[f"m{m}"]["amp"] for t in tab], [t[f"m{m}"]["err"] for t in tab], fmt="o-", ms=3, lw=.8)
    vdv_s = np.array([0.153, 1.477, 2.525, 3.691, 5.219, 7.534, 11.069, 17.547])*60
    vdv_q = np.array([1.000, 0.934, 0.877, 0.848, 0.850, 0.836, 0.866, 0.926])
    ax[0].plot(vdv_s, 1-vdv_q, "k*", ms=9, label="van de Ven+06 MGE q'_j at R = sigma_j (Geyer+83)")
    ax[0].set_ylim(-0.1, 0.4); ax[0].set_ylabel("eps (small-eps, from |Z2|)"); ax[0].legend(fontsize=7); ax[0].axhline(0, c="k", lw=.5)
    ax[1].axhline(PA_KIN, c="k", ls="--", lw=.8); ax[1].set_ylabel("PA_2 (deg E of N)"); ax[1].set_ylim(0, 180)
    ax[2].set_ylabel("|Z1| (lopsided)"); ax[3].set_ylabel("|Z3|"); ax[3].set_xscale("log"); ax[3].set_xlabel("R (arcsec)")
    for a in ax: a.grid(alpha=.3)
    plt.tight_layout(); p = ROOT/"plots/shape_rings_20261115.png"; plt.savefig(p, dpi=100); print(p)
    for name, (s, tab) in samples.items():
        print(f"\n{name} ({len(s['x'])} stars)")
        for t in tab:
            print(f"  {t['r_lower']:6.0f}-{t['r_upper']:6.0f} N {t['n']:6d} f {t['field_frac']:.2f} gam {t['gamma']:5.2f} "
                  f"eps {t['eps']:6.3f}+-{t['eps_err']:.3f} PA {t['m2']['pa']:5.1f}+-{t['m2']['pa_err']:4.1f} "
                  f"|Z1| {t['m1']['amp']:.3f} PA1 {t['m1']['pa']:5.0f} |Z3| {t['m3']['amp']:.3f}")


if __name__ == "__main__":
    main()
