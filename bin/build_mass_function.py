#!/usr/bin/env python3
"""Present-day stellar mass function of omega Cen from oMEGACat photometry, and the
implied luminous mass per counted bright star (an external prior on M_star).

Steps
1. Sample: all catalogue stars with measured F625W and F814W (no quality flags,
   which are magnitude- and radius-dependent), in radial zones about the pixel
   centre, restricted per zone to magnitudes where the luminosity function is
   still rising smoothly (completeness limits set by hand from the LF turnover).
2. Isochrone: MIST v1.2 (minimint, prefix ~/data/isochrones/minimint_hst) in
   ACS_WFC_F625W and ACS_WFC_F814W, (m-M)_0 from the adopted distance, A_V from
   E(B-V) = 0.12 (Harris) applied through the bolometric-correction A_V axis.
   The effective [Fe/H] is chosen by matching the main-sequence ridge line
   (median colour per magnitude); the metallicity spread is a systematic.
3. Mass-luminosity relation on the main sequence -> dN/dm per zone from the LF
   of stars within the ridge-line colour window; power-law slope fitted by
   Poisson likelihood in mass bins.
4. Mass budget per zone: main-sequence stars (observed range + extrapolation to
   0.1 Msun with stated slopes), evolved stars counted above the turnoff (mass
   ~ m_TO), white dwarfs from the IMF above the turnoff (continuous with the
   observed slope, or Kroupa -2.3) with mean mass 0.55 Msun. NS/BH progenitors
   are reported, not added (they belong to the free remnant component).
5. Output: luminous mass that follows the light per star with F625W < 19 in
   each zone (the count product's selection), with variants, compared with
   the value implied by the dynamical fit (1/a in Msun per star).

Everything is written to results/mass_function/ (JSON) and plots/.

Usage: python bin/build_mass_function.py [--age 12.5] [--distance 5.43] [--ebv 0.12]
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar

ARCSEC_PER_RAD = 206264.806
BANDS = ["ACS_WFC_F625W", "ACS_WFC_F814W"]
# zone edges [arcsec] and per-zone faint limit in F625W for the LF display. The mass-function
# slope is fitted to a COMMON limit in all zones: with zone-dependent limits the outer zones
# reached magnitudes where completeness already falls and the slope flattened spuriously
# (JOURNAL 2026-09-23). F625W < 21.0 is ~0.52 Msun and well inside the complete range everywhere.
ZONES = [(10., 30., 21.5), (30., 60., 21.5), (60., 100., 22.0), (100., 175., 22.0), (175., 250., 22.5)]
MF_MAGLIM = 21.0
M_MIN_EXTRAP = 0.10          # hydrogen-burning limit for the extrapolation
WD_MASS = 0.55               # mean white-dwarf mass for an old population (IFMR)
FEH_GRID = np.array([-2.0, -1.8, -1.6, -1.53, -1.4, -1.2, -1.0])


def load_catalogue():
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits")[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.)
    y = 0.04*(np.asarray(t["y"], float)-15000.)
    f6 = np.asarray(t["f625w"], float); f8 = np.asarray(t["f814w"], float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(f6) & np.isfinite(f8)
    return np.hypot(x, y)[ok], f6[ok], f8[ok]


class Isochrone:
    """MIST main sequence in the two ACS bands with distance and extinction applied."""

    def __init__(self, age_gyr, feh, distance_kpc, ebv, prefix=None):
        import minimint
        prefix = prefix or os.path.expanduser("~/data/isochrones/minimint_hst")+"/"
        self.ii = minimint.Interpolator(BANDS, data_prefix=prefix)
        self.age, self.feh, self.dm0, self.av = age_gyr, feh, 5*np.log10(distance_kpc*1e3/10.), 3.1*ebv
        mass = np.geomspace(.09, 1.2, 1200)
        res = self.ii(mass, np.log10(age_gyr*1e9), feh)
        good = np.isfinite(res[BANDS[0]]) & (res["phase"] <= 0)           # main sequence (phase 0), alive
        self.mass = mass[good]
        # extinction from the BC tables: A_band = BC(Av=0) - BC(Av)
        p0 = np.column_stack((res["logteff"][good], res["logg"][good], np.full(good.sum(), feh), np.zeros(good.sum())))
        p1 = p0.copy(); p1[:, 3] = self.av
        bc0, bc1 = self.ii.bolomInt(p0), self.ii.bolomInt(p1)
        self.m625 = res[BANDS[0]][good]+self.dm0+(bc0[BANDS[0]]-bc1[BANDS[0]])
        self.m814 = res[BANDS[1]][good]+self.dm0+(bc0[BANDS[1]]-bc1[BANDS[1]])
        self.a625 = float(np.median(bc0[BANDS[0]]-bc1[BANDS[0]])); self.a814 = float(np.median(bc0[BANDS[1]]-bc1[BANDS[1]]))
        self.m_to = float(self.mass.max())
        order = np.argsort(self.m625)                                        # F625W decreases with mass on the MS
        self._m_of_mag = (self.m625[order], self.mass[order])

    def colour_at(self, mag):
        return np.interp(mag, self._m_of_mag[0], (self.m625-self.m814)[np.argsort(self.m625)])

    def mass_at(self, mag):
        return np.interp(mag, *self._m_of_mag)

    def mag_at(self, mass):
        order = np.argsort(self.mass)
        return np.interp(mass, self.mass[order], self.m625[order])


def ridge(f6, col, lo=19., hi=21.5, step=.25):
    """Observed main-sequence ridge: median colour and 16-84% half-width per magnitude bin."""
    edges = np.arange(lo, hi+1e-9, step)
    mid, med, hw = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        k = (f6 >= a) & (f6 < b)
        if k.sum() > 200:
            q = np.percentile(col[k], [16, 50, 84])
            mid.append(.5*(a+b)); med.append(q[1]); hw.append(.5*(q[2]-q[0]))
    return np.array(mid), np.array(med), np.array(hw)


def powerlaw_fit(mass_lo, mass_hi, counts):
    """Poisson fit of dN/dm = k m^alpha to counts in mass bins (alpha profiled by 1D search)."""
    lo, hi, n = map(np.asarray, (mass_lo, mass_hi, counts))

    def negloglike(alpha):
        integ = (hi**(alpha+1)-lo**(alpha+1))/(alpha+1) if abs(alpha+1) > 1e-9 else np.log(hi/lo)
        k = n.sum()/integ.sum()
        mu = k*integ
        return -np.sum(n*np.log(mu)-mu)
    r = minimize_scalar(negloglike, bounds=(-5., 3.), method="bounded")
    alpha = r.x
    # curvature -> 1 sigma
    h = 1e-3
    curv = (negloglike(alpha+h)-2*negloglike(alpha)+negloglike(alpha-h))/h**2
    integ = (hi**(alpha+1)-lo**(alpha+1))/(alpha+1)
    return float(alpha), float(1/np.sqrt(max(curv, 1e-12))), float(n.sum()/integ.sum())


def mass_integral(k, alpha, a, b):
    """int_a^b k m^{alpha} m dm."""
    return k*(b**(alpha+2)-a**(alpha+2))/(alpha+2)


def number_integral(k, alpha, a, b):
    return k*(b**(alpha+1)-a**(alpha+1))/(alpha+1) if abs(alpha+1) > 1e-9 else k*np.log(b/a)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--age", type=float, default=12.5, help="Gyr")
    parser.add_argument("--distance", type=float, default=5.43, help="kpc")
    parser.add_argument("--ebv", type=float, default=.12)
    parser.add_argument("--tag", default="20260923")
    args = parser.parse_args()
    R, f6, f8 = load_catalogue()
    col = f6-f8
    inzone = (R >= ZONES[0][0]) & (R < ZONES[-1][1])

    # ---- 1. metallicity from the ridge line (all zones together, 19-21.5)
    mid, med, hw = ridge(f6[inzone], col[inzone])
    fits = []
    for feh in FEH_GRID:
        iso = Isochrone(args.age, feh, args.distance, args.ebv)
        resid = med-iso.colour_at(mid)
        fits.append((float(np.sqrt(np.mean(resid**2))), float(feh), float(np.mean(resid))))
    fits.sort()
    rms_best, feh_best, _ = fits[0]
    iso = Isochrone(args.age, feh_best, args.distance, args.ebv)
    print(f"ridge-line fit at {args.age:g} Gyr, (m-M)0 {iso.dm0:.3f}, A_V {iso.av:.3f} (A_625 {iso.a625:.3f}, A_814 {iso.a814:.3f}):")
    for rms, feh, mean in fits:
        print(f"   [Fe/H] {feh:+.2f}: colour rms {rms:.4f} mag, mean offset {mean:+.4f}")
    print(f"   adopted [Fe/H] = {feh_best:+.2f}; turnoff mass {iso.m_to:.3f} Msun; F625W(0.5 Msun) = {iso.mag_at(.5):.2f}")

    # ---- 2. per-zone luminosity and mass functions
    zones_out = []
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    ax_cmd, ax_lf, ax_mf, ax_ratio, ax_alpha, ax_budget = axes.ravel()
    for zi, (lo, hi, maglim) in enumerate(ZONES):
        kz = (R >= lo) & (R < hi)
        # ridge-line window: +-3 half-widths (interpolated in magnitude) around the isochrone colour
        hw_at = np.interp(f6, mid, hw, left=hw[0], right=hw[-1])
        on_ms = kz & (np.abs(col-iso.colour_at(f6)) < 3*hw_at) & (f6 > iso.mag_at(iso.m_to)-.05)
        mag_to = iso.mag_at(iso.m_to)
        # mass bins from the faint limit up to the turnoff
        m_faint = float(iso.mass_at(MF_MAGLIM))
        mbins = np.linspace(m_faint, iso.m_to, 8)
        counts = []
        for a, b in zip(mbins[:-1], mbins[1:]):
            mag_hi, mag_lo = iso.mag_at(a), iso.mag_at(b)        # fainter mag for lower mass
            counts.append(int(np.sum(on_ms & (f6 >= mag_lo) & (f6 < mag_hi))))
        counts = np.array(counts)
        alpha, alpha_err, k = powerlaw_fit(mbins[:-1], mbins[1:], counts)
        # observed MS mass in the fitted range
        m_ms_obs = mass_integral(k, alpha, m_faint, iso.m_to)
        # extrapolations below m_faint to 0.1 Msun
        variants = {"fitted slope": alpha, "Kroupa (-1.3 below 0.5)": None, "flat (alpha=0)": 0.}
        m_extrap = {}
        for name, al in variants.items():
            if al is None:   # Kroupa: continue fitted slope down to 0.5, then -1.3 (or the fitted slope if m_faint < 0.5)
                if m_faint > .5:
                    part = mass_integral(k, alpha, .5, m_faint)
                    k2 = k*.5**(alpha+1.3)
                    part += mass_integral(k2, -1.3, M_MIN_EXTRAP, .5)
                else:
                    k2 = k*m_faint**(alpha+1.3)
                    part = mass_integral(k2, -1.3, M_MIN_EXTRAP, m_faint)
                m_extrap[name] = part
            else:
                k2 = k*m_faint**(alpha-al)
                m_extrap[name] = mass_integral(k2, al, M_MIN_EXTRAP, m_faint)
        # evolved stars: brighter than the turnoff, redder than the MS extension (excludes HB/blue stragglers roughly)
        evolved = kz & (f6 < mag_to) & (col > .45)
        n_evolved = int(evolved.sum())
        # white dwarfs from the IMF above the turnoff, continuous at m_TO
        wd = {}
        for name, al_imf in (("IMF continuous with fitted slope", alpha), ("Kroupa/Salpeter IMF (-2.3)", -2.3)):
            k_imf = k*iso.m_to**(alpha-al_imf)
            wd[name] = dict(n_wd=number_integral(k_imf, al_imf, iso.m_to, 8.), n_ns=number_integral(k_imf, al_imf, 8., 25.),
                            n_bh=number_integral(k_imf, al_imf, 25., 100.))
        n19 = int(np.sum(kz & (f6 < 19.)))
        area = np.pi*(hi**2-lo**2)
        row = dict(zone_arcsec=[lo, hi], maglim_lf=maglim, mf_maglim=MF_MAGLIM, m_faint=m_faint, n_ms_fit=int(counts.sum()), alpha=alpha, alpha_err=alpha_err,
                   k=k, mass_ms_observed=m_ms_obs, mass_ms_extrapolated=m_extrap, n_evolved=n_evolved,
                   mass_evolved=n_evolved*iso.m_to, white_dwarfs=wd, n_f625w_lt19=n19, area_arcsec2=area,
                   mbins=mbins.tolist(), counts=counts.tolist())
        # luminous mass following the light per bright star, two bracketing variants
        for name in variants:
            for wname, w in wd.items():
                tot = m_ms_obs+m_extrap[name]+row["mass_evolved"]+w["n_wd"]*WD_MASS
                row.setdefault("mass_per_n19", {})[f"{name} | {wname}"] = tot/n19
        zones_out.append(row)
        # plots
        c = plt.get_cmap("viridis")(zi/(len(ZONES)-1))
        mids = .5*(mbins[:-1]+mbins[1:]); dm = np.diff(mbins)
        ax_mf.errorbar(mids, counts/dm/area*3600, np.sqrt(np.maximum(counts, 1))/dm/area*3600, fmt="o", color=c, ms=4,
                       label=f"{lo:g}-{hi:g}\": alpha {alpha:+.2f}+-{alpha_err:.2f}")
        mm = np.linspace(m_faint, iso.m_to, 50); ax_mf.plot(mm, k*mm**alpha/area*3600, color=c, lw=1)
        medges = np.arange(16., 25.01, .25)
        ax_lf.plot(.5*(medges[1:]+medges[:-1]), np.histogram(f6[kz], medges)[0]/area*3600/.25, color=c, label=f"{lo:g}-{hi:g}\"", drawstyle="steps-mid")
        ax_lf.axvline(maglim, color=c, ls=":", lw=.8)
        if zi == 0:
            ax_lf.axvline(MF_MAGLIM, color="black", ls="--", lw=1, label=f"MF fit limit {MF_MAGLIM:g}")
        ax_alpha.errorbar(np.sqrt(lo*hi)*args.distance*1e3/ARCSEC_PER_RAD, alpha, alpha_err, fmt="o", color=c)
        print(f"zone {lo:g}-{hi:g}\": N_MS(fit) {counts.sum()}, alpha {alpha:+.2f}+-{alpha_err:.2f}, m_faint {m_faint:.2f}, N_evolved {n_evolved}, N(F625W<19) {n19}, "
              f"M/N19 range {min(row['mass_per_n19'].values()):.1f}-{max(row['mass_per_n19'].values()):.1f} Msun")
    # CMD with isochrone
    sel = inzone & (np.random.default_rng(1).random(inzone.size) < .1)
    ax_cmd.scatter(col[sel], f6[sel], s=.3, color="0.6", alpha=.5)
    ax_cmd.plot(iso.m625-iso.m814, iso.m625, color="tab:red", lw=1.5, label=f"MIST {args.age:g} Gyr, [Fe/H] {feh_best:+.2f}, E(B-V) {args.ebv}")
    ax_cmd.errorbar(med, mid, xerr=hw, fmt=".", color="tab:blue", ms=3, label="ridge (median, 16-84%)")
    ax_cmd.set(xlim=(0, 1.6), ylim=(25, 16), xlabel="F625W - F814W", ylabel="F625W", title="(a) CMD (10% of stars, 10-250\") and isochrone"); ax_cmd.legend(fontsize=7)
    ax_lf.set(yscale="log", xlabel="F625W", ylabel="stars per arcmin$^2$ per mag", title="(b) luminosity functions (dotted: faint limits)"); ax_lf.legend(fontsize=7)
    ax_mf.set(xscale="log", yscale="log", xlabel=r"mass [$M_\odot$]", ylabel=r"dN/dm per arcmin$^2$", title="(c) present-day mass functions"); ax_mf.legend(fontsize=7)
    ax_alpha.set(xscale="log", xlabel="R [pc]", ylabel=r"slope $\alpha$ (dN/dm $\propto m^\alpha$)", title="(d) mass-function slope vs radius"); ax_alpha.axhline(-2.3, color="0.5", ls=":", lw=.8)
    ax_alpha.text(ax_alpha.get_xlim()[0], -2.3, " Salpeter/Kroupa -2.3", fontsize=7, va="bottom")
    # ratio panel: M/N19 per zone, all variants
    for zi, row in enumerate(zones_out):
        vals = list(row["mass_per_n19"].values()); xx = np.sqrt(row["zone_arcsec"][0]*row["zone_arcsec"][1])*args.distance*1e3/ARCSEC_PER_RAD
        ax_ratio.plot([xx]*len(vals), vals, "o", ms=4, color=plt.get_cmap("viridis")(zi/(len(ZONES)-1)))
    ax_ratio.set(xscale="log", xlabel="R [pc]", ylabel=r"luminous mass per star with F625W<19 [$M_\odot$]", title="(e) M(follows light)/N(F625W<19): all variants")
    # budget stacked for the middle zone
    row = zones_out[2]; base_key = "fitted slope | Kroupa/Salpeter IMF (-2.3)"
    parts = [("MS observed", row["mass_ms_observed"]), ("MS extrapolated (fitted slope)", row["mass_ms_extrapolated"]["fitted slope"]),
             ("evolved", row["mass_evolved"]), ("white dwarfs (IMF -2.3)", row["white_dwarfs"]["Kroupa/Salpeter IMF (-2.3)"]["n_wd"]*WD_MASS)]
    ax_budget.bar([p[0] for p in parts], [p[1]/row["n_f625w_lt19"] for p in parts])
    ax_budget.set(ylabel=r"$M_\odot$ per star with F625W<19", title=f"(f) mass budget, zone {row['zone_arcsec'][0]:g}-{row['zone_arcsec'][1]:g}\"")
    ax_budget.tick_params(axis="x", labelsize=7)
    fig.suptitle(f"omega Cen present-day mass function (oMEGACat, unflagged sample), MIST {args.age:g} Gyr")
    plot = ROOT/"plots"/f"mass_function_{args.tag}.png"
    fig.savefig(plot, dpi=140)
    out = ROOT/"results/mass_function"; out.mkdir(parents=True, exist_ok=True)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), age_gyr=args.age, distance_kpc=args.distance, ebv=args.ebv,
                  feh_adopted=feh_best, ridge_fits=fits, dm0=iso.dm0, av=iso.av, a625=iso.a625, a814=iso.a814, m_turnoff=iso.m_to,
                  wd_mass=WD_MASS, m_min_extrap=M_MIN_EXTRAP, zones=zones_out, plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (out/f"mass_function_{args.tag}.json").write_text(json.dumps(record, indent=1, default=float))
    print(plot)


if __name__ == "__main__":
    main()
