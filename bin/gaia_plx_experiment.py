#!/usr/bin/env python3
"""EXPERIMENT (2026-11-19; does not replace the production Gaia profile): Gaia EDR3 outer dispersions with parallax
as a third mixture dimension (adversarial review 5, scratchpad attack_gaia).

Per star the mixture likelihood becomes
    (1-f) N2(mu | cluster) g_c(plx) + f F(mu_abs) g_f(plx),
g_c = N(plx; plx_c, (k e_plx)^2) with plx_c, k measured from secure members at 300-900" (absorbs the parallax zero
point, which is common to cluster and field at fixed G), and g_f = the observed parallax density of FIELD stars
(template annulus, r > 1 deg, same G range and PM-error ceiling as the ring, absolute PM within --pmwin mas/yr of the
cluster's systemic PM, i.e. the field population that competes with the cluster). Because g_c does not depend on the
fitted parameters, the production estimator (outer_profile.dispersion_2d) is used unchanged with the field density
multiplied per star by g_f/g_c (exact). The field template is re-fetched from WSDB with parallaxes into a NEW file
(ocen_field_template_dr3_plx.ecsv); its proper-motion columns are checked to equal the production template.

Steps: (1) fetch/build the parallax template (once); (2) PM-only refit of every ring (must reproduce the product);
(3) PM+parallax fit; (4) parallax posterior-predictive check in ring 8; (5) injections in ring 8 (field = template
stars with their PM AND parallax, cluster with known sigma) for bias/pulls.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src python bin/gaia_plx_experiment.py [--ninj 100]
Output: results/gaia_plx/plx_experiment.json, plots/gaia_plx_experiment_20261119.png
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src"))
from astropy.table import Table
from scipy.ndimage import gaussian_filter1d
from ocen_dm.cluster import OCEN_DEC, OCEN_DISTANCE_KPC, OCEN_RA, OCEN_VSYS_KMS
from ocen_dm.paths import processed_dir
from ocen_dm.kinematics import outer_gaia as G
from ocen_dm.kinematics.outer_profile import dispersion_2d, load_members
from ocen_dm.kinematics.vb2021_replication import published_profile
from ocen_dm.selection import field_template as FT
from ocen_dm.kinematics.perspective import depth_dispersion, systemic_pm_field, systemic_velocity_vector

TPL = processed_dir()/"kinematics"/"ocen_field_template_dr3_plx.ecsv"
OUT = ROOT/"results/gaia_plx"
QUERY_PLX = FT.QUERY.replace("phot_g_mean_mag, ruwe,", "phot_g_mean_mag, ruwe, parallax, parallax_error, bp_rp,")


def build_template():
    """Same query, cuts and frame as field_template.build_product, plus parallax columns; new file."""
    import sqlutilpy as sqlutil
    q = QUERY_PLX.format(ra0=OCEN_RA, dec0=OCEN_DEC, r_in=FT.R_INNER_DEG, r_out=FT.R_OUTER_DEG, g_max=20.5)
    t0 = time.time(); d = sqlutil.get(q, asDict=True); print(f"fetched {len(d['ra'])} rows in {time.time()-t0:.1f}s")
    mu_sys = load_members().mu_sys
    v = systemic_velocity_vector(OCEN_RA, OCEN_DEC, mu_sys[0], mu_sys[1], OCEN_VSYS_KMS, OCEN_DISTANCE_KPC)
    ea, ed = systemic_pm_field(d["ra"], d["dec"], v, OCEN_DISTANCE_KPC)
    pa, pd = d["pmra"] - ea, d["pmdec"] - ed
    err = 0.5*(d["pmra_error"] + d["pmdec_error"])
    keep = (d["ruwe"] < 1.4) & (d["astrometric_excess_noise"] < 1.0) & (err < 1.5)
    x = (d["ra"] - OCEN_RA)*np.cos(np.radians(OCEN_DEC))*3600.0; y = (d["dec"] - OCEN_DEC)*3600.0
    t = Table({"r_arcsec": np.hypot(x, y)[keep], "sys_a": ea[keep], "sys_d": ed[keep], "mu_a": pa[keep],
               "mu_d": pd[keep], "pm_error": err[keep], "g_mag": d["phot_g_mean_mag"][keep],
               "parallax": d["parallax"][keep], "parallax_error": d["parallax_error"][keep], "bp_rp": d["bp_rp"][keep]})
    t.meta.update(product="ocen_field_template_dr3_plx (experiment)", query=q, n_fetched=int(len(keep)),
                  n_kept=int(keep.sum()), note="as ocen_field_template_dr3 plus parallax, parallax_error, bp_rp")
    t.write(TPL, format="ascii.ecsv", overwrite=True)
    return t


def check_template(t):
    """Same rows as the production template. Rows cannot be matched one to one (the production file has no ra/dec, and
    the cluster centre constant moved by ~0.001" since it was built), so compare each column's sorted values."""
    ref = FT.load_field_template()
    assert len(ref) == len(t), (len(ref), len(t))
    d = {c: float(np.max(np.abs(np.sort(np.asarray(ref[c], float)) - np.sort(np.asarray(t[c], float)))))
         for c in ("r_arcsec", "mu_a", "mu_d", "g_mag", "pm_error")}
    print(f"template check: {len(t)} rows as production; max |diff| of sorted columns {d}")
    assert d["g_mag"] < 1e-6 and d["pm_error"] < 1e-6 and d["mu_a"] < 2e-3 and d["mu_d"] < 2e-3 and d["r_arcsec"] < 0.01
    return d


class Ring:
    """Inputs of one production ring (outer_gaia.build_edr3_profile, raw errors)."""
    def __init__(self, base, cat, keep, err, lo, hi):
        self.lo, self.hi = lo, hi
        self.m = keep & (base.r_arcsec >= lo) & (base.r_arcsec < hi)
        gm = base.g_mag[self.m]; self.g_range = (float(gm.min())-0.5, float(gm.max())+0.5)
        self.err_max = float(np.max(err[self.m]))
        self.dens = FT.field_density_2d(err_max=self.err_max, g_range=self.g_range)
        self.sd = depth_dispersion(G._composite_tracer(OCEN_DISTANCE_KPC), base.r_arcsec[self.m]*OCEN_DISTANCE_KPC*1e3/206264.806,
                                   float(np.hypot(*base.mu_sys)), OCEN_DISTANCE_KPC)
        self.plx = np.asarray(cat["parallax"], float)[self.m]; self.eplx = np.asarray(cat["parallax_error"], float)[self.m]


def field_plx_density(tpl, ring, mu_sys, pmwin, bw=0.05):
    """Observed parallax density of template field stars matched to the ring (G, PM error, r > 1 deg) with absolute
    PM within pmwin of the systemic PM. Returns (callable density per mas, n stars used)."""
    t = tpl
    sel = ((np.asarray(t["r_arcsec"]) >= 3600) & (np.asarray(t["pm_error"]) <= ring.err_max)
           & (np.asarray(t["g_mag"]) >= ring.g_range[0]) & (np.asarray(t["g_mag"]) <= ring.g_range[1])
           & (np.hypot(np.asarray(t["mu_a"]) + np.asarray(t["sys_a"]) - mu_sys[0],
                       np.asarray(t["mu_d"]) + np.asarray(t["sys_d"]) - mu_sys[1]) < pmwin)
           & np.isfinite(np.asarray(t["parallax"], float)))
    p = np.asarray(t["parallax"], float)[sel]
    e = np.arange(-3, 6.0001, 0.01); h, _ = np.histogram(p, e); c = 0.5*(e[1:]+e[:-1])
    g = gaussian_filter1d(h.astype(float), bw/0.01) + 1e-6; g /= g.sum()*0.01
    return (lambda x: np.interp(x, c, g, left=g[0], right=g[-1])), int(sel.sum())


def cluster_plx(base, cat, keep):
    """plx_c and error scale k from secure members (membership_prob > 0.9) at 300-900"."""
    m = keep & (base.r_arcsec >= 300) & (base.r_arcsec < 900) & (np.asarray(cat["membership_prob"]) > 0.9)
    p, e = np.asarray(cat["parallax"], float)[m], np.asarray(cat["parallax_error"], float)[m]
    w = 1/e**2; pc = float(np.sum(w*p)/np.sum(w))
    z = (p - pc)/e; k = float(1.4826*np.median(np.abs(z - np.median(z))))   # robust scale
    return pc, k, int(m.sum())


def fit(base, ring, dens, mask=None):
    m = ring.m if mask is None else mask
    return dispersion_2d(base, m, dens, depth_var=ring.sd**2 if mask is None else ring.sd_m**2, field_at=base.absolute_pm)


class Ratio:
    """field density x g_f(plx)/g_c(plx), per star, in the order dispersion_2d passes the masked stars."""
    def __init__(self, dens, ratio): self.dens, self.ratio = dens, ratio
    def __call__(self, a, d): return self.dens(a, d)*self.ratio


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ninj", type=int, default=100)
    ap.add_argument("--pmwin", type=float, default=2.0); ap.add_argument("--refetch", action="store_true")
    a = ap.parse_args(); OUT.mkdir(parents=True, exist_ok=True)
    tpl = build_template() if (a.refetch or not TPL.exists()) else Table.read(TPL)
    check_template(tpl)
    D = OCEN_DISTANCE_KPC
    cat = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv"); qf = np.asarray(cat["quality_flag"], int)
    rp, sp = published_profile(); base = load_members(exact=True, distance_kpc=D)
    err = 0.5*(base.err_r + base.err_t)
    keep = ((qf & G.QUALITY_BIT) > 0) & (err < 0.4*np.interp(base.r_arcsec, rp, sp))
    pc, k, nsec = cluster_plx(base, cat, keep)
    print(f"cluster parallax from {nsec} secure members at 300-900\": plx_c {pc:.4f} mas, error scale k {k:.3f}")
    prod = G.load_edr3_profile(); rows = []
    for i, (lo, hi) in enumerate(zip(G.DEFAULT_EDGES[:-1], G.DEFAULT_EDGES[1:])):
        ring = Ring(base, cat, keep, err, lo, hi)
        if ring.m.sum() < 25:
            continue
        o0 = fit(base, ring, ring.dens)
        gf, nf = field_plx_density(tpl, ring, base.mu_sys, a.pmwin)
        gc = np.exp(-0.5*((ring.plx-pc)/(k*ring.eplx))**2)/(np.sqrt(2*np.pi)*k*ring.eplx)
        o1 = fit(base, ring, Ratio(ring.dens, gf(ring.plx)/np.maximum(gc, 1e-300)))
        pr = prod[len(rows)]
        row = dict(ring=i, lo=lo, hi=hi, n=int(ring.m.sum()), n_field_plx=nf,
                   sr_prod=float(pr["sigma_pmr"]), st_prod=float(pr["sigma_pmt"]),
                   sr0=o0["sigma_r"], st0=o0["sigma_t"], esr0=o0["sigma_r_err"], est0=o0["sigma_t_err"], f0=o0["f"],
                   sr1=o1["sigma_r"], st1=o1["sigma_t"], esr1=o1["sigma_r_err"], est1=o1["sigma_t_err"], f1=o1["f"],
                   ncl0=float(np.sum(o0["p_member"])), ncl1=float(np.sum(o1["p_member"])))
        rows.append(row)
        print(f"ring {lo:5.0f}-{hi:5.0f} N {row['n']:6d} f {o0['f']:.3f}->{o1['f']:.3f}  sigma_R {o0['sigma_r']:.4f}+-{o0['sigma_r_err']:.4f} -> "
              f"{o1['sigma_r']:.4f}+-{o1['sigma_r_err']:.4f}  sigma_T {o0['sigma_t']:.4f}+-{o0['sigma_t_err']:.4f} -> "
              f"{o1['sigma_t']:.4f}+-{o1['sigma_t_err']:.4f}  (product {row['sr_prod']:.4f}/{row['st_prod']:.4f}; field-plx stars {nf})", flush=True)
        assert abs(o0["sigma_r"]-row["sr_prod"]) < 2e-3 and abs(o0["sigma_t"]-row["st_prod"]) < 2e-3
        last = (ring, o1, gf, gc)
    # parallax PPC in the outermost ring
    ring, o1, gf, gc = last
    p1 = np.asarray(o1["p_member"])
    # restricted to the stars the field parallax model describes: absolute PM within pmwin of the systemic PM
    s8 = base.select(ring.m); aa, dd = s8.absolute_pm
    win = np.hypot(aa - base.mu_sys[0], dd - base.mu_sys[1]) < a.pmwin
    e = np.arange(-1.0, 2.0001, 0.05); cc = 0.5*(e[1:]+e[:-1]); obs = np.histogram(ring.plx[win], e)[0]
    ncw = float(np.sum(p1[win])); nfw = float(win.sum()) - ncw
    expf = nfw*gf(cc)*0.05
    expc = sum(pi*np.exp(-0.5*((cc-pc)/(k*s))**2)/(np.sqrt(2*np.pi)*k*s) for pi, s in zip(p1[win], ring.eplx[win]))*0.05
    print(f"ring 8 parallax PPC window: {win.sum()} stars, {ncw:.1f} members, {nfw:.1f} field")
    chi = (obs-expf-expc)/np.sqrt(np.maximum(expf+expc, 1)); print(f"ring 8 parallax PPC: chi2 {np.sum(chi**2):.1f} / {len(cc)}")
    # injections in the outermost ring: field drawn from matched template stars (PM and parallax jointly)
    rng = np.random.default_rng(7)
    t = tpl
    pool = ((np.asarray(t["r_arcsec"]) >= 3600) & (np.asarray(t["pm_error"]) <= ring.err_max)
            & (np.asarray(t["g_mag"]) >= ring.g_range[0]) & (np.asarray(t["g_mag"]) <= ring.g_range[1])
            & np.isfinite(np.asarray(t["parallax"], float)))
    ia = np.flatnonzero(pool)
    idx = np.flatnonzero(ring.m); nst = len(idx); ncl = int(round(np.sum(p1)))
    truth = (o1["sigma_r"], o1["sigma_t"])
    inj = []
    s0 = base.select(ring.m)
    for it in range(a.ninj):
        b = base.select(ring.m)                      # copy positions/errors of the real ring stars
        cl = rng.random(nst) < ncl/nst
        # field stars: PM and parallax from template stars (absolute PM -> residual at this star's systemic field)
        j = rng.choice(ia, nst)
        fa_abs = np.asarray(t["mu_a"])[j] + np.asarray(t["sys_a"])[j]; fd_abs = np.asarray(t["mu_d"])[j] + np.asarray(t["sys_d"])[j]
        phi = np.asarray(b.phi); cp, sp_ = np.sin(phi), np.cos(phi)
        vr = rng.normal(o1["mean_r"], truth[0], nst); vt = rng.normal(o1["mean_t"], truth[1], nst)
        ca, cd = vr*cp - vt*sp_, vr*sp_ + vt*cp
        ea, ed = np.asarray(b.err_a), np.asarray(b.err_d)
        na, nd = rng.normal(0, ea), rng.normal(0, ed)
        mua = np.where(cl, ca + na, fa_abs - np.asarray(b.sys_a)); mud = np.where(cl, cd + nd, fd_abs - np.asarray(b.sys_d))
        plx = np.where(cl, rng.normal(pc, k*ring.eplx), np.asarray(t["parallax"], float)[j])
        from dataclasses import replace as _rep
        b = _rep(b, mu_a=mua, mu_d=mud, mu_r=mua*cp + mud*sp_, mu_t=-mua*sp_ + mud*cp, err_corr=np.zeros(nst))
        gci = np.exp(-0.5*((plx-pc)/(k*ring.eplx))**2)/(np.sqrt(2*np.pi)*k*ring.eplx)
        mm = np.ones(nst, bool)
        q0 = dispersion_2d(b, mm, ring.dens, depth_var=0.0, field_at=b.absolute_pm)
        q1 = dispersion_2d(b, mm, Ratio(ring.dens, gf(plx)/np.maximum(gci, 1e-300)), depth_var=0.0, field_at=b.absolute_pm)
        inj.append([q0["sigma_r"], q0["sigma_t"], q0["sigma_r_err"], q0["sigma_t_err"],
                    q1["sigma_r"], q1["sigma_t"], q1["sigma_r_err"], q1["sigma_t_err"]])
    inj = np.array(inj)
    res = {}
    for lab, c in (("PM", 0), ("PM+plx", 4)):
        pr_ = (inj[:, c]-truth[0])/inj[:, c+2]; pt_ = (inj[:, c+1]-truth[1])/inj[:, c+3]
        res[lab] = dict(bias_R=float(np.mean(inj[:, c]-truth[0])), bias_T=float(np.mean(inj[:, c+1]-truth[1])),
                        pull_mean_R=float(pr_.mean()), pull_sd_R=float(pr_.std()), pull_mean_T=float(pt_.mean()), pull_sd_T=float(pt_.std()),
                        rms_R=float(np.std(inj[:, c])), rms_T=float(np.std(inj[:, c+1])))
        print(f"injections ({a.ninj}) {lab}: bias R {res[lab]['bias_R']:+.4f} T {res[lab]['bias_T']:+.4f}; pulls R {pr_.mean():+.2f}+-{pr_.std():.2f} "
              f"T {pt_.mean():+.2f}+-{pt_.std():.2f}; scatter R {res[lab]['rms_R']:.4f} T {res[lab]['rms_T']:.4f}")
    json.dump(dict(plx_c=pc, k=k, n_secure=nsec, pmwin=a.pmwin, rows=rows, ppc_chi2=float(np.sum(chi**2)), ppc_nbin=len(cc),
                   injections=res, truth=truth, note="EXPERIMENT; production product unchanged"),
              open(OUT/"plx_experiment.json", "w"), indent=1)
    # plot
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    r = np.array([0.5*(x["lo"]+x["hi"]) for x in rows])
    for c, (k0, k1, e0, e1, lab) in enumerate((("sr0", "sr1", "esr0", "esr1", "R"), ("st0", "st1", "est0", "est1", "T"))):
        v0 = np.array([x[k0] for x in rows]); v1 = np.array([x[k1] for x in rows])
        s0_ = np.array([x[e0] for x in rows]); s1_ = np.array([x[e1] for x in rows])
        ax[0].errorbar(r*(0.985 if c == 0 else 1.0), v0, s0_, fmt="o", mfc="none", color=f"C{c}", label=f"sigma_{lab} PM only (production)")
        ax[0].errorbar(r*1.015, v1, s1_, fmt="s", color=f"C{c}", ms=4, label=f"sigma_{lab} PM + parallax")
        ax[1].plot(r, (v1-v0)/s0_, "o-", color=f"C{c}", label=f"sigma_{lab}")
    ax[0].set(xscale="log", xlabel="R (arcsec)", ylabel="dispersion (mas/yr)", title="Gaia rings: PM-only vs PM+parallax mixture")
    ax[0].legend(fontsize=7)
    ax[1].axhline(0, c="k", lw=.5); ax[1].set(xscale="log", xlabel="R (arcsec)", ylabel="(PM+plx - PM) / error", title="shift per ring")
    ax[1].legend(fontsize=7)
    ax[2].errorbar(cc, obs, np.sqrt(np.maximum(obs, 1)), fmt="ko", ms=3, label="observed")
    ax[2].plot(cc, expf, "C7-", label="field (template parallaxes)"); ax[2].plot(cc, expc, "C0-", label="cluster N(plx_c, k e)")
    ax[2].plot(cc, expf+expc, "r-", lw=2, label=f"total: chi2 {np.sum(chi**2):.0f}/{len(cc)}")
    ax[2].set(xlabel="parallax (mas)", ylabel="stars / 0.05 mas", title=f"ring {ring.lo:.0f}-{ring.hi:.0f}\": parallax check"); ax[2].legend(fontsize=7)
    fig.tight_layout(); f = ROOT/"plots/gaia_plx_experiment_20261119.png"; fig.savefig(f, dpi=100); print(f)


if __name__ == "__main__":
    main()
