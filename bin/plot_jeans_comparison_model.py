#!/usr/bin/env python3
"""The rung-2 Jeans comparison model (no DM, free black hole) against the DF base model.

Reads only saved products: the Jeans run's posterior (posterior.ecsv) and its 300 saved
mass-profile draws (profiles.npz); the DF base model is rebuilt once from its saved best fit.
Figures (plots/jeans_<name>.png):
  components  enclosed mass by component (Jeans posterior median and 16-84 per cent;
              DF base model dashed), circular velocity, and the dark (non-stellar) share
              M_dark(<r)/M_tot(<r) with dark = remnants + black hole (Jeans) or remnants (DF);
  beta        anisotropy beta(r): Jeans posterior median, 16-84 and 2-98 per cent, its
              transition radii, and the DF base model; radial ranges of the data shaded;
  corner      posterior of the black-hole mass and the parameters it trades against;
  centre      the central bins (R < 2 pc): data, Jeans and DF predictions and residuals.
Numbers go to results/plot_data/jeans_<name>.json.

Usage: python bin/plot_jeans_comparison_model.py --name rung2_K1_turnover \
    --jeans results/fits/rung2_K1_turnover \
    --df results/df/dftwo_counts_area_20261013::observed_no_halo_start3
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.compact_recovery import mass_profiles, refined_config
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model, model_config_from_dict
from ocen_dm.kinematics.run_io import read_data_snapshot

HST_MUSE = (0., 9.5)       # pc, projected range of the HST/MUSE bins
GAIA = (9.4, 63.)          # pc, projected range of the Gaia bins
JC, DC = "C1", "C0"        # Jeans, DF colours


def beta_turnover(r, b0, bmid, binf, r1, dr):
    """TurnoverAnisotropy (src/ocen_dm/kinematics/anisotropy.py): weights f_i = r^2/(r^2+r_i^2)."""
    f1 = r**2/(r**2+r1**2)
    f2 = r**2/(r**2+(r1+dr)**2)
    return b0*(1-f1)+bmid*(f1-f2)+binf*f2


def q(a, axis=0):
    return np.percentile(a, [2, 16, 50, 84, 98], axis=axis)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--name", required=True)
    p.add_argument("--jeans", required=True)
    p.add_argument("--df", required=True, help="BATCH::JOB")
    a = p.parse_args()
    jd = Path(a.jeans)
    js = json.loads((jd/"summary.json").read_text())
    post = Table.read(jd/"posterior.ecsv")
    prof = dict(np.load(jd/"profiles.npz"))
    batch, job = a.df.split("::")
    ds = json.loads((Path(batch)/"fits"/job/"summary.json").read_text())
    cfg = refined_config(model_config_from_dict(ds["best"]["config"]))
    dfm = build_df_model(cfg)

    r = prof["r"]                                       # pc, 0.1-500
    rd = r[r < cfg.numerics.r_max/100*0.999]
    dmp = {k: np.asarray(v) for k, v in mass_profiles(dfm, rd).items() if isinstance(v, list)}
    dfbeta = dfm.intrinsic_moments(rd)["beta"]
    G = 4.300917e-3                                     # pc (km/s)^2 / Msun

    # ---------------------------------------------------------------- components
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.2), constrained_layout=True)
    for key, lab, c in (("M_stars", "stars", "C0"), ("M_remnants", "remnants", "C3"), ("M_imbh", "black hole", "C4"),
                        ("M_total", "total", "k")):
        lo2, lo, med, hi, hi2 = q(prof[key])
        ax[0].fill_between(r, lo, hi, color=c, alpha=0.25, lw=0)
        ax[0].plot(r, med, color=c, lw=2 if key == "M_total" else 1.4, label=f"Jeans: {lab}")
    for key, lab, c in (("stars", "stars", "C0"), ("remnants", "remnants", "C3"), ("total", "total", "k")):
        ax[0].plot(rd, dmp[key], color=c, ls="--", lw=1.2, label=f"DF: {lab}")
    ax[0].set(xscale="log", yscale="log", ylim=(1e3, 5e6), ylabel=r"$M(<r)$ [M$_\odot$]",
              title="(a) enclosed mass (Jeans: median, 16-84%; DF dashed)")
    lo2, lo, med, hi, hi2 = q(prof["v_circ"])
    ax[1].fill_between(r, lo, hi, color=JC, alpha=0.3, lw=0)
    ax[1].plot(r, med, color=JC, label="Jeans (with black hole)")
    ax[1].plot(rd, np.sqrt(G*dmp["total"]/rd), color=DC, ls="--", label="DF base model")
    ax[1].set(xscale="log", ylabel=r"$v_c$ [km/s]", title="(b) circular velocity")
    dark = (prof["M_remnants"]+prof["M_imbh"])/prof["M_total"]
    lo2, lo, med, hi, hi2 = q(dark)
    ax[2].fill_between(r, lo, hi, color=JC, alpha=0.3, lw=0)
    ax[2].plot(r, med, color=JC, label="Jeans: (remnants + black hole) / total")
    bhf = q(prof["M_imbh"]/prof["M_total"])[2]
    ax[2].plot(r, bhf, color="C4", lw=1, label="Jeans: black hole / total")
    ax[2].plot(rd, dmp["remnants"]/dmp["total"], color=DC, ls="--", label="DF: remnants / total")
    ax[2].set(xscale="log", ylim=(0, 1), ylabel="fraction of enclosed mass", title="(c) non-stellar share of M(<r)")
    for x in ax:
        x.set_xlabel("r [pc]")
        x.axvspan(r[0], 0.1, color="0.92", lw=0)
        x.axvspan(63, r[-1], color="0.92", lw=0)
        x.set_xlim(r[0], r[-1])
        x.legend(fontsize=7.5)
    fig.suptitle("Jeans comparison model (rung 2, no DM, free black hole) against the DF base model: mass components")
    fig.savefig(ROOT/"plots"/f"jeans_{a.name}_components.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- beta
    rr = np.geomspace(0.01, 200, 300)
    B = np.array([beta_turnover(rr, *[float(row[k]) for k in ("beta_0", "beta_mid", "beta_inf", "r_beta", "delta_r_beta")])
                  for row in post[np.random.default_rng(1).choice(len(post), 3000, replace=False)]])
    b2, b16, b50, b84, b98 = q(B)
    ml = js["parameters"]
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    ax.axvspan(*HST_MUSE, color="C2", alpha=0.08, lw=0, label="HST/MUSE projected range")
    ax.axvspan(*GAIA, color="C0", alpha=0.08, lw=0, label="Gaia projected range")
    ax.fill_between(rr, b2, b98, color=JC, alpha=0.15, lw=0)
    ax.fill_between(rr, b16, b84, color=JC, alpha=0.35, lw=0)
    ax.plot(rr, b50, color=JC, lw=2, label="Jeans turnover: median, 16-84%, 2-98%")
    ax.plot(rd, dfbeta, color=DC, lw=2, ls="--", label="DF base model")
    r1 = ml["r_beta"]["p50"]
    r2 = r1+ml["delta_r_beta"]["p50"]
    for x, lab in ((r1, "$r_1$"), (r2, "$r_2$")):
        ax.axvline(x, color=JC, ls=":", lw=1)
        ax.text(x*1.05, 0.55, lab, color=JC)
    ax.axhline(-0.5, color="0.4", ls="-.", lw=0.8)
    ax.text(0.012, -0.47, r"$\beta_0\leq-1/2$ required with a black hole and a cored tracer", fontsize=8, color="0.3")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set(xscale="log", xlim=(0.01, 200), ylim=(-1.05, 0.7), xlabel="r [pc]", ylabel=r"$\beta = 1-\sigma_t^2/\sigma_r^2$",
           title="Anisotropy: Jeans turnover profile against the DF base model")
    ax.legend(fontsize=8, loc="lower left")
    fig.savefig(ROOT/"plots"/f"jeans_{a.name}_beta.png", dpi=140)
    plt.close(fig)

    # ---------------------------------------------------------------- corner
    import corner
    names = ["M_bh", "beta_0", "M_rem", "a_rem", "M_star", "r_beta"]
    labels = [r"$\log_{10}M_{\rm BH}$", r"$\beta_0$", r"$\log_{10}M_{\rm rem}$", r"$a_{\rm rem}$ [pc]",
              r"$\log_{10}M_\star$", r"$r_1$ [pc]"]
    X = np.column_stack([np.log10(post[n]) if n in ("M_bh", "M_rem", "M_star") else np.asarray(post[n]) for n in names])
    fig = corner.corner(X, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".3g",
                        title_kwargs=dict(fontsize=9), label_kwargs=dict(fontsize=10), smooth=1.0, color=JC)
    fig.suptitle(f"Jeans comparison model: posterior ({len(post)} samples)", y=1.01)
    fig.savefig(ROOT/"plots"/f"jeans_{a.name}_corner.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    corr = np.corrcoef(X.T)

    # ---------------------------------------------------------------- centre
    # Jeans predictions at this run's own maximum-likelihood point (chi2 replayed against the saved value);
    # before 2026-10-16 this read the rung2_K1_turnover column of a ladder table for every run.
    from ocen_dm.kinematics.report import load_run_metadata, problem_for
    jrun = load_run_metadata(jd)
    jprob = problem_for(jrun)
    jx = np.array([jrun["summary"]["parameters"][n]["ml"] for n in jprob.family.names])
    jpred = jprob.predict(jx)
    assert np.isclose(sum(v[0] for v in jprob.chi2(jx).values()), jrun["summary"]["chi2_ml_total"], atol=.011)
    od = Path(batch)/"observed"
    rec = json.loads((od/"problem.json").read_text())
    counts = [CountProfile.from_dict(json.loads((od/f"{n}.json").read_text())) for n in rec["counts"]]
    prob = DFJointProblem(read_data_snapshot(od, rec["data_snapshot"]), None, counts)
    ev = prob.evaluate(cfg, dfm)
    pc = 5.43e3/206264.806
    fig, ax = plt.subplots(2, 3, figsize=(17, 7.5), sharex=True, constrained_layout=True,
                           gridspec_kw=dict(height_ratios=[2.2, 1]))
    central = {}
    for i, (name, lab, unit) in enumerate((("hst_pm_radial_ours", "HST PM radial", "mas/yr"),
                                           ("hst_pm_tangential_ours", "HST PM tangential", "mas/yr"),
                                           ("muse_los_dispersion", "MUSE line of sight", "km/s"))):
        prf = [x for x in prob.data.profiles if x.name == name][0]
        R = prf.r*pc
        k = R < 2.
        jm = np.asarray(jpred[name], float)
        dm = ev["predictions"][name]
        ax[0, i].errorbar(R[k], prf.value[k], yerr=[prf.err_lo[k], prf.err_hi[k]], fmt="o", color="k", ms=4, lw=0.8, label="data")
        ax[0, i].plot(R[k], jm[k], color=JC, marker="s", ms=3, label="Jeans + black hole")
        ax[0, i].plot(R[k], dm[k], color=DC, marker="D", ms=3, ls="--", label="DF base model")
        rj = (prf.value-jm)/np.where(jm > prf.value, prf.err_hi, prf.err_lo)
        rdf = (prf.value-dm)/np.where(dm > prf.value, prf.err_hi, prf.err_lo)
        ax[1, i].plot(R[k], rj[k], color=JC, marker="s", ms=3)
        ax[1, i].plot(R[k], rdf[k], color=DC, marker="D", ms=3, ls="--")
        ax[1, i].axhline(0, color="0.5", lw=0.8)
        ax[0, i].set(title=f"{lab} (R < 2 pc)", ylabel=f"dispersion [{unit}]")
        ax[1, i].set(xscale="log", xlabel="R [pc]", ylabel="(data-model)/err")
        ax[0, i].legend(fontsize=8)
        central[name] = dict(n=int(k.sum()), chi2_jeans=float(np.sum(rj[k]**2)), chi2_df=float(np.sum(rdf[k]**2)))
    fig.suptitle("Central bins: the region the black hole affects")
    fig.savefig(ROOT/"plots"/f"jeans_{a.name}_centre.png", dpi=140)
    plt.close(fig)

    out = dict(jeans=str(jd), df=a.df, n_posterior=len(post), profile_draws=int(prof["M_total"].shape[0]),
               parameters={k: {kk: v[kk] for kk in ("ml", "p16", "p50", "p84", "prior")} for k, v in ml.items()},
               r1_r2_pc=[r1, r2], posterior_correlations=dict(names=names, matrix=np.round(corr, 3).tolist()),
               central_chi2=central,
               enclosed_mass_median={str(x): {k: float(np.interp(np.log(x), np.log(r), np.median(prof[k], axis=0)))
                                              for k in ("M_stars", "M_remnants", "M_imbh", "M_total")} for x in (0.3, 1, 3, 10, 30, 100)},
               beta_median={str(x): float(np.interp(np.log(x), np.log(rr), b50)) for x in (0.03, 0.3, 1, 3, 10, 30, 100)},
               logz=js["logz"], chi2_ml={k: v["chi2"] for k, v in js["chi2_ml"].items()})
    (ROOT/"results/plot_data"/f"jeans_{a.name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("r1_r2_pc", "posterior_correlations", "central_chi2", "enclosed_mass_median", "beta_median")}, indent=1))


if __name__ == "__main__":
    main()
