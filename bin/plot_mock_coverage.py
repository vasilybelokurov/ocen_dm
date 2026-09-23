#!/usr/bin/env python3
"""Recovery of mass components in the realistic mocks (coverage test).

Collects every finished mock fit under results/df/cov_*_<TAG> and
results/df/mock{A,B}_*_<TAG> (truth rho20 = 0 or injected), and plots:
 (a) recovered rho20 against the truth, one point per noise seed, with the
     profiled r_s as colour; (b) recovered total enclosed mass relative to the
     truth as a function of radius (one curve per mock); (c) recovered DM mass
     inside 20 and 63 pc against the truth; (d) fractional errors of the
     stellar and remnant masses. The observed-data rho20 profile crossings
     (Delta = 1, 3.84) are marked in (a) for scale.

Usage: python bin/plot_mock_coverage.py [--tag 20260923]
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def m_at(prof, key, r):
    return np.interp(r, prof["r_pc"], prof[key])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    args = parser.parse_args()
    mocks = []
    for b in sorted(glob.glob(str(ROOT/f"results/df/cov_*_{args.tag}"))+glob.glob(str(ROOT/f"results/df/mock*_{args.tag}"))):
        if not os.path.exists(f"{b}/mocks/mock/mock.json"):
            continue
        t = read(f"{b}/mocks/mock/mock.json")
        fits = [read(f) for f in glob.glob(f"{b}/fits/*/summary.json")]
        fits = [s for s in fits if "best" in s]
        if not fits:
            continue
        s = min(fits, key=lambda s: s["refined_score"])      # best converged start of this mock
        c = s["best"]["config"]; tc = t["config"]
        mocks.append(dict(name=os.path.basename(b), seed=t.get("seed"), truth_rho20=tc["matter"]["rho20"], truth_rs=tc["matter"]["r_s"],
                          fit_rho20=c["matter"]["rho20"], fit_rs=c["matter"]["r_s"],
                          mstar_err=c["M_star"]/tc["M_star"]-1, mrem_err=c["matter"]["M_rem"]/tc["matter"]["M_rem"]-1,
                          D=c["distance_kpc"], truth_D=tc["distance_kpc"], n_starts=len(fits),
                          truth_prof=t["profiles"], fit_prof=s["refined_profiles"],
                          chi2_kin=s["gates"]["chi2_kinematic"], nominal_chi2_kin=t["nominal_chi2_kin"]))
    if not mocks:
        raise SystemExit("no finished mocks")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    # (a) rho20 recovery
    ax = axes[0, 0]
    truths = sorted({m["truth_rho20"] for m in mocks})
    for m in mocks:
        ax.scatter(m["truth_rho20"]+np.random.default_rng(m["seed"] or 0).uniform(-.04, .04), m["fit_rho20"],
                   c=[np.log10(m["fit_rs"])], vmin=np.log10(5), vmax=np.log10(500), cmap=cmap, s=70, edgecolor="black", zorder=3)
    lim = max(1.6, max(m["fit_rho20"] for m in mocks)+.2)
    ax.plot([0, lim], [0, lim], color="0.5", ls="--", lw=1, label="recovered = truth")
    for y, txt in ((0.5, r"observed profile $\Delta=1$"), (1.7, r"observed profile $\Delta=3.84$")):
        ax.axhline(y, color="tab:red", ls=":", lw=1); ax.text(lim, y, txt, ha="right", va="bottom", fontsize=7, color="tab:red")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(np.log10(5), np.log10(500)))
    fig.colorbar(sm, ax=ax, label=r"$\log_{10}$ profiled $r_s$ [pc] (bounds 5-500)")
    ax.set(xlabel=r"truth $\rho_{20}$ [$M_\odot$ pc$^{-3}$]", ylabel=r"recovered $\rho_{20}$", xlim=(-.15, max(truths)+.4), ylim=(-.05, lim),
           title=f"(a) rho20 recovery, {len(mocks)} mocks (one best start each)")
    ax.legend(fontsize=7, loc="upper left")
    # (b) total mass profile ratio
    ax = axes[0, 1]
    r = np.geomspace(.5, 80., 64)
    for m in mocks:
        ratio = m_at(m["fit_prof"], "total", r)/m_at(m["truth_prof"], "total", r)
        ax.plot(r, ratio, color="tab:blue" if m["truth_rho20"] == 0 else "tab:orange", alpha=.8,
                label=("truth: no halo" if m["truth_rho20"] == 0 else f"truth: rho20 = {m['truth_rho20']:g}") )
    h, l = ax.get_legend_handles_labels(); uniq = dict(zip(l, h)); ax.legend(uniq.values(), uniq.keys(), fontsize=7)
    ax.axhline(1, color="black", lw=.8); ax.axvspan(63, 80, color="0.92")
    for y in (.95, 1.05): ax.axhline(y, color="0.6", ls=":", lw=.8)
    ax.set(xscale="log", xlabel="r [pc]", ylabel="recovered / true total enclosed mass", title="(b) total mass profile recovery (grey: beyond data)")
    # (c) DM mass inside 20 and 63 pc
    ax = axes[1, 0]
    for m in mocks:
        for rr, mk in ((20., "o"), (63., "s")):
            ax.scatter(m_at(m["truth_prof"], "halo", rr)+1., m_at(m["fit_prof"], "halo", rr)+1., marker=mk, s=60, edgecolor="black",
                       color="tab:blue" if m["truth_rho20"] == 0 else "tab:orange", label=f"M_DM(<{rr:g} pc)" if m is mocks[0] else None)
    ax.plot([1, 3e6], [1, 3e6], color="0.5", ls="--", lw=1)
    ax.set(xscale="log", yscale="log", xlim=(1, 3e6), ylim=(1, 3e6), xlabel=r"true $M_{\rm DM}(<r)+1$ [$M_\odot$]", ylabel=r"recovered $M_{\rm DM}(<r)+1$",
           title="(c) DM mass inside 20 pc (o) and 63 pc (s)")
    ax.legend(fontsize=7)
    # (d) stellar / remnant errors
    ax = axes[1, 1]
    x = np.arange(len(mocks))
    ax.bar(x-.2, [100*m["mstar_err"] for m in mocks], .4, label=r"$M_\star$")
    ax.bar(x+.2, [100*m["mrem_err"] for m in mocks], .4, label=r"$M_{\rm rem}$")
    ax.set_xticks(x); ax.set_xticklabels([f"{'A' if m['truth_rho20']==0 else 'B'}{m['seed']}" for m in mocks], fontsize=7)
    ax.axhline(0, color="black", lw=.8)
    ax.set(ylabel="recovered / truth - 1 [%]", title="(d) stellar and remnant mass errors (A: no halo; B: injected halo)")
    ax.legend(fontsize=7)
    fig.suptitle("Realistic mocks (published kinematic noise, Poisson counts): what the fit recovers")
    plot = ROOT/"plots"/f"mock_coverage_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  mocks=[{k: v for k, v in m.items() if k not in ("truth_prof", "fit_prof")} for m in mocks],
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"mock_coverage_{args.tag}.json").write_text(json.dumps(record))
    for m in mocks:
        print(f"{m['name']:28s} truth {m['truth_rho20']:.2f} -> {m['fit_rho20']:.3f} (r_s {m['fit_rs']:.0f})  M* {100*m['mstar_err']:+.1f}%  Mrem {100*m['mrem_err']:+.1f}%  "
              f"Mtot(<20) {100*(m_at(m['fit_prof'],'total',20.)/m_at(m['truth_prof'],'total',20.)-1):+.1f}%  Mtot(<63) {100*(m_at(m['fit_prof'],'total',63.)/m_at(m['truth_prof'],'total',63.)-1):+.1f}%")
    print(plot)


if __name__ == "__main__":
    main()
