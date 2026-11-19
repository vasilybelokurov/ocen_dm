#!/usr/bin/env python3
"""Plot the dark-mass limit products (docs/PLAN_2026-11-19_limits.md):
(a) null bootstrap: t = 2 Delta lnpost (NFW - K1) on K1 mocks, vs the real-data value and the Chernoff 2.71;
(b) injections: q = 2[max lnpost - lnpost(M_true, r_s = 300)] for injected NFW halos (one-sided);
(c) spherical profile q(M(<50 pc)) per r_s, base and rotation-marginalised.
Missing inputs are skipped. Usage: PYTHONPATH=src:bin python bin/plot_limits.py  ->  plots/dm_limits_20261119.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; R = ROOT/"results/limits"


def jl(p):
    return [json.loads(l) for l in open(p)] if p.exists() else []


def main():
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    t = np.array([r["t"] for r in jl(R/"boot_null.jsonl")])
    K = json.loads((R/"sph_local_K1_inner.json").read_text())["best"]["lnpost"]
    N = json.loads((R/"sph_local_NFW_inner.json").read_text())["best"]["lnpost"]
    treal = 2*(N - K)
    if len(t):
        ax[0].hist(t, bins=np.arange(0, max(t.max(), treal)+0.5, 0.25), color="0.6")
        ax[0].axvline(treal, color="C3", label=f"real data t = {treal:.2f}")
        ax[0].axvline(np.percentile(t, 95), color="k", ls="--", label=f"null 95% = {np.percentile(t, 95):.2f}")
        ax[0].axvline(2.71, color="C0", ls=":", label="Chernoff 95% = 2.71")
        ax[0].set(xlabel=r"$t = 2\Delta\ln p$ (NFW $-$ K1)", ylabel="mocks",
                  title=f"null (K1 truth), {len(t)} mocks; {int((t >= treal).sum())} with t >= real")
        ax[0].legend(fontsize=8)
    for i, M in enumerate(("1e6", "2e6")):
        q = np.array([r["q"] for r in jl(R/f"boot_inject_{float(M):.0e}.jsonl")])
        if len(q):
            ax[1].hist(q, bins=np.arange(0, max(q.max(), 3)+0.5, 0.25), histtype="step", lw=2, color=f"C{i}",
                       label=f"M100 = {M}: 95% = {np.percentile(q, 95):.2f} ({len(q)} mocks)")
    ax[1].axvline(2.71, color="k", ls=":", label="2.71")
    ax[1].set(xlabel=r"$q(M_{\rm true})$, $r_s$ = 300 pc", ylabel="mocks", title="NFW injections (coverage)")
    ax[1].legend(fontsize=8)
    for f, ls in (("sph_profile.json", "-"), ("sph_profile_rotmarg.json", "--")):
        if not (R/f).exists():
            continue
        d = json.loads((R/f).read_text())
        for j, (rs, rows) in enumerate(d["grid"].items()):
            m = [r["M50"] for r in rows]; q = [2*(d["lnpost_max"] - r["lnpost"]) for r in rows]
            ax[2].plot(m, q, ls, color=f"C{j}", label=f"r_s = {float(rs):g} pc" + (" (rot. marg.)" if ls == "--" else ""))
    ax[2].axhline(2.71, color="k", ls=":")
    ax[2].set(xscale="log", xlabel=r"$M_{\rm DM}(<50$ pc) [M$_\odot$]", ylabel=r"$2[\max\ln p - \ln p]$",
              ylim=(-0.3, 15), title="spherical profile likelihood")
    ax[2].legend(fontsize=7)
    fig.tight_layout(); out = ROOT/"plots/dm_limits_20261119.png"; fig.savefig(out, dpi=120); print(out)


if __name__ == "__main__":
    main()
