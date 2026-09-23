#!/usr/bin/env python3
"""Plot the noiseless mocks of a compact-DF recovery batch against the real data.

Top row: photometry, LOS dispersion, and radial/tangential PM dispersions for
the no-DM and cored-DM mocks, with the adopted errors, and the observed values
they replace (grey). Photometric zero points are profiled in the fit, so the
observed profile is shifted by its median offset from the no-DM mock.
Bottom row: (cored - no-DM)/adopted error, i.e. how separable the two
equilibria are bin by bin. Reads saved files only.

Usage: python bin/plot_compact_mocks.py results/df/compact_recovery_20260923_v2
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
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
KINDS = [("phot", r"$\mu$ [mag, arbitrary zero point]"), ("los", r"$\sigma_{\rm los}$ [km s$^{-1}$]"),
         ("pmr", r"$\sigma_{\rm pm,R}$ [mas yr$^{-1}$]"), ("pmt", r"$\sigma_{\rm pm,T}$ [mas yr$^{-1}$]")]
CASES = {"no_dm": ("tab:blue", "o", "no-DM mock"), "cored_dm": ("tab:orange", "s", "cored-DM mock")}


def read(path):
    return json.loads(Path(path).read_text())


def load(directory):
    """Return {kind: [(name, r_arcsec, value, err_lo, err_hi), ...]} for one snapshot directory."""
    data = {k: [] for k, _ in KINDS}
    for p in read(directory/"data_snapshot.json")["profiles"]:
        data[p["kind"]].append((p["name"], *(np.array(p[k], float) for k in ("r", "value", "err_lo", "err_hi"))))
    ph = read(directory/"photometry.json")
    s = np.array(ph["sigma_mag"], float)
    data["phot"].append(("photometry", np.array(ph["r_arcsec"], float), np.array(ph["mu"], float), s, s))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    out = parser.parse_args().batch.resolve()
    manifest = read(out/"batch.json")
    mocks = {c: load(out/"mocks"/c) for c in CASES}
    observed = load(out/"inputs")
    distance_kpc = read(out/"mocks/no_dm/mock.json")["config"]["distance_kpc"]
    pc_per_arcsec = distance_kpc*1e3*np.pi/(180*3600)

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), sharex="col", constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 1.4]))
    record = dict(batch=str(out), created_utc=datetime.now(timezone.utc).isoformat(),
                  source_commit=manifest["source_commit"], pc_per_arcsec=pc_per_arcsec, panels={})
    for col, (kind, ylabel) in enumerate(KINDS):
        ax, dx = axes[0, col], axes[1, col]
        panel = record["panels"][kind] = {}
        for i, (name, r, v, lo, hi) in enumerate(observed[kind]):
            if kind == "phot":
                v = v + np.median(mocks["no_dm"][kind][i][2] - v)
            ax.errorbar(r*pc_per_arcsec, v, yerr=[lo, hi], fmt=".", color="0.65", ms=3, lw=.6,
                        zorder=1, label="observed" if i == 0 else None)
            panel.setdefault("observed", {})[name] = dict(r_arcsec=r.tolist(), value=v.tolist())
        for case, (color, marker, label) in CASES.items():
            for i, (name, r, v, lo, hi) in enumerate(mocks[case][kind]):
                off = 1.03 if case == "cored_dm" else 1.0  # tiny radial offset so bars do not overlap
                ax.errorbar(r*pc_per_arcsec*off, v, yerr=[lo, hi], fmt=marker, color=color, ms=3, lw=.7,
                            mfc="none", zorder=2, label=label if i == 0 else None)
                panel.setdefault(case, {})[name] = dict(r_arcsec=r.tolist(), value=v.tolist(),
                                                       err_lo=lo.tolist(), err_hi=hi.tolist())
        for (name, r, v0, lo, hi), (_, _, v1, _, _) in zip(mocks["no_dm"][kind], mocks["cored_dm"][kind]):
            # Asymmetric errors: use the side the cored value lies on.
            err = np.where(v1 >= v0, hi, lo)
            z = (v1 - v0)/err
            if kind == "phot":
                z = z - np.average(z, weights=1/err**2)  # zero point is profiled in the fit
            dx.plot(r*pc_per_arcsec, z, ".-", ms=4, lw=.8, label=name)
            panel.setdefault("cored_minus_nodm_sigma", {})[name] = z.tolist()
        dx.axhline(0, color="black", lw=.8)
        ax.set(xscale="log", ylabel=ylabel, title={"phot": "surface brightness", "los": "line of sight",
                                                   "pmr": "PM radial", "pmt": "PM tangential"}[kind])
        if kind == "phot":
            ax.invert_yaxis()
        dx.set(xscale="log", xlabel="Projected radius [pc]", ylabel="(cored - no DM) / error")
        ax.legend(fontsize=7)
        dx.legend(fontsize=6)
    fig.suptitle(f"Noiseless compact DF mocks ({out.name}); errors are the adopted diagnostic errors")
    plot = ROOT/"plots"/(out.name+"_mocks.png")
    fig.savefig(plot, dpi=160)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/(out.name+"_mocks.json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
