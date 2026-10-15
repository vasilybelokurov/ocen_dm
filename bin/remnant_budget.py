#!/usr/bin/env python3
"""Present-day remnant mass budget of an old, metal-poor cluster from the IMF, against our fitted remnant components.

Assumptions (all stated, the uncertain ones bracketed):
  IMF      Kroupa (2001): dN/dm ~ m^-1.3 (0.08-0.5 Msun), m^-2.3 (0.5-100 Msun)
  turnoff  0.79 Msun (our MIST fit at 12.5 Gyr gives 0.779-0.789)
  WDs      progenitors 0.79-8 Msun, IFMR m_WD = 0.109 m_i + 0.394 (Kalirai et al. 2008)
  NSs      progenitors 8-20 Msun, 1.3 Msun each; retention r_NS = 10% (default) or 30%
  BHs      progenitors > 20 Msun, m_BH = f m_i capped at 40 Msun, f = 0.25 or 0.4 (low Z); retention r_BH 0-100%
  escape   none: no preferential loss of low-mass stars (an evolved cluster that lost stars has HIGHER remnant fractions)
"Heavier than the tracer" = remnants that can be more concentrated than the ~0.8 Msun stars that trace the light:
BHs + NSs + WDs heavier than 0.8 Msun.

Outputs: plots/remnant_budget.png, results/plot_data/remnant_budget.json.
Usage: python bin/remnant_budget.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
M_TO, M_UP = 0.79, 100.
M_WD_TRACER = 0.8                       # WDs heavier than this can be more concentrated than the tracers
MODELS = {"DF base model": dict(frac=2.73e5/3.26e6, q=2.73e5/3.0e6, color="C0"),
          "Jeans model (remnants + black hole)": dict(frac=(5.36e5+4.25e4)/3.24e6, q=(5.36e5+4.25e4)/2.66e6, color="C1")}
PUBLISHED_BH = {"Zocchi+2019": 0.045, "Baumgardt+2019": 0.046, "Dickson+2024": 0.059, "Dickson+2026a": 0.083}


def imf(m):
    return np.where(m < 0.5, (m/0.5)**-1.3, (m/0.5)**-2.3)


def integ(f, a, b):
    return quad(f, a, b, limit=200)[0]


def budget(r_ns, r_bh, f_bh):
    stars = integ(lambda m: m*imf(m), 0.08, M_TO)
    wd = integ(lambda m: (0.109*m+0.394)*imf(m), M_TO, 8.)
    m_cut = (M_WD_TRACER-0.394)/0.109
    wd_heavy = integ(lambda m: (0.109*m+0.394)*imf(m), m_cut, 8.)
    ns = r_ns*1.3*integ(imf, 8., 20.)
    bh = r_bh*integ(lambda m: min(f_bh*m, 40.)*imf(m), 20., M_UP)
    tot = stars+wd+ns+bh
    heavy = wd_heavy+ns+bh
    return dict(stars=stars/tot, wd=wd/tot, ns=ns/tot, bh=bh/tot, all_remnants=(wd+ns+bh)/tot,
                heavier_than_tracer=heavy/tot, heavier_over_light_following=heavy/(tot-heavy),
                mean_wd_mass=wd/integ(imf, M_TO, 8.), wd_heavy=wd_heavy/tot)


def main():
    rows = []
    for f_bh in (0.25, 0.4):
        for r_ns in (0.1, 0.3):
            for r_bh in (0., 0.5, 1.):
                rows.append(dict(f_bh=f_bh, r_ns=r_ns, r_bh=r_bh, **budget(r_ns, r_bh, f_bh)))
    r = np.linspace(0, 1, 41)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
    for f_bh, ls in ((0.25, "-"), (0.4, "--")):
        b = [budget(0.1, x, f_bh) for x in r]
        ax[0].plot(r, [x["all_remnants"] for x in b], color="k", ls=ls, label=f"all remnants (m_BH = {f_bh} m_i)")
        ax[0].plot(r, [x["wd"] for x in b], color="0.5", ls=ls, label=f"white dwarfs (m_BH = {f_bh} m_i)")
        ax[0].plot(r, [x["bh"] for x in b], color="C3", ls=ls, label=f"black holes (m_BH = {f_bh} m_i)")
        ax[1].plot(r, [x["heavier_than_tracer"] for x in b], color="C3", ls=ls,
                   label=f"IMF budget heavier than the tracers (m_BH = {f_bh} m_i)")
    for name, mdl in MODELS.items():
        ax[1].axhline(mdl["frac"], color=mdl["color"], lw=2, label=f"{name}: {mdl['frac']:.1%}")
    label_y = {"Zocchi+2019": 0.038, "Baumgardt+2019": 0.050, "Dickson+2024": 0.062, "Dickson+2026a": 0.086}
    for name, v in PUBLISHED_BH.items():
        ax[1].plot(1.02, v, marker="<", color="0.3", ms=7, clip_on=False)
        ax[1].text(1.04, label_y[name], f"{name} ({v:.1%})", fontsize=7, va="center", transform=ax[1].get_yaxis_transform())
    ax[0].set(xlabel="black-hole retention fraction", ylabel="fraction of present total mass", ylim=(0, 0.4),
              title="Kroupa IMF at 12.5 Gyr, no stellar escape (NS retention 10%)")
    ax[1].set(xlabel="black-hole retention fraction", ylabel="fraction of present total mass", ylim=(0, 0.3),
              title="Mass that can be more concentrated than the ~0.8 Msun tracers")
    ax[0].legend(fontsize=8)
    ax[1].legend(fontsize=8, loc="upper left")
    ax[1].text(0.98, 0.06, "published omega Cen\nblack-hole fractions", fontsize=7, ha="right", transform=ax[1].transAxes, color="0.3")
    fig.savefig(ROOT/"plots/remnant_budget.png", dpi=140)
    (ROOT/"results/plot_data").mkdir(parents=True, exist_ok=True)
    (ROOT/"results/plot_data/remnant_budget.json").write_text(json.dumps(dict(rows=rows, models=MODELS, published_bh=PUBLISHED_BH), indent=1))
    for x in rows:
        print(f"f_BH={x['f_bh']} r_NS={x['r_ns']:.0%} r_BH={x['r_bh']:.0%}: stars {x['stars']:.1%} WD {x['wd']:.1%} NS {x['ns']:.2%} "
              f"BH {x['bh']:.1%} all {x['all_remnants']:.1%} | heavier-than-tracer {x['heavier_than_tracer']:.1%} "
              f"(/light-following {x['heavier_over_light_following']:.3f}; heavy WDs {x['wd_heavy']:.1%})")
    print("mean WD mass", round(rows[0]["mean_wd_mass"], 3))


if __name__ == "__main__":
    main()
