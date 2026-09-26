import json, numpy as np, sys
sys.path.insert(0, "bin/nbody")
from run_nbody import cluster_centre, bound_mask, OCEN_TODAY
from analyse_tails import to_sky
out = {}
for m in ("A_nodm", "B_dm_phot"):
    ics = np.load(f"results/nbody/{m}/ics.npz"); sp = ics["species"]; mass = ics["mass"].astype(float)
    s = np.load(f"results/nbody/{m}/orbit/snap_today.npz"); pos, vel = s["pos"].astype(float), s["vel"].astype(float)
    lum = sp <= 1
    c, vc = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0))
    bound = bound_mask(pos, vel, mass, 0.3, c, vc)
    st = sp == 0; rem = sp == 1
    # sky coordinates relative to the model cluster
    l, b, d, vr, pml, pmb = to_sky(pos[st], vel[st]); lc, bc, dc, vrc, pmlc, pmbc = [x[0] for x in to_sky(c[None], vc[None])]
    dl = ((l-lc+180) % 360-180)*np.cos(np.radians(bc)); db = b-bc; sep = np.hypot(dl, db)
    ub = ~bound[st]; Nb = np.sum(~ub)
    row = dict(bound_star_mass=mass[st][~ub].sum(), unbound_star_mass=mass[st][ub].sum())
    row["tail_over_cluster_all"] = ub.sum()/Nb
    for R in (2, 5, 10, 20):
        k = ub & (sep < R)
        row[f"tail_over_cluster_within_{R}deg"] = k.sum()/Nb
        if k.sum() > 30:
            row[f"sigma_vlos_tail_{R}deg"] = float(np.std(vr[k]-vrc))
            row[f"width_deg_{R}deg"] = float(np.std(np.where(abs(dl[k]) > abs(db[k]), db[k], dl[k])))
    # recent mass-loss rate (last 500 Myr) per unit bound stellar mass
    rows = [json.loads(x) for x in open(f"results/nbody/{m}/orbit/diagnostics.jsonl")]
    t = np.array([r["t_myr"] for r in rows]); ms = np.array([r["bound_mass"]["stars"] for r in rows])
    i0, i1 = np.argmin(abs(t-1455.6)), np.argmin(abs(t-1955.6))
    row["mass_loss_rate_frac_per_Gyr_last500"] = (ms[i0]-ms[i1])/ms[i1]/0.5
    # remnants: fraction of remnants among unbound vs bound (luminous-to-dark)
    row["remnant_to_star_bound"] = mass[rem & bound].sum()/mass[st & bound].sum()
    row["remnant_to_star_tails"] = mass[rem & ~bound].sum()/mass[st & ~bound].sum()
    # bound cluster: projected LOS dispersion profile (line of sight = Sun direction) and outer density
    sun = np.array([-8178., 0, 0])  # Baumgardt frame: Sun at X = -R0? check: X points from Sun to GC -> Sun at X = -8.178 kpc
    los = (c-sun)/np.linalg.norm(c-sun)
    rel = pos-c; dv = vel-vc
    Rp = np.linalg.norm(rel-np.outer(rel@los, los), axis=1); vl = dv@los
    prof = {}
    for a, bb in ((1, 5), (5, 15), (15, 30), (30, 60), (60, 100), (100, 200)):
        k = st & (Rp >= a) & (Rp < bb) & (np.abs(rel@los) < 500)
        prof[f"{a}-{bb}pc"] = (int(k.sum()), float(np.std(vl[k])) if k.sum() > 30 else None, float(np.mean(bound[k])) if k.sum() else None)
    row["los_profile"] = prof
    rh = np.sort(np.linalg.norm(rel[st & bound], axis=1)); row["r_half_bound_stars"] = rh[len(rh)//2]
    out[m] = row
for k in out["A_nodm"]:
    if k == "los_profile": continue
    a, b = out["A_nodm"][k], out["B_dm_phot"].get(k)
    print(f"{k:38s} A {a:12.4g}   B {b:12.4g}   A/B {a/b if b else float('nan'):.2f}")
print("projected LOS dispersion of stars (N, sigma km/s, bound fraction), |depth| < 500 pc:")
for r in out["A_nodm"]["los_profile"]:
    a, b = out["A_nodm"]["los_profile"][r], out["B_dm_phot"]["los_profile"][r]
    print(f"   R {r:9s} A {a}   B {b}")
json.dump(out, open("results/plot_data/nbody_discriminants_today.json", "w"), default=float, indent=1)
