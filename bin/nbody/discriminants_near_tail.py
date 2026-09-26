import numpy as np, sys
sys.path.insert(0, "bin/nbody")
from run_nbody import cluster_centre, bound_mask
from analyse_tails import to_sky
for m in ("A_nodm", "B_dm_phot"):
    ics = np.load(f"results/nbody/{m}/ics.npz"); sp = ics["species"]; mass = ics["mass"].astype(float)
    s = np.load(f"results/nbody/{m}/orbit/snap_today.npz"); pos, vel = s["pos"].astype(float), s["vel"].astype(float)
    lum = sp <= 1; c, vc = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0)); bound = bound_mask(pos, vel, mass, 0.3, c, vc)
    st = sp == 0
    l, b, d, vr, *_ = to_sky(pos[st], vel[st]); lc, bc, dc, vrc = [x[0] for x in to_sky(c[None], vc[None])[:4]]
    dl = ((l-lc+180) % 360-180)*np.cos(np.radians(bc)); db = b-bc; sep = np.hypot(dl, db); ub = ~bound[st]
    for R in (2, 5, 10):
        k = ub & (sep < R); dv = vr[k]-vrc; dd = d[k]-dc
        mad = 1.4826*np.median(np.abs(dv-np.median(dv)))
        near = k & (np.abs(d-dc) < 0.5)                       # same wrap: within 0.5 kpc in distance
        dvn = vr[near]-vrc
        print(f"{m} <{R:2d} deg: N {k.sum():5d}  sigma_vlos std {np.std(dv):5.1f}  robust {mad:5.1f} km/s | same-wrap (|dd|<0.5 kpc) N {near.sum():5d} sigma {np.std(dvn):5.1f} robust {1.4826*np.median(np.abs(dvn-np.median(dvn))):5.1f} ; N_tail/N_bound {near.sum()/np.sum(~ub):.2e}")
