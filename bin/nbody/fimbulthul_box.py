#!/usr/bin/env python3
"""Anatomy of the model debris in the Fimbulthul box (l -52..-32, b 30..40 deg).

For the unbound star particles in the box: heliocentric distance, v_helio, proper motions,
energy offset from the cluster (dE < 0 leading, > 0 trailing), time offset along the
cluster orbit (6D match, as in analyse_tails), and the time each particle became unbound
(first snapshot in which it is unbound, from the 10-Myr snapshots, using the per-snapshot
cluster centre from diagnostics.jsonl). Also applies an I19-like selection: proper
motions within 1 mas/yr of the I19 STREAMFINDER members' trend (their Table 1, first
10 rows: mu_alpha* -9.5..-18.4, mu_delta -9.7..-11.5 mas/yr) and parallax-distance 3-5 kpc.

Usage: python bin/nbody/fimbulthul_box.py
"""
import glob, json, os, sys
from pathlib import Path
for k in ("OMP_NUM_THREADS",): os.environ.setdefault(k, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import Galactocentric, CartesianDifferential, SkyCoord, ICRS
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"bin"/"nbody"))
from run_nbody import cluster_centre, bound_mask, mw_potential, AGAMA_T_MYR
from analyse_tails import R_SUN_KPC, Z_SUN_KPC, V_SUN, V_SCALE
from scipy.spatial import cKDTree

FR = Galactocentric(galcen_distance=R_SUN_KPC*u.kpc, z_sun=Z_SUN_KPC*u.kpc, galcen_v_sun=CartesianDifferential(*[v*u.km/u.s for v in V_SUN]))
I19 = dict(ra=np.array([200.478245, 204.456441, 205.286684, 200.220408, 210.247730, 211.095707, 213.725755, 213.443432, 212.585609, 206.290828]),
           pmra=np.array([-10.610, -14.119, -15.146, -9.488, -16.371, -16.720, -16.934, -18.403, -16.477, -15.188]),
           pmdec=np.array([-9.710, -10.936, -11.223, -9.677, -10.787, -11.470, -11.491, -10.558, -11.187, -11.159]))


def sky(pos, vel):
    c = SkyCoord(x=-pos[:, 0]*u.pc, y=pos[:, 1]*u.pc, z=pos[:, 2]*u.pc, v_x=-vel[:, 0]*u.km/u.s, v_y=vel[:, 1]*u.km/u.s, v_z=vel[:, 2]*u.km/u.s, frame=FR)
    g = c.transform_to("galactic"); e = c.transform_to(ICRS())
    l = g.l.deg; return (np.where(l > 180, l-360, l), g.b.deg, g.distance.kpc, g.radial_velocity.to_value(u.km/u.s),
                        e.ra.deg, e.pm_ra_cosdec.to_value(u.mas/u.yr), e.pm_dec.to_value(u.mas/u.yr))


def main():
    agama, pot = mw_potential("McMillan17")
    fig, ax = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True); ax = ax.ravel()
    rec = {}
    for m, lab, col in (("A_nodm", "A: no DM", "tab:blue"), ("B_dm_phot", "B: DM halo", "tab:orange")):
        ics = np.load(ROOT/f"results/nbody/{m}/ics.npz"); sp = ics["species"]; mass = ics["mass"].astype(float)
        s = np.load(ROOT/f"results/nbody/{m}/orbit/snap_today.npz"); pos, vel = s["pos"].astype(float), s["vel"].astype(float)
        lum = sp <= 1; c, vc = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0)); bound = bound_mask(pos, vel, mass, 0.3, c, vc)
        idx = np.where((sp == 0) & ~bound)[0]
        l, b, d, vr, ra, pmra, pmdec = sky(pos[idx], vel[idx])
        box = (l > -52) & (l < -32) & (b > 30) & (b < 40)
        I = idx[box]
        # energy offset and orbit time offset
        E = 0.5*np.sum(vel[I]**2, 1)+pot.potential(pos[I]*1e-3); Ec = 0.5*vc@vc+pot.potential((c*1e-3)[None])[0]
        tr = []
        for sgn in (-1, 1):
            t, o = agama.orbit(potential=pot, ic=np.concatenate((c*1e-3, vc)), time=sgn*1500/AGAMA_T_MYR, trajsize=6001); o = np.asarray(o)
            tr.append((np.asarray(t)*AGAMA_T_MYR, o[:, :3]*1e3, o[:, 3:]))
        tt = np.concatenate((tr[0][0][::-1], tr[1][0][1:])); pp = np.vstack((tr[0][1][::-1], tr[1][1][1:])); vv = np.vstack((tr[0][2][::-1], tr[1][2][1:]))
        _, j = cKDTree(np.hstack((pp, V_SCALE*vv))).query(np.hstack((pos[I], V_SCALE*vel[I]))); dt = tt[j]
        # stripping time: first 10-Myr snapshot in which the particle is unbound
        diag = {round(r["t_myr"], 1): r for r in (json.loads(x) for x in open(ROOT/f"results/nbody/{m}/orbit/diagnostics.jsonl"))}
        t_esc = np.full(len(I), np.nan); pending = np.ones(len(I), bool)
        for f in sorted(glob.glob(str(ROOT/f"results/nbody/{m}/orbit/snap_0*.npz"))):
            z = np.load(f); T = float(z["t_myr"])
            if T > 1956 or not pending.any(): continue
            r = diag[round(T, 1)]; cc = np.array(r["centre_pc"]); vcc = np.array(r["vcentre_kms"])
            P = z["pos"][I].astype(float); V = z["vel"][I].astype(float)
            # energy in a Plummer approximation of the bound cluster (M_bound, a = 1.3 r_half): cheap per-snapshot test
            Mb = sum(r["bound_mass"].values()); a = 1.3*r["r_half_stars_pc"]
            R = np.linalg.norm(P-cc, axis=1); e = 0.5*np.sum((V-vcc)**2, 1)-4.30091727e-3*Mb/np.sqrt(R**2+a**2)
            new = pending & (e > 0); t_esc[new] = T; pending &= ~new
        dE = E-Ec
        sel_pm = np.zeros(box.sum(), bool)
        for k2, (raa, pa, pdd) in enumerate(zip(ra[box], pmra[box], pmdec[box])):
            ref_a = np.interp(raa, np.sort(I19["ra"]), I19["pmra"][np.argsort(I19["ra"])]); ref_d = np.interp(raa, np.sort(I19["ra"]), I19["pmdec"][np.argsort(I19["ra"])])
            sel_pm[k2] = np.hypot(pa-ref_a, pdd-ref_d) < 1.0
        sel = sel_pm & (d[box] > 3) & (d[box] < 5)
        rec[m] = dict(n_box=int(box.sum()), n_I19like=int(sel.sum()),
                      v_I19like=(float(np.median(vr[box][sel])), float(np.std(vr[box][sel]))) if sel.sum() > 2 else None,
                      leading_frac=float(np.mean(dE < 0)), dt_pct=np.percentile(dt, [5, 50, 95]).tolist(), t_esc_pct=np.nanpercentile(t_esc, [5, 50, 95]).tolist(),
                      cluster_lb=[float(x[0]) for x in sky(c[None], vc[None])[:2]])
        ax[0].scatter(d[box], vr[box], s=12, color=col, alpha=.7, label=f"{lab} ({box.sum()})")
        ax[0].scatter(d[box][sel], vr[box][sel], s=40, facecolor="none", edgecolor="black")
        sc = ax[1 if m == "A_nodm" else 2].scatter(dt, vr[box], c=t_esc, cmap="viridis", vmin=0, vmax=1956, s=14)
        ax[1 if m == "A_nodm" else 2].set(xlabel="time offset along the orbit [Myr] (>0 leading)", ylabel=r"$v_{\rm helio}$ [km/s]",
                                           title=f"{lab}: box debris by orbit phase (colour: escape time)")
        ax[3].scatter(dE, vr[box], s=12, color=col, alpha=.7, label=lab)
        ax[4].scatter(pmra[box], pmdec[box], s=12, color=col, alpha=.6, label=lab)
        ax[5].hist(t_esc, np.linspace(0, 1960, 50), histtype="step", color=col, lw=2, label=lab)
    fig.colorbar(sc, ax=ax[2], label="escape time t [Myr] (present = 1956)")
    for a in (ax[0], ax[3]):
        a.axhspan(199.7-5.4, 199.7+5.4, color="0.85", zorder=0)
    ax[0].set(xlabel="heliocentric distance [kpc]", ylabel=r"$v_{\rm helio}$ [km/s]", title="(a) box debris; circles: I19-like PM + distance 3-5 kpc; band: I19")
    ax[0].axvline(4.1, color="black", ls=":"); ax[0].legend(fontsize=8)
    ax[3].axvline(0, color="black", lw=.7); ax[3].set(xlabel=r"$\Delta E$ from the cluster [km$^2$ s$^{-2}$] (<0 leading)", ylabel=r"$v_{\rm helio}$ [km/s]", title="(d) energy offset"); ax[3].legend(fontsize=8)
    ax[4].scatter(I19["pmra"], I19["pmdec"], marker="*", s=90, color="black", label="I19 members (Table 1, 10 rows)")
    ax[4].set(xlabel=r"$\mu_{\alpha*}$ [mas/yr]", ylabel=r"$\mu_\delta$ [mas/yr]", title="(e) proper motions of box debris"); ax[4].legend(fontsize=8)
    ax[5].set(xlabel="escape time [Myr]", ylabel="particles", title="(f) when the box debris left the cluster"); ax[5].legend(fontsize=8)
    fig.suptitle("Which debris lands in the Fimbulthul box, and why its velocities differ between models")
    fig.savefig(ROOT/"plots/nbody_fimbulthul_box_today.png", dpi=130)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
