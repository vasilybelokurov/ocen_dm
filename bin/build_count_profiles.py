#!/usr/bin/env python3
"""Build the two Poisson star-count profiles used as the tracer-density constraint.

1. HST (oMEGACat): all catalogue stars with a measured F625W < 19 (no quality
   flag), in log annuli from 2 to 250 arcsec about the pixel centre
   (15000, 15000; 0.04 arcsec/px). The authors' `selection_hq_f625w` flag is
   NOT applied: the fraction of F625W<19 stars carrying it rises from 0.81-0.83
   at 8-22 arcsec to 0.94-0.95 beyond 175 arcsec (crowding-dependent PSF-quality
   cuts), which would imprint a 12-15% radial incompleteness on the profile.
   Without the flag the ratio N(F625W<20)/N(F625W<19) is 2.0-2.1 in every
   annulus (2026-09-23), so the bright sample is complete. Footprint coverage
   per annulus comes from the occupancy of 3-arcsec cells by any catalogue star
   (>= 0.99 inside 250 arcsec; set to 1 inside 30 arcsec). Bins stop at 250
   arcsec: beyond it the mosaic edge lowers coverage and F814W availability.
   Field contamination is < 0.5% at 250 arcsec, so no field term is fitted.
2. Gaia EDR3 (Vasiliev & Baumgardt 2021 catalogue of the field): all stars with
   G < 17 from 300 to 2400 arcsec, members and field together, full annulus
   areas. Completeness check: N(G<16)/N(G<15) and N(G<17)/N(G<16) are flat
   (1.8-2.0) from 300 arcsec out, whereas G<18 is depleted inside ~500 arcsec.
   The uniform field density is fitted (`fit_field=True`); membership
   probabilities are not used.

Radial nodes are 8 equal-count quantiles of the counted stars' radii per bin,
so the model is averaged over the stars actually counted. Products go to
data/processed/kinematics/ as ECSV (table) and JSON (fit input).

Usage: python bin/build_count_profiles.py
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.paths import processed_dir, raw_dir

N_NODES = 8


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def nodes(r):
    q = (np.arange(N_NODES)+.5)/N_NODES
    return np.quantile(np.asarray(r, float), q)


def hst_counts(edges, mag_cut=19., cell=3., full_inside=30.):
    path = raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits"
    t = fits.open(path)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.)
    y = 0.04*(np.asarray(t["y"], float)-15000.)
    ok = np.isfinite(x) & np.isfinite(y)
    R = np.hypot(x, y)
    f = np.asarray(t["f625w"], float)
    sel = ok & np.isfinite(f) & (f < mag_cut)
    # footprint coverage from cell occupancy by any star
    ix = np.floor(x[ok]/cell).astype(np.int64)+10**5
    iy = np.floor(y[ok]/cell).astype(np.int64)+10**5
    cells = np.unique(ix*10**6+iy)
    cx = (cells//10**6-10**5+.5)*cell
    cy = (cells % 10**6-10**5+.5)*cell
    Rc = np.hypot(cx, cy)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        annulus = np.pi*(hi**2-lo**2)
        cov = 1. if hi <= full_inside else min(np.sum((Rc >= lo) & (Rc < hi))*cell**2/annulus, 1.)
        k = sel & (R >= lo) & (R < hi)
        rows.append((lo, hi, int(k.sum()), cov, annulus*cov, nodes(R[k])))
    meta = dict(catalogue=str(path.relative_to(ROOT)), catalogue_sha256=sha256(path), centre="pixel (15000, 15000), 0.04 arcsec/px",
                selection=f"all stars with measured F625W < {mag_cut:g} (no quality flag)", n_selected=int(sel.sum()),
                coverage_method=f"occupancy of {cell:g}-arcsec cells by any catalogue star; 1 inside {full_inside:g} arcsec")
    return rows, meta


def gaia_counts(edges, mag_cut=17.):
    path = processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv"
    t = Table.read(path)
    x = (np.asarray(t["ra"], float)-OCEN_RA)*np.cos(np.radians(OCEN_DEC))*3600.
    y = (np.asarray(t["dec"], float)-OCEN_DEC)*3600.
    R = np.hypot(x, y)
    g = np.asarray(t["g_mag"], float)
    sel = np.isfinite(g) & (g < mag_cut)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = sel & (R >= lo) & (R < hi)
        rows.append((lo, hi, int(k.sum()), 1., np.pi*(hi**2-lo**2), nodes(R[k])))
    meta = dict(catalogue=str(path.relative_to(ROOT)), catalogue_sha256=sha256(path),
                centre=f"({OCEN_RA}, {OCEN_DEC}) deg", selection=f"all catalogue stars with G < {mag_cut:g} (members and field)",
                n_selected=int(sel.sum()), coverage_method="full annuli (Gaia all-sky)")
    return rows, meta


def build(name, rows, meta, source, fit_field, note):
    lo, hi, n, cov, area, nd = (np.array([r[i] for r in rows]) for i in range(6))
    profile = CountProfile(name, lo, hi, nd, n.astype(int), area, source, meta["selection"], fit_field, note,
                           dict(meta, built_utc=datetime.now(timezone.utc).isoformat(), coverage=cov.tolist()))
    out = processed_dir()/"kinematics"
    table = Table(dict(r_lower=lo, r_upper=hi, r_median=profile.r_median, n_stars=n.astype(int), coverage=cov,
                       area_arcsec2=area, density_per_arcmin2=n/area*3600., poisson_err_per_arcmin2=np.sqrt(np.maximum(n, 1))/area*3600.),
                  meta=dict(meta, product=name, fit_field=fit_field, note=note))
    table.write(out/f"{name}.ecsv", format="ascii.ecsv", overwrite=True)
    (out/f"{name}.json").write_text(json.dumps(profile.to_dict()))
    print(f"{name}: {len(rows)} bins, {int(n.sum())} stars; written to {out}")
    for r in rows:
        print("   %7.1f-%7.1f  N=%6d  cov %.3f" % (r[0], r[1], r[2], r[3]))
    return profile


def main():
    hst_edges = np.geomspace(2., 250., 21)
    build("ocen_counts_hst_f625w19", *hst_counts(hst_edges),
          source="oMEGACat catalogue (Haeberle et al. 2025)", fit_field=False,
          note="All F625W<19 detections, no quality flag (the flag is 12-15% radially incomplete); N(<20)/N(<19) = 2.0-2.1 at all radii; field < 0.5% at 250 arcsec, no field term")
    gaia_edges = np.array([300., 350., 400., 450., 500., 600., 700., 800., 1000., 1300., 1600., 2000., 2400.])
    build("ocen_counts_gaia_g17", *gaia_counts(gaia_edges),
          source="Vasiliev & Baumgardt 2021 EDR3 catalogue of the omega Cen field", fit_field=True,
          note="G<17 complete from 300 arcsec (N(<16)/N(<15), N(<17)/N(<16) flat); uniform field density fitted")


if __name__ == "__main__":
    main()
