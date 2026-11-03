"""Marginalising the DF likelihood over the cluster centre (docs/CENTRE_MARGINALISATION.md).

The centre c = (dRA*, dDec) [arcsec, east and north positive] is a nuisance parameter with an
isotropic Gaussian prior about the reporting origin (AvdM10 = oMEGACat pixel 15000, 15000 =
(OCEN_RA, OCEN_DEC) to 0.002"). For each centre on a Gauss-Hermite grid the data that depend on
the centre are rebuilt with the reference bin edges and selections, and

    ln L_m(theta) = logsumexp_k [ ln w_k + ln L(theta | d(c_k)) ],

with ln L the full log-likelihood of DFJointProblem.evaluate ("loglike_profiled": split-normal
kinematics + Poisson counts with the Gamma term, count amplitudes and field profiled per centre).
The deviance objective is never summed across centres: its saturated term changes with the data.

Rebuilt per centre (all other bins are shared):
* HST F625W < 19 counts, every bin, coverage about c (same cell-occupancy method as
  bin/build_count_profiles.py; the cell grid stays fixed on the sky);
* Gaia G < 17 counts, every bin (full annuli);
* hst_pm_{radial,tangential}_ours, bins with r_upper <= R_KIN (our mixture estimator, nodes
  re-measured);
* muse_los_dispersion, bins with r_upper <= R_KIN: published value + [ours(c) - ours(AvdM10)],
  with ours = ML Gaussian with mean v0 + a sin(phi) + b cos(phi) and per-star errors on
  selection_hq_pm_and_los stars (this selection reproduces the published n_stars of all 29 bins
  about AvdM10). Published errors and the published rotation term are kept.
"""
from __future__ import annotations

from dataclasses import replace
from functools import lru_cache

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from .cluster import OCEN_DEC, OCEN_RA
from .kinematics.counts import CountProfile
from .kinematics.likelihood import KinematicData, _hst_ours
from .paths import processed_dir, raw_dir

__all__ = ["R_KIN", "centre_grid", "hst_counts", "gaia_counts", "hst_pm_table", "muse_sigma",
           "CentreData", "marginal_loglike"]

#: kinematic bins with r_upper <= R_KIN (arcsec) are re-centred
R_KIN = 15.0
N_NODES = 8                      # sub-rings per count bin (as in bin/build_count_profiles.py)
HST_NAME, GAIA_NAME = "ocen_counts_hst_f625w19", "ocen_counts_gaia_g17"


def centre_grid(sigma: float, n: int):
    """Product Gauss-Hermite nodes (K, 2) in arcsec and weights (K,) summing to 1 for N(0, sigma^2 I)."""
    if sigma == 0:
        return np.zeros((1, 2)), np.ones(1)
    x, w = np.polynomial.hermite_e.hermegauss(n)          # weight exp(-x^2/2)
    w = w/w.sum()
    X, Y = np.meshgrid(x*sigma, x*sigma, indexing="ij")
    W = np.outer(w, w)
    return np.column_stack([X.ravel(), Y.ravel()]), W.ravel()


# ------------------------------------------------------------------ counts ---
def _quadrature(lo, hi, covered_area=None):
    """Copy of bin/build_count_profiles.quadrature (pinned by tests/test_centre_marg.py)."""
    e = np.geomspace(lo, hi, N_NODES+1)
    ring = np.pi*(e[1:]**2-e[:-1]**2)
    a = ring if covered_area is None else np.array([covered_area(r1, r2) for r1, r2 in zip(e[:-1], e[1:])])
    if a.sum() <= 0:
        raise ValueError(f"no covered area in [{lo}, {hi})")
    return np.sqrt(.5*(e[1:]**2+e[:-1]**2)), a/a.sum(), float(a.sum())


@lru_cache(maxsize=1)
def _hst_catalogue(cell=3., sub=10, r_keep=300.):
    """Star positions (arcsec about pixel 15000,15000), F625W, and the occupied-cell point set."""
    with fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits", memmap=True) as h:
        t = h[1].data
        x = -0.04*(np.asarray(t["x"], float)-15000.)
        y = 0.04*(np.asarray(t["y"], float)-15000.)
        f = np.asarray(t["f625w"], float)
        v = np.asarray(t["vlos"], float)
        ev = np.asarray(t["vlos_err"], float)
        hq_los = np.asarray(t["selection_hq_pm_and_los"]) > 0
    ok = np.isfinite(x) & np.isfinite(y)
    ix = np.floor(x[ok]/cell).astype(np.int64)+10**5
    iy = np.floor(y[ok]/cell).astype(np.int64)+10**5
    cells = np.unique(ix*10**6+iy)
    cx = (cells//10**6-10**5+.5)*cell
    cy = (cells % 10**6-10**5+.5)*cell
    off = ((np.arange(sub)+.5)/sub-.5)*cell
    px = (cx[:, None, None]+off[None, :, None]+0*off[None, None, :]).ravel()
    py = (cy[:, None, None]+0*off[None, :, None]+off[None, None, :]).ravel()
    near = np.hypot(px, py) < r_keep
    sel = ok & np.isfinite(f) & (f < 19.)
    los = ok & hq_los & np.isfinite(v) & np.isfinite(ev)
    return dict(x=x[sel], y=y[sel], px=px[near], py=py[near], dA=(cell/sub)**2,
                lx=x[los], ly=y[los], v=v[los], ev=ev[los])


def hst_counts(template: CountProfile, offset=(0., 0.), full_inside=30.) -> CountProfile:
    """HST F625W < 19 counts in ``template``'s bins about the shifted centre (coverage about it too)."""
    c = _hst_catalogue()
    dx, dy = offset
    R = np.hypot(c["x"]-dx, c["y"]-dy)
    Rp = np.sort(np.hypot(c["px"]-dx, c["py"]-dy))
    if template.r_upper.max() > 299.:
        raise ValueError("coverage point set is truncated at 300 arcsec")

    def covered(r1, r2):
        if r2 <= full_inside:
            return np.pi*(r2**2-r1**2)
        return min((np.searchsorted(Rp, r2)-np.searchsorted(Rp, r1))*c["dA"], np.pi*(r2**2-r1**2))

    nodes, weights, area, n = [], [], [], []
    for lo, hi in zip(template.r_lower, template.r_upper):
        nd, wt, a = _quadrature(lo, hi, covered)
        nodes.append(nd); weights.append(wt); area.append(a)
        n.append(int(np.count_nonzero((R >= lo) & (R < hi))))
    return replace(template, r_nodes=np.array(nodes), node_weights=np.array(weights),
                   area_arcsec2=np.array(area), counts=np.array(n, int),
                   meta=dict(template.meta, centre_offset_arcsec=[float(dx), float(dy)]))


@lru_cache(maxsize=1)
def _gaia_catalogue():
    t = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    x = (np.asarray(t["ra"], float)-OCEN_RA)*np.cos(np.radians(OCEN_DEC))*3600.
    y = (np.asarray(t["dec"], float)-OCEN_DEC)*3600.
    g = np.asarray(t["g_mag"], float)
    sel = np.isfinite(g) & (g < 17.)
    return x[sel], y[sel]


def gaia_counts(template: CountProfile, offset=(0., 0.)) -> CountProfile:
    """Gaia G < 17 counts in ``template``'s bins about the shifted centre (full annuli: geometry unchanged)."""
    x, y = _gaia_catalogue()
    R = np.hypot(x-offset[0], y-offset[1])
    n = [int(np.count_nonzero((R >= lo) & (R < hi))) for lo, hi in zip(template.r_lower, template.r_upper)]
    return replace(template, counts=np.array(n, int),
                   meta=dict(template.meta, centre_offset_arcsec=[float(offset[0]), float(offset[1])]))


# -------------------------------------------------------------- kinematics ---
def hst_pm_table(edges, offset=(0., 0.)) -> Table:
    """Our HST PM dispersion product (flagged stars, mixture estimator) in ``edges`` about the shifted centre."""
    from .kinematics.hst_profile import hst_profile, load_hst_sample
    s = load_hst_sample(require_flag=True, offset=offset)
    return hst_profile(edges_arcsec=edges, require_flag=True, min_stars=0, sample=s)


def _ml_sigma(v, ev, phi):
    """ML dispersion for v ~ N(v0 + a sin phi + b cos phi, sigma^2 + ev^2)."""
    X = np.column_stack([np.ones_like(v), np.sin(phi), np.cos(phi)])

    def nll(lns):
        var = np.exp(2*lns)+ev**2
        W = 1/var
        beta = np.linalg.solve(X.T@(X*W[:, None]), X.T@(W*v))
        r = v-X@beta
        return 0.5*np.sum(np.log(var)+r*r*W)

    s0 = np.log(max(np.std(v), 1.))
    res = minimize_scalar(nll, bounds=(s0-4, s0+2), method="bounded", options=dict(xatol=1e-10))
    return float(np.exp(res.x))


def muse_sigma(r_lower, r_upper, offset=(0., 0.)) -> np.ndarray:
    """Our MUSE LOS dispersion estimate per bin about the shifted centre (km/s)."""
    c = _hst_catalogue()
    x, y = c["lx"]-offset[0], c["ly"]-offset[1]
    R, phi = np.hypot(x, y), np.arctan2(x, y)
    out = []
    for lo, hi in zip(r_lower, r_upper):
        m = (R >= lo) & (R < hi)
        out.append(_ml_sigma(c["v"][m], c["ev"][m], phi[m]))
    return np.array(out)


# ------------------------------------------------------------------ data ---
class CentreData:
    """Builds the per-centre KinematicData and count products from a frozen run's inputs."""

    def __init__(self, data: KinematicData, counts, r_kin: float = R_KIN):
        self.data, self.counts, self.r_kin = data, tuple(counts), float(r_kin)
        prof = {p.name: p for p in data.profiles}
        self.hst = {k: prof[k] for k in ("hst_pm_radial_ours", "hst_pm_tangential_ours") if k in prof}
        if self.hst:
            p = next(iter(self.hst.values()))
            self.n_hst = int(np.count_nonzero(p.r_upper <= self.r_kin+1e-9))
            self.hst_edges = np.r_[p.r_lower[:self.n_hst], p.r_upper[self.n_hst-1]]
        self.muse = prof.get("muse_los_dispersion")
        if self.muse is not None:
            self.n_muse = int(np.count_nonzero(self.muse.r_upper <= self.r_kin+1e-9))
            self._muse0 = muse_sigma(self.muse.r_lower[:self.n_muse], self.muse.r_upper[:self.n_muse])

    def _splice(self, p, new, n):
        """Replace the first n bins of profile p by those of profile new."""
        cat = lambda a, b: None if a is None else np.concatenate([np.asarray(b)[:n], np.asarray(a)[n:]])
        return replace(p, r=cat(p.r, new.r), value=cat(p.value, new.value), err_lo=cat(p.err_lo, new.err_lo),
                       err_hi=cat(p.err_hi, new.err_hi),
                       streaming2=cat(p.streaming2, new.streaming2) if p.streaming2 is not None else None,
                       r_nodes=None if p.r_nodes is None else np.vstack([new.r_nodes[:n], p.r_nodes[n:]]))

    def kinematics(self, offset) -> KinematicData:
        out = []
        tab = hst_pm_table(self.hst_edges, offset) if self.hst else None
        if tab is not None and len(tab) != self.n_hst:
            raise ValueError("HST PM rebuild returned a different number of bins")
        for p in self.data.profiles:
            if p.name in self.hst:
                out.append(self._splice(p, _hst_ours(p.kind, table=tab), self.n_hst))
            elif self.muse is not None and p.name == self.muse.name:
                d = muse_sigma(p.r_lower[:self.n_muse], p.r_upper[:self.n_muse], offset)-self._muse0
                value = np.asarray(p.value, float).copy()
                value[:self.n_muse] += d
                out.append(replace(p, value=value))
            else:
                out.append(p)
        return KinematicData(tuple(out))

    def count_products(self, offset):
        out = []
        for c in self.counts:
            if c.name.startswith(HST_NAME):
                out.append(hst_counts(c, offset))
            elif c.name.startswith(GAIA_NAME):
                out.append(gaia_counts(c, offset))
            else:
                raise ValueError(f"no per-centre builder for {c.name}")
        return out


def marginal_loglike(lnl, weights):
    """ln sum_k w_k L_k for per-centre log-likelihoods ``lnl`` (any leading shape, centres last)."""
    with np.errstate(divide="ignore"):
        lw = np.log(np.asarray(weights, float))
    return logsumexp(np.asarray(lnl, float)+lw, axis=-1)
