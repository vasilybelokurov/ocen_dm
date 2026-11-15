#!/usr/bin/env python3
"""JAM mass models with an optional dark-matter halo, and fast MGEs of the non-Gaussian components (2026-11-14).

jam_models(jeans, D, qobs, qrem, qdm=1.0) is jam_mock_rotation_test.jam_models plus a halo: every non-Gaussian,
non-point component of the mass model becomes a projected NNLS MGE as before, but
  * the Plummer remnants use the analytic surface density Sigma = M a^2 / (pi (a^2 + R^2)^2) instead of 100 adaptive
    quadratures (same sigma grid, fit radii and normalisation as backends.mge_from_component, so the MGE is the same
    to the quadrature accuracy);
  * a TruncatedGNFW halo is projected with a fixed 256-node Gauss-Legendre rule in u (r = R cosh u), vectorised over
    radius, on wider grids (fit radii 0.005-5000 pc, sigmas 0.005-4000 pc, 56 Gaussians: r_s reaches 1000 pc and the NFW cusp needs small sigmas);
    it is ROUND by default (qdm = 1), a choice, not a measurement.
Remnants take qrem, stars qobs (sigma_major = sigma/sqrt(q'), mass preserved), as in the tuned models.
qobs / qrem may be a number, an array (one per Gaussian) or a QProfile (2026-11-15): the measured projected-shape
profile q'(R) from the star counts (bin/fit_shape_counts.py), evaluated per Gaussian.

Self-test: PYTHONPATH=src:bin python bin/jam_halo_models.py
  1. remnant MGE (analytic) vs backends.mge_from_component (quadrature): max relative difference of enclosed mass
  2. halo MGE enclosed mass vs the exact TruncatedGNFW, 0.1-1000 pc, for r_s = 3, 30, 300, 1000 pc: < 0.5%
  3. timing per call
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
from scipy import optimize

import jam_mock_rotation_test as J
from ocen_dm.kinematics.backends import _split_mass_model
from ocen_dm.mass_models.dark_matter import TruncatedGNFW
from ocen_dm.mass_models.stellar import Plummer

REM_SIG = np.geomspace(0.05, 500., 40)          # as jam_mock_rotation_test.jam_models
REM_FIT = np.geomspace(0.05, 2000., 100)        # as backends.mge_from_component
DM_SIG = np.geomspace(0.005, 4000., 56)
DM_FIT = np.geomspace(0.005, 5000., 140)
_XG, _WG = np.polynomial.legendre.leggauss(256)


class QProfile:
    """Measured q'(R) for JAM (docs/JAM_FIT_PLAN.md step 4f/5). Each Gaussian of circular width s [arcsec] gets the
    isophotal q' of the star-count fit at major-axis radius R = s/sqrt(q') (= its sigma_major; van de Ven et al. 2006
    prescription: q'_j = observed flattening at R = sigma'_j), solved by fixed-point iteration. Zones:
      R < r_in            : qin (free in the JAM fit; no shape constraint in the centre)
      r_in <= R <= r_out  : measured (linear interpolation in log R)
      R > r_out           : outer = 'hold' (q'(r_out)), 'round' (1.0) or 'measured' (measured to 2000", then held)
    Picklable (plain arrays), so it travels to worker processes."""

    def __init__(self, R, q, qin, r_in=60., r_out=1300., outer="hold", src=""):
        self.R, self.q = np.asarray(R, float), np.asarray(q, float)
        self.qin, self.r_in, self.r_out, self.outer, self.src = float(qin), float(r_in), float(r_out), outer, src

    @classmethod
    def from_json(cls, path, qin=0.9, **kw):
        import json
        d = json.loads(Path(path).read_text())
        return cls(d["R_iso"], d["q_iso"], qin, src=Path(path).name, **kw)

    def with_qin(self, qin):
        return QProfile(self.R, self.q, qin, self.r_in, self.r_out, self.outer, self.src)

    def at_R(self, R):
        R = np.asarray(R, float)
        qm = np.interp(np.log(R), np.log(self.R), self.q)
        if self.outer == "hold":
            qo = np.interp(np.log(self.r_out), np.log(self.R), self.q)
        elif self.outer == "round":
            qo = 1.0
        else:
            qo = np.interp(np.log(np.minimum(R, 2000.)), np.log(self.R), self.q)
        return np.where(R < self.r_in, self.qin, np.where(R > self.r_out, qo, qm))

    def __call__(self, s_arcsec):
        s = np.asarray(s_arcsec, float); q = np.full(s.shape, 0.9)
        for _ in range(30):
            q = self.at_R(s/np.sqrt(q))
        return q

    def describe(self):
        return dict(src=self.src, qin=self.qin, r_in=self.r_in, r_out=self.r_out, outer=self.outer)


def qvals(q, s_arcsec):
    """q' per Gaussian: q may be a number, an array or a QProfile (callable of the circular width in arcsec)."""
    return q(s_arcsec) if callable(q) else np.broadcast_to(np.asarray(q, float), np.shape(s_arcsec))


def projected(component, R, r_max):
    """Sigma(R) = 2 int_0^umax rho(R cosh u) R cosh u du (Msun/pc^2), analytic for a Plummer."""
    R = np.asarray(R, float)
    if isinstance(component, Plummer):
        a = component.scale
        return component.mass*a**2/(np.pi*(a**2+R**2)**2)
    umax = np.arccosh(np.maximum(r_max/R, 1.0+1e-9))
    u = 0.5*umax[:, None]*(_XG[None, :]+1)
    r = R[:, None]*np.cosh(u)
    rho = np.asarray(component.density(r.ravel()), float).reshape(r.shape)
    return 2*0.5*umax*np.sum(_WG[None, :]*rho*r, axis=1)


def mge_masses(component, sig, rfit):
    """NNLS projected MGE (masses per sigma), renormalised to the mass inside rfit.max() (backends convention)."""
    S = projected(component, rfit, 10*rfit.max())
    basis = np.exp(-0.5*(rfit[:, None]/sig[None, :])**2)/(2*np.pi*sig[None, :]**2)
    w = 1/np.maximum(S, 1e-300)
    m, _ = optimize.nnls(basis*w[:, None], S*w, maxiter=200*len(sig))
    rm = rfit.max()
    target = float(np.ravel(component.enclosed_mass(rm))[0])
    return m*target/_mge_enclosed(m, sig, np.array([rm]))[0]


def jam_models(jeans, D, qobs, qrem, qdm=1.0):
    """(lum, pot, mbh) for JamPy with stars (qobs), remnants (qrem) and a halo (qdm)."""
    pc = D*1e3/J.ARCSEC_PER_RAD
    gaussians, mbh, other = _split_mass_model(jeans.mass)
    stars = gaussians[0]
    ms, ss = np.asarray(stars.masses)*stars.mass_to_light, np.asarray(stars.sigmas)
    parts = [J.flatten(ms, ss, qvals(qobs, ss/pc), pc)]
    want = float(np.sum(ms))
    for c in other:
        if isinstance(c, TruncatedGNFW):
            m = mge_masses(c, DM_SIG, DM_FIT); q, sg = qdm, DM_SIG
        else:
            m = mge_masses(c, REM_SIG, REM_FIT); q, sg = qrem, REM_SIG
        k = m > 0
        parts.append(J.flatten(m[k], sg[k], qvals(q, sg[k]/pc), pc))
        want += float(np.sum(m[k]))
    pot = tuple(np.concatenate([pp[i] for pp in parts]) for i in range(3))
    got = J.mge_total(*pot, pc)
    assert abs(got/want-1) < 1e-10, (got, want)
    sl, sg, qq = parts[0]
    lum = (np.r_[sl, 1e-10*sl.max()], np.r_[sg, 5*sg.max()], np.r_[qq, qq[-1]])   # see J.jam_models
    return lum, pot, mbh


def _oblate_enclosed(m, sig, q, r):
    """Mass within spheres of radius r of oblate 3D Gaussians rho ~ exp(-(R^2 + z^2/q^2)/(2 sig^2)), total masses m,
    major-axis sigmas sig, intrinsic axis ratios q (< 1, or 1):
    M(<r) = m [erf(r/(sqrt2 sig q)) - exp(-r^2/(2 sig^2)) erf(r e/(sqrt2 sig q))/e], e = sqrt(1-q^2);
    the e -> 0 limit is the spherical Gaussian."""
    from scipy.special import erf
    r = np.asarray(r, float)[:, None]; m, sig, q = (np.broadcast_to(np.asarray(v, float), np.shape(m)) for v in (m, sig, q))
    e = np.sqrt(np.clip(1-q**2, 0, None))
    x = r/(np.sqrt(2)*sig*q)
    sph = np.sqrt(2/np.pi)*(r/sig)*np.exp(-0.5*(r/sig)**2)
    tail = np.where(e > 1e-6, np.exp(-0.5*(r/sig)**2)*erf(x*np.where(e > 1e-6, e, 1.))/np.where(e > 1e-6, e, 1.), sph)
    return np.sum(m*(erf(x)-tail), axis=1)


def enclosed_flat(jeans, qobs, qrem, inc, r, qdm=1.0, D=None):
    """Spherical enclosed mass M(<r) [Msun] per component of the flattened JAM model built by jam_models: each
    Gaussian deprojected at inclination inc (deg), intrinsic q = sqrt(q'^2 - cos^2 i)/sin i, major sigma
    s/sqrt(q') (J.flatten). r in pc. Returns {component class name: array}; the BH is 'BH'. D [kpc] is needed only
    when qobs/qrem is a QProfile (widths converted to arcsec)."""
    pc = None if D is None else D*1e3/J.ARCSEC_PER_RAD
    if pc is None and (callable(qobs) or callable(qrem)):
        raise ValueError("enclosed_flat: D is required with a QProfile")
    ci, si = np.cos(np.radians(inc)), np.sin(np.radians(inc))
    qint = lambda qp: np.sqrt(np.clip(qp**2-ci**2, 1e-12, None))/si
    gaussians, mbh, other = _split_mass_model(jeans.mass)
    stars = gaussians[0]
    ms, ss = np.asarray(stars.masses)*stars.mass_to_light, np.asarray(stars.sigmas)
    qs = qvals(qobs, ss/pc if pc else ss)
    out = {"stars": _oblate_enclosed(ms, ss/np.sqrt(qs), qint(qs), r)}
    for c in other:
        if isinstance(c, TruncatedGNFW):
            m, q, sg = mge_masses(c, DM_SIG, DM_FIT), qdm, DM_SIG
        else:
            m, q, sg = mge_masses(c, REM_SIG, REM_FIT), qrem, REM_SIG
        k = m > 0
        qk = qvals(q, sg[k]/pc if pc else sg[k])
        out[type(c).__name__] = out.get(type(c).__name__, 0)+_oblate_enclosed(m[k], sg[k]/np.sqrt(qk), qint(qk), r)
    out["BH"] = np.full(len(np.atleast_1d(r)), float(mbh))
    return out


def _mge_enclosed(m, sig, r):
    from scipy.special import erf
    x = r[:, None]/sig[None, :]
    return np.sum(m[None, :]*(erf(x/np.sqrt(2))-np.sqrt(2/np.pi)*x*np.exp(-0.5*x**2)), axis=1)


def selftest():
    from ocen_dm.kinematics.backends import mge_from_component
    from ocen_dm.mass_models.remnants import RemnantPlummer
    r = np.geomspace(0.1, 1000., 60)
    print("1. remnant MGE, analytic vs quadrature (enclosed mass, max |rel|):")
    for a in (0.5, 2., 5., 15.):
        c = RemnantPlummer(1e6, a)
        ma = mge_masses(c, REM_SIG, REM_FIT)
        g = mge_from_component(c, REM_SIG)
        mq = np.zeros_like(REM_SIG); mq[np.isin(REM_SIG, g.sigmas)] = g.masses
        ea, eq, ex = _mge_enclosed(ma, REM_SIG, r), _mge_enclosed(mq, REM_SIG, r), np.ravel(c.enclosed_mass(r))
        print(f"   a_rem {a:5.1f} pc: analytic vs quad {np.max(np.abs(ea/eq-1)):.1e}; vs exact {np.max(np.abs(ea/ex-1)):.1e}")
    print("2. halo MGE vs exact TruncatedGNFW (r_t 1000 pc), enclosed mass 0.1-1000 pc, max |rel| (pass < 5e-3):")
    ok = True
    for rs in (3., 30., 300., 1000.):
        h = TruncatedGNFW(rho_s=1.0, r_s=rs, gamma=1.0, r_t=1000.)
        m = mge_masses(h, DM_SIG, DM_FIT)
        e = np.max(np.abs(_mge_enclosed(m, DM_SIG, r)/np.ravel(h.enclosed_mass(r))-1))
        ok &= e < 5e-3
        print(f"   r_s {rs:6.0f} pc: {e:.1e}")
    t0 = time.perf_counter()
    for _ in range(20):
        mge_masses(TruncatedGNFW(rho_s=1.0, r_s=300., gamma=1.0, r_t=1000.), DM_SIG, DM_FIT)
        mge_masses(RemnantPlummer(1e6, 5.), REM_SIG, REM_FIT)
    print(f"3. halo + remnant MGE: {(time.perf_counter()-t0)/20*1e3:.1f} ms per call")
    print("PASS" if ok else "FAIL")


if __name__ == "__main__":
    selftest()
