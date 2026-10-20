"""Recovery tests for the joint centre fit on small synthetic catalogues (2026-10-20)."""
import numpy as np

from ocen_dm.centre_fit import CountModel, DispersionModel, Stars, fit, symmetry_centre


class FlatDepth:
    def __call__(self, x, y):
        return np.zeros(np.shape(x))


def _mock_positions(n, c, rmax, a=120., b=1.0, rng=None):
    rng = rng or np.random.default_rng(0)
    out = []
    while sum(len(o) for o in out) < n:
        x, y = rng.uniform(-rmax, rmax, (2, 4*n))
        keep = (np.hypot(x, y) < rmax) & (rng.random(4*n) < (1+((x-c[0])**2+(y-c[1])**2)/a**2)**(-b))
        out.append(np.column_stack((x[keep], y[keep])))
    p = np.vstack(out)[:n]
    return p[:, 0], p[:, 1]


def test_counts_recover_injected_centre():
    c = (1.5, -1.0)
    x, y = _mock_positions(200000, c, 150., rng=np.random.default_rng(1))
    s = Stars(x, y, np.full(len(x), 20.))
    r = fit(counts=CountModel(s, FlatDepth(), 150., grid=2.0))
    assert abs(r.centre[0]-c[0]) < 4*r.centre_err[0] and abs(r.centre[1]-c[1]) < 4*r.centre_err[1]
    assert r.centre_err[0] < 1.0


def test_dispersion_recovers_injected_centre():
    rng = np.random.default_rng(2); c = (-2.0, 1.0)
    x, y = rng.uniform(-150, 150, (2, 200000)); k = np.hypot(x, y) < 150; x, y = x[k], y[k]
    R2 = (x-c[0])**2+(y-c[1])**2; sig = 0.8*(1+R2/120.**2)**(-0.3)
    w = R2/(R2+200.**2); sR, sT = sig*np.sqrt(1+0.3*w), sig*np.sqrt(1-0.3*w)
    phi = np.arctan2(y-c[1], x-c[0]); vR, vT = rng.normal(0, sR), rng.normal(0, sT)
    mag = rng.uniform(18., 23., len(x)); e = 0.03*np.exp(0.4*(mag-20.))
    vx = vR*np.cos(phi)-vT*np.sin(phi)+rng.normal(0, e); vy = vR*np.sin(phi)+vT*np.cos(phi)+rng.normal(0, e)
    s = Stars(x, y, mag, vx, vy, e, e, np.ones(len(x), bool))
    r = fit(disp=DispersionModel(s, FlatDepth(), 150.))
    assert abs(r.centre[0]-c[0]) < 4*r.centre_err[0] and abs(r.centre[1]-c[1]) < 4*r.centre_err[1]


def test_symmetry_centre_not_truncated():
    c = (4.5, -3.0)
    x, y = _mock_positions(60000, c, 120., a=60., b=1.5, rng=np.random.default_rng(3))
    cc, edge = symmetry_centre(x, y, 5., 55.)
    assert not edge and np.hypot(cc[0]-c[0], cc[1]-c[1]) < 2.0
