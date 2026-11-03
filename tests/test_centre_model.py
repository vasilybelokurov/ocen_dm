"""Model-side centre marginalisation (src/ocen_dm/centre_model.py): geometry and factors."""
import numpy as np

from ocen_dm.centre_model import ModelMoments, PointSet, effective_sample_size, projected_variance


def _ps(seed=0, n=20000, lo=(2., 10.), hi=(10., 40.)):
    rng = np.random.default_rng(seed)
    r = np.sqrt(rng.uniform(lo[0]**2, hi[-1]**2, n))
    p = rng.uniform(0, 2*np.pi, n)
    return PointSet.from_positions(r*np.sin(p), r*np.cos(p), lo, hi)


def test_isotropic_constant_moments_give_unit_factors():
    mom = ModelMoments.from_functions(s2r=lambda r: 0*r+2., s2t=lambda r: 0*r+2., s2los=lambda r: 0*r+3.,
                                      Sigma=lambda r: 0*r+1.)
    ps = _ps()
    for comp in ("r", "t", "los"):
        for c in ((0., 0.), (3., -2.)):
            v = projected_variance(mom, ps.x, ps.y, c, comp)
            assert np.allclose(ps.mean(v), 3. if comp == "los" else 2., rtol=1e-12)


def test_rotation_of_axes():
    mom = ModelMoments.from_functions(s2r=lambda r: 0*r+1., s2t=lambda r: 0*r+3.)
    x, y = np.array([10.]), np.array([0.])
    # about c = (0, 10) the star's radial direction is (1, -1)/sqrt2: 45 deg from the data axis
    assert np.isclose(projected_variance(mom, x, y, (0., 10.), "r")[0], 2.)
    assert np.isclose(projected_variance(mom, x, y, (0., 10.), "t")[0], 2.)
    assert np.isclose(projected_variance(mom, x, y, (0., 0.), "r")[0], 1.)
    assert np.isclose(projected_variance(mom, x, y, (0., 0.), "t")[0], 3.)


def test_annuli_quadrature_exact_for_r2():
    lo, hi = np.array([0., 5.]), np.array([5., 20.])
    ps = PointSet.annuli(lo, hi)
    c = np.array([1.3, -0.7])
    m = ps.mean((ps.x-c[0])**2+(ps.y-c[1])**2)
    expect = 0.5*(hi**2+np.maximum(lo, 1e-3*hi)**2)+c @ c
    assert np.allclose(m, expect, rtol=1e-10)


def test_bin_assignment_and_drop():
    ps = PointSet.from_positions(np.array([1., 3., 12., 50.]), np.zeros(4), [2., 10.], [10., 40.])
    assert ps.idx.tolist() == [0, 1] and ps.x.tolist() == [3., 12.]


def test_ess():
    ess, amax = effective_sample_size([0., 0., 0., 0.], [.25]*4)
    assert np.isclose(ess, 4.) and np.isclose(amax, .25)
    ess, amax = effective_sample_size([0., -50., -50.], [1/3]*3)
    assert np.isclose(ess, 1., atol=1e-12)


def _cs(lo=(0., 4.), hi=(4., 9.)):
    from ocen_dm.centre_model import ContinuousSet
    rng = np.random.default_rng(1)
    radii = [np.sqrt(rng.uniform(max(a, 0)**2, b**2, 500)) for a, b in zip(lo, hi)]
    cs = ContinuousSet(lo, hi, radii)
    cs.gamma[:] = 0.
    return cs


def test_continuous_exact_for_r2():
    cs = _cs()
    c = (1.7, -2.2)
    m = cs.mean(lambda x, y, cc: (x-cc[0])**2+(y-cc[1])**2, c)
    lo = np.array([0., 4.]); hi = np.array([4., 9.])
    expect = 0.5*(hi**2+np.maximum(lo, 1e-3*hi)**2)+c[0]**2+c[1]**2
    assert np.allclose(m, expect, rtol=2e-3)


def test_continuous_cusp_is_smooth_and_converged():
    """<1/|x - c|> over a disc/ring: no spikes along a line of centres; agrees with a 4x finer quadrature."""
    from ocen_dm.centre_model import ContinuousSet
    cs = _cs()
    fine = ContinuousSet(cs.lo, cs.hi, [np.array([1., 2.]), np.array([5., 6.])], n_panel=192, n_theta=768)
    fine.gamma[:] = 0.
    f = lambda x, y, cc: 1/np.maximum(np.hypot(x-cc[0], y-cc[1]), 1e-6)
    xs = np.linspace(0., 6., 61)
    a = np.array([cs.mean(f, (x, 0.3)) for x in xs])
    for x, row in zip(xs[::10], a[::10]):
        assert np.allclose(row, fine.mean(f, (x, 0.3)), rtol=2e-3)
    d2 = np.abs(np.diff(a, 2, axis=0))
    assert d2.max() < 0.02*np.abs(a).max()


def test_gamma_fit():
    from ocen_dm.centre_model import ContinuousSet
    rng = np.random.default_rng(2)
    lo, hi, g = 3., 12., 0.8                       # p(R) ~ R^(1-g)
    a = 2-g
    r = (lo**a+rng.uniform(size=200000)*(hi**a-lo**a))**(1/a)
    assert abs(ContinuousSet._fit_gamma(r, lo, hi)-g) < 0.02
