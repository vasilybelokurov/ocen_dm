"""Fixed-potential DF approximation tests, separate from mass inference.

Cache actions on a velocity grid in the known potential, then optimize the
same DoublePowerLaw shapes as PositiveDFModel with analytic derivatives.
Component amplitudes here are light-density fractions at an anchor radius,
not gravitating mass fractions. No candidate Poisson solve is performed.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.special import expit, softmax

from .likelihood import ARCSEC_PER_RAD, KinematicData, ProfileLikelihood
from .jeans import KMS_PER_MASYR_KPC
from .positive_df import agama_pc


PARAMETER_NAMES = ("ln_J0", "slope_in", "slope_out", "ln_steepness", "h_r", "g_r")
SHAPE_LOWER = np.array([np.log(5.), -3., 3.6, np.log(.3), .15, .15])
SHAPE_UPPER = np.array([np.log(3000.), 2.5, 25., np.log(8.), 2.85, 2.85])


def capacity_bounds(count, wide=False):
    if count not in (1, 2, 3):
        raise ValueError("capacity ladder supports one, two or three DFs")
    lower = np.array([np.log(.5), -6., 3.3, np.log(.15), .01, .01]) if wide else SHAPE_LOWER
    upper = np.array([np.log(1e5), 2.8, 60., np.log(12.), 2.99, 2.99]) if wide else SHAPE_UPPER
    return (np.r_[np.tile(lower, count), np.full(count-1, -10.)],
            np.r_[np.tile(upper, count), np.full(count-1, 10.)])


def dpl_shape(actions, theta, derivatives=True):
    """Unnormalized native DoublePowerLaw shape and log-shape derivatives.

The arbitrary normalization J0^-3 cancels when anchor-normalizing a component.
AGAMA's (2 pi)^-3 factor is likewise irrelevant here. Verified against native
DF evaluation and finite-difference derivatives in tests.
    """
    actions = np.asarray(actions)
    jr, angular = actions[..., 0], actions[..., 1]+np.abs(actions[..., 2])
    logj, inner, outer, logeta, hr, gr = theta
    eta = np.exp(logeta)
    h = hr*jr+(3-hr)*angular/2
    g = gr*jr+(3-gr)*angular/2
    a, b = eta*(logj-np.log(h)), eta*(np.log(g)-logj)
    la, lb = np.logaddexp(0., a), np.logaddexp(0., b)
    f = np.exp(-3*logj+inner/eta*la-outer/eta*lb)
    if not derivatives:
        return f
    sa, sb = expit(a), expit(b)
    difference = jr-angular/2
    grad = np.stack((-3+inner*sa+outer*sb, la/eta, -lb/eta,
                     inner/eta*(a*sa-la)-outer/eta*(b*sb-lb),
                     -inner*sa*difference/h, -outer*sb*difference/g), axis=-1)
    return f, grad


class FixedPotentialProjector:
    """Cached action/velocity integrals with differentiable positive projection."""

    def __init__(self, potential, radii, radial_nodes=180, velocity_nodes=32,
                 projection_nodes=80, r_min=1e-4, r_max=1e5):
        self.potential = potential
        self.R = np.asarray(radii, float)
        if (np.any(~np.isfinite(self.R)) or np.any(self.R < 10*r_min)
                or np.any(self.R > r_max/100)):
            raise ValueError("projection radii need a wider integration grid")
        if min(radial_nodes, velocity_nodes, projection_nodes) < 8:
            raise ValueError("insufficient quadrature resolution")
        self.settings = dict(radial_nodes=radial_nodes, velocity_nodes=velocity_nodes,
                             projection_nodes=projection_nodes, r_min=r_min, r_max=r_max)
        self.r = np.geomspace(r_min, r_max, radial_nodes)
        self.anchor = -1
        self.anchor_radius = 1.
        integration_r = np.r_[self.r, self.anchor_radius]
        agama = agama_pc()
        af = agama.ActionFinder(potential)
        x, w = np.polynomial.legendre.leggauss(velocity_nodes)
        v, mu, w = (x+1)/2, (x+1)/2, w/2
        escape = np.sqrt(-2*potential.potential(np.column_stack((integration_r, integration_r*0, integration_r*0))))
        pv = np.zeros((radial_nodes+1, velocity_nodes, velocity_nodes, 6))
        pv[..., 0] = integration_r[:, None, None]
        speed = escape[:, None, None]*v[None, :, None]
        pv[..., 3] = speed*mu[None, None, :]
        pv[..., 4] = speed*np.sqrt(1-mu[None, None, :]**2)
        self.actions = af(pv.reshape(-1, 6)).reshape(radial_nodes+1, -1, 3)
        if np.any(~np.isfinite(self.actions)):
            raise ValueError("nonfinite cached actions")
        phase_weight = 4*np.pi*escape[:, None, None]**3*v[None, :, None]**2*w[None, :, None]*w[None, None, :]
        self.weights = np.stack((phase_weight, phase_weight*pv[..., 3]**2,
                                 phase_weight*pv[..., 4]**2/2), axis=-1).reshape(radial_nodes+1, -1, 3)
        x, w = np.polynomial.legendre.leggauss(projection_nodes)
        umax = np.arccosh(r_max/self.R)
        rr = self.R[:, None]*np.cosh(umax[:, None]*(x+1)/2)
        self.los_weights = rr*umax[:, None]*w
        self.q = (self.R[:, None]/rr)**2
        # Local cubic interpolation of log(moment) on four log-radius nodes.
        # It preserves positive moments and has explicit derivatives.
        t = (np.log(rr)-np.log(r_min))/np.log(r_max/r_min)*(radial_nodes-1)
        base = np.clip(np.floor(t).astype(int)-1, 0, radial_nodes-4)
        self.indices = base[..., None]+np.arange(4)
        d = t-base
        self.interpolation_weights = np.stack(
            [np.prod([(d-j)/(i-j) for j in range(4) if j != i], axis=0)
             for i in range(4)], axis=-1)

    def intrinsic(self, x, count, derivatives=True):
        if len(x) != 7*count-1 or np.any(~np.isfinite(x)):
            raise ValueError("invalid mixture parameter vector")
        theta = np.asarray(x[:6*count]).reshape(count, 6)
        fractions = softmax(np.r_[0., x[6*count:]])
        components, jacobians, anchors = [], [], []
        for p in theta:
            if derivatives:
                f, dlog = dpl_shape(self.actions, p)
            else:
                f = dpl_shape(self.actions, p, False)
            moment = np.einsum("rq,rqc->rc", f, self.weights)
            anchor = moment[self.anchor, 0]
            anchors.append(anchor)
            if derivatives:
                dm = np.einsum("rq,rqk,rqc->rck", f, dlog, self.weights)
                dm = (dm-moment[..., None]*dm[self.anchor, 0]/anchor)/anchor
                jacobians.append(dm)
            components.append(moment/anchor)
        total = np.einsum("i,irc->rc", fractions, components)
        if np.any(~np.isfinite(total)) or np.any(total <= 0):
            raise ValueError("invalid fixed-potential moments")
        if not derivatives:
            return total[:-1]
        jac = np.zeros(total.shape+(len(x),))
        for i in range(count):
            jac[..., 6*i:6*(i+1)] = fractions[i]*jacobians[i]
        for i in range(1, count):
            jac[..., 6*count+i-1] = fractions[i]*(components[i]-total)
        return total[:-1], jac[:-1]

    def project(self, x, count, derivatives=True):
        if derivatives:
            moment, jac = self.intrinsic(x, count)
            dlog = jac/moment[..., None]
        else:
            moment = self.intrinsic(x, count, False)
        logm = np.log(moment)
        interp = np.zeros(self.q.shape+(3,))
        if derivatives:
            di = np.zeros(interp.shape+(len(x),))
        for j in range(4):
            weight = self.interpolation_weights[..., j, None]
            interp += weight*logm[self.indices[..., j]]
            if derivatives:
                di += weight[..., None]*dlog[self.indices[..., j]]
        interp = np.exp(interp)
        # Output order: surface density, LOS pressure, PM radial, PM tangential.
        rho, pr, pt = np.moveaxis(interp, -1, 0)
        q, w = self.q, self.los_weights
        out = np.stack((np.sum(w*rho, axis=1), np.sum(w*((1-q)*pr+q*pt), axis=1),
                        np.sum(w*(q*pr+(1-q)*pt), axis=1), np.sum(w*pt, axis=1)), axis=-1)
        if not derivatives:
            return out
        di *= interp[..., None]
        drho, dpr, dpt = np.moveaxis(di, -2, 0)
        q, w = q[..., None], w[..., None]
        dj = np.stack((np.sum(w*drho, axis=1), np.sum(w*((1-q)*dpr+q*dpt), axis=1),
                       np.sum(w*(q*dpr+(1-q)*dpt), axis=1), np.sum(w*dpt, axis=1)), axis=1)
        return out, dj

    def native_df(self, x, count):
        """Rebuild exactly the fitted shape in AGAMA for independent checks."""
        agama = agama_pc()
        weights = softmax(np.r_[0., x[6*count:]])
        dfs = []
        for i, theta in enumerate(np.asarray(x[:6*count]).reshape(count, 6)):
            f = dpl_shape(self.actions, theta, False)
            anchor = np.sum(f[self.anchor]*self.weights[self.anchor, :, 0])
            logj, inner, outer, logeta, hr, gr = theta
            dfs.append(agama.DistributionFunction(type="DoublePowerLaw", norm=(2*np.pi)**3*weights[i]/anchor,
                       J0=np.exp(logj), slopeIn=inner, slopeOut=outer, steepness=np.exp(logeta),
                       coefJrIn=hr, coefJzIn=(3-hr)/2, coefJrOut=gr, coefJzOut=(3-gr)/2))
        return dfs[0] if count == 1 else agama.DistributionFunction(*dfs)


class CapacityObservations:
    """Noise-free mock values at saved observational bins; real error scales.

Streaming is zero by construction. Photometry has the pilot's explicitly
adopted error scale, so the combined score is diagnostic, not a likelihood
calibrated for model selection. All observations sample the same light DF.
    """

    def __init__(self, template, photometry, truth, distance_kpc=5.43):
        self.data = KinematicData(tuple(replace(p, streaming2=None) for p in template.profiles))
        self.photo = photometry
        self.distance_kpc = distance_kpc
        self.likelihood = ProfileLikelihood(self.data)
        self.nodes = [self.likelihood._nodes(p, distance_kpc*1000/ARCSEC_PER_RAD)
                      for p in self.data.profiles]
        self.radii = np.concatenate([r for r, _ in self.nodes]+[
            photometry.r_arcsec*distance_kpc*1000/ARCSEC_PER_RAD])
        truth._check_inside_tracer(self.radii)
        expected = truth.projected_moments(self.radii)
        self.truth_projected = np.column_stack([expected[k] for k in ("Sigma", "los", "pmr", "pmt")])
        self.truth = self.reduce(self.truth_projected)
        self.errors = np.concatenate([(p.err_lo+p.err_hi)/2 for p in self.data.profiles]+[photometry.sigma_mag])
        self.nkin = self.data.n_points
        self.photo_weights = 1/photometry.sigma_mag**2

    def reduce(self, projected, jac=None):
        values, derivs, start = [], [], 0
        for p, (r, w) in zip(self.data.profiles, self.nodes):
            sl = slice(start, start+len(r))
            start += len(r)
            m = projected[sl]
            kind = {"los": 1, "pmr": 2, "pmt": 3}[p.kind]
            s2 = m[:, kind]/m[:, 0]
            if jac is not None:
                dm = jac[sl]
                ds2 = (dm[:, kind]-s2[:, None]*dm[:, 0])/m[:, 0, None]
            if p.r_nodes is not None:
                s2 = s2.reshape(p.n, -1).mean(axis=1)
                if jac is not None:
                    ds2 = ds2.reshape(p.n, -1, jac.shape[-1]).mean(axis=1)
            elif w is not None:
                weight = m[:, 0]*r*w
                denominator = weight.reshape(p.n, -1).sum(axis=1)
                s2 = (weight*s2).reshape(p.n, -1).sum(axis=1)/denominator
                if jac is not None:
                    numerator_deriv = (dm[:, kind]*(r*w)[:, None]).reshape(p.n, -1, jac.shape[-1]).sum(axis=1)
                    denominator_deriv = (dm[:, 0]*(r*w)[:, None]).reshape(p.n, -1, jac.shape[-1]).sum(axis=1)
                    ds2 = (numerator_deriv-s2[:, None]*denominator_deriv)/denominator[:, None]
            conversion = 1. if p.kind == "los" else KMS_PER_MASYR_KPC*self.distance_kpc
            sigma = np.sqrt(s2)/conversion
            values.append(sigma)
            if jac is not None:
                derivs.append(ds2/(2*np.sqrt(s2)[:, None]*conversion))
        values.append(-2.5*np.log10(projected[start:, 0]))
        if jac is not None:
            derivs.append(-2.5/np.log(10)*jac[start:, 0]/projected[start:, 0, None])
            return np.concatenate(values), np.concatenate(derivs)
        return np.concatenate(values)

    def residual(self, projected, jac=None):
        if jac is None:
            prediction = self.reduce(projected)
        else:
            prediction, deriv = self.reduce(projected, jac)
        residual = prediction-self.truth
        residual[self.nkin:] -= np.average(residual[self.nkin:], weights=self.photo_weights)
        if jac is not None:
            deriv[self.nkin:] -= np.average(deriv[self.nkin:], axis=0, weights=self.photo_weights)
            return residual/self.errors, deriv/self.errors[:, None]
        return residual/self.errors

    def to_dict(self):
        return dict(distance_kpc=self.distance_kpc, projected_radii_pc=self.radii.tolist(),
                    truth_projected=self.truth_projected.tolist(), truth_values=self.truth.tolist(),
                    errors=self.errors.tolist(), n_kinematic=self.nkin,
                    streaming="none: nonrotating truth and fit", noise="none",
                    kinematic_errors="mean of saved lower/upper errors; diagnostic weighting",
                    photometric_errors=self.photo.error_model,
                    population_selection="common luminosity-weighted mixture for all profiles")
