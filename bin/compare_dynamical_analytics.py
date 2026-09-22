"""Compare the completed pilot with energy-cut and pericentre approximations.

Read-only with respect to every simulation. Figures are PNGs in plots/; the
JSON in results/plot_data/ records inputs, definitions, checks and unrounded values.
No analytical energy input is interpreted as a fractional density loss.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar, brentq
from scipy.special import roots_legendre
import yaml

from ocen_dm.dynamics import G, TIME_UNIT_MYR
from ocen_dm.dynamics.experiment import load_prepared
from ocen_dm.dynamics.models import (
    agama_module, density_parameters, density_log_derivatives, signed_df_grid, xyz,
)
from ocen_dm.dynamics.diagnostics import tidal_tensor
from ocen_dm.tails.truncated_equilibrium import kepler_truncation_factor, shock_heating

ROOT = Path(__file__).resolve().parents[1]
PASSAGE = ROOT / "results/dynamics/remnant_passage_20260921T114939Z"
ANALYSIS = ROOT / "results/dynamics/batch_analysis_20260921_v2"
RADII = np.array([3., 10., 20., 35., 50., 70.])
COLOURS = {"cusp": "#2869a4", "core": "#c25e24"}


def digest(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def density_integral(psi, binding_cut, f_binding, order=256):
    """4 pi integral v^2 f(psi-v^2/2) dv, with binding energy >= cut.

    This is an isotropic, fixed-potential DF deletion, without renormalization
    or re-equilibration. Units follow the supplied potential and DF.
    """
    psi = np.atleast_1d(psi)
    nodes, weights = roots_legendre(order)
    vmax = np.sqrt(2*np.maximum(psi-binding_cut, 0))
    speed = vmax[:, None]*(nodes+1)/2
    value = f_binding(psi[:, None]-.5*speed**2)
    return 4*np.pi*np.sum(speed**2*value*weights, axis=1)*vmax/2


def shell_retention(radii_pc, point_density, point_cut_density, order=24):
    nodes, weights = roots_legendre(order)
    lo, hi = np.asarray(radii_pc)/1000*np.exp(-.1), np.asarray(radii_pc)/1000*np.exp(.1)
    rr = lo[:, None]+(hi-lo)[:, None]*(nodes+1)/2
    numerator = np.sum(rr**2*point_cut_density(rr.ravel()).reshape(rr.shape)*weights, axis=1)
    denominator = np.sum(rr**2*point_density(rr.ravel()).reshape(rr.shape)*weights, axis=1)
    return numerator/denominator


def main():
    agama = agama_module()
    output = ROOT/"plots"
    data_output = ROOT/"results/plot_data"
    data_output.mkdir(parents=True, exist_ok=True)
    audit = json.loads((ANALYSIS/"summary.json").read_text())
    hashes = {}

    def record(path):
        key = str(path.relative_to(ROOT))
        value = digest(path)
        if key in audit["inputs_sha256"]:
            assert value == audit["inputs_sha256"][key], f"Audited input changed: {key}"
        hashes[key] = value

    def read_diagnostics(prepared, engine):
        directory = prepared/f"evolution_{engine}"
        manifest = yaml.safe_load((directory/"run.yaml").read_text())
        assert manifest["status"] == "complete"
        path = directory/"diagnostics.ecsv"
        assert digest(path) == manifest["diagnostics_sha256"]
        record(path)
        record(directory/"run.yaml")
        return Table.read(path)

    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), models={},
                  definitions=dict(
                      density="Final-10-Myr median tide/isolation shell-density ratio; bands are temporal p16-p84",
                      energy="Frozen tidal minus matched frozen isolation specific energy in initial satellite potential; initial shell membership retained, including escapers",
                      heating_proxy="2 mean(delta E)/(G M_n/r); not delta E divided by true initial binding energy",
                      energy_error="Particle-sampling standard error of the endpoint mean for this fixed potential and realization; not model or convergence uncertainty",
                      truncation="Delete f(E) for E>Phi(r_t) in initial spherical potential, no renormalization or self-consistent re-equilibration",
                      shock="One scalar pericentre estimate, not a disc-shock model; characteristic energy G M_n/(2r)"))
    for shape in COLOURS:
        prepared = PASSAGE/f"remnant_{shape}_passage"
        isolated = PASSAGE/f"remnant_{shape}_isolation_reference"
        config, particles, track, manifest = load_prepared(prepared)
        _, control_particles, control_track, _ = load_prepared(isolated)
        for field in ("xv", "mass", "particle_id", "component"):
            assert np.array_equal(particles[field], control_particles[field])
        for directory in (prepared, isolated):
            for name in ("run.yaml", "config.yaml", "initial_local.npz", "orbit.npz",
                         "satellite_initial.ini", "nucleus.ini", "dm_potential.ini", "equilibrium.json"):
                record(directory/name)
        record(prepared/"host.ini")
        potential = agama.Potential(str(prepared/"satellite_initial.ini"))
        host = agama.Potential(str(prepared/"host.ini"))
        equilibrium = json.loads((prepared/"equilibrium.json").read_text())
        component = config.components[0]
        assert component.beta == 0 and not config.nucleus.live
        density = agama.Density(**density_parameters(component, equilibrium["components"]["dm"]["density_norm_msun_kpc3"]))

        i = np.argmin(np.linalg.norm(track.xv[:, :3], axis=1))
        peri_time = minimize_scalar(lambda tm: np.linalg.norm(track.at(tm)[:3]),
                                   bounds=(track.time[i-1], track.time[i+1]), method="bounded",
                                   options={"xatol": 1e-14}).x
        phase = track.at(peri_time)
        radius = np.linalg.norm(phase[:3])
        speed = np.linalg.norm(phase[3:])
        radial_acc = -np.dot(host.force(phase[:3], t=peri_time), phase[:3])/radius
        equivalent_mass = radial_acc*radius**2/G
        # Initial-potential tidal scale: forecast does not use a depleted mass profile.
        def tidal_radius(tm):
            phase = track.at(tm)
            xx, vv = phase[:3], phase[3:]
            omega = np.cross(xx, vv)/np.dot(xx, xx)
            tensor = tidal_tensor(host, xx, tm)+np.dot(omega, omega)*np.eye(3)-np.outer(omega, omega)
            lam = np.linalg.eigvalsh(tensor)[-1]
            return brentq(lambda rr: -potential.force([rr, 0, 0])[0]/rr-lam, .001, 1.)
        trial_times = np.linspace(peri_time-.005, peri_time+.005, 101)
        trial_r = np.array([tidal_radius(tm) for tm in trial_times])
        j = np.argmin(trial_r)
        rt_time = minimize_scalar(tidal_radius, bounds=(trial_times[j-1], trial_times[j+1]),
                                  method="bounded", options={"xatol": 1e-12}).x
        rt = tidal_radius(rt_time)
        binding_cut = -float(potential.potential([rt, 0, 0]))
        df_checks = []
        # Independent signed inversion; double the energy grid and quadrature.
        for n_energy, order in [(512, 192), (1024, 384)]:
            df = signed_df_grid(potential, density, lambda rr: density_log_derivatives(component, rr),
                                0., 1e-5, 1., n_energy=n_energy, quadrature_order=order)
            assert np.min(df["signed_over_absolute"]) >= 0
            ee, ff = df["binding_energy_kms2"][::-1], df["f_energy"][::-1]
            positive = ff > 0
            log_f = PchipInterpolator(np.log(ee[positive]), np.log(ff[positive]), extrapolate=False)

            def f_binding(binding):
                inside = (binding >= ee[positive][0]) & (binding <= ee[positive][-1])
                out = np.zeros_like(binding)
                out[inside] = np.exp(log_f(np.log(binding[inside])))
                return out

            def rho(rr):
                return density_integral(-potential.potential(xyz(rr)), 0., f_binding, order)

            def rho_cut(rr):
                return density_integral(-potential.potential(xyz(rr)), binding_cut, f_binding, order)

            smooth_r = np.geomspace(2.5, 75, 180)
            retention = shell_retention(smooth_r, rho, rho_cut)
            bins = shell_retention(RADII, rho, rho_cut)
            checks_r = np.geomspace(2.2, 85, 100)/1000
            recovery = np.max(np.abs(rho(checks_r)/density.density(xyz(checks_r))-1))
            df_checks.append(dict(n_energy=n_energy, quadrature_order=order,
                                  max_density_recovery_error=float(recovery), retained=bins.tolist()))
        quadrature_difference = np.max(np.abs(np.array(df_checks[0]["retained"])-bins))
        assert recovery < .003 and quadrature_difference < .001

        observed = {}
        for engine in ("live", "frozen"):
            tide, control = read_diagnostics(prepared, engine), read_diagnostics(isolated, engine)
            assert np.allclose(tide["elapsed_myr"], control["elapsed_myr"], atol=1e-6, rtol=0)
            late = tide["elapsed_myr"] >= tide["elapsed_myr"][-1]-10
            observed[engine] = []
            for rr in RADII:
                key = f"dm_{rr:g}pc_rho_msun_pc3"
                ratio = np.asarray(tide[key])/control[key]
                p16, median, p84 = np.percentile(ratio[late], [16, 50, 84])
                observed[engine].append(dict(radius_pc=rr, median=median, p16=p16, p84=p84,
                                             initial_shell_neff=float(tide[f"dm_{rr:g}pc_shell_neff"][0])))

        initial = particles["xv"]
        radii = np.linalg.norm(initial[:, :3], axis=1)*1000
        mass = particles["mass"]
        assert np.all(mass == mass[0]), "Sampling-error formula assumes equal particle weights"
        e0 = potential.potential(initial[:, :3])+.5*np.sum(initial[:, 3:]**2, axis=1)
        apocentre_pc = potential.Rperiapo(initial)[:, 1]*1000
        assert np.all(np.isfinite(apocentre_pc))
        snapshots = [p/"evolution_frozen/snapshots.h5" for p in (prepared, isolated)]
        for path in snapshots:
            record(path)
        with h5py.File(snapshots[0]) as tidal, h5py.File(snapshots[1]) as iso:
            for f in (tidal, iso):
                assert f.attrs["snapshots_written"] == len(f["time"])
                assert np.array_equal(f["particle_id"][:], particles["particle_id"])
                assert np.array_equal(f["mass"][:], mass)
            times = tidal["time"][:]
            assert np.allclose(times-track.time[0], iso["time"][:]-control_track.time[0], atol=1e-12)
            late = np.flatnonzero(times >= times[-1]-10/TIME_UNIT_MYR)
            changes, iso_changes = [], []
            for n in late:
                xx = tidal["xv"][n]-track.at(times[n])
                zz = iso["xv"][n]-control_track.at(iso["time"][n])
                et = potential.potential(xx[:, :3])+.5*np.sum(xx[:, 3:]**2, axis=1)
                ei = potential.potential(zz[:, :3])+.5*np.sum(zz[:, 3:]**2, axis=1)
                changes.append(et-ei)
                iso_changes.append(ei-e0)
        changes = np.asarray(changes)
        iso_changes = np.asarray(iso_changes)
        energy_rows = []
        for rr in RADII:
            shell = (radii >= rr*np.exp(-.1)) & (radii < rr*np.exp(.1))
            weights = mass[shell]
            de = changes[-1, shell]
            # Reproduce the previously audited median before subtracting the
            # small isolation error. This catches frame/energy mismatches.
            assert np.isclose(np.median((et-e0)[shell]),
                              tide[f"dm_initial_{rr:g}pc_median_delta_reference_energy_kms2"][-1],
                              rtol=1e-7, atol=1e-7)
            shot = shock_heating(radii[shell], equivalent_mass, radius*1000, speed,
                                 m_nuc=config.nucleus.mass_msun)
            impulse = .5*shot["dv"]**2
            predicted = np.average(impulse*shot["A"], weights=weights)
            uncorrected = np.average(impulse, weights=weights)
            avg = np.average(de, weights=weights)
            means = np.average(changes[:, shell], weights=weights, axis=1)
            neff = weights.sum()**2/np.sum(weights**2)
            sem = np.sqrt(np.average((de-avg)**2, weights=weights)/(neff-1))
            local = apocentre_pc[shell] < 2*rr
            late_particle = np.mean(changes[:, shell], axis=0)
            positive_energy = et[shell] >= 0
            energy_rows.append(dict(radius_pc=rr, n_particles=int(shell.sum()), mean_delta_e=avg,
                                   median_delta_e=float(np.median(de)), sampling_sem=sem,
                                   late_mean_p16=float(np.percentile(means, 16)),
                                   late_mean_median=float(np.median(means)),
                                   late_mean_p84=float(np.percentile(means, 84)),
                                   predicted_delta_e=predicted, impulse_delta_e=uncorrected,
                                   measured_to_predicted=float(avg/predicted),
                                   mean_initial_binding=float(np.average(-e0[shell], weights=weights)),
                                   fractional_binding_change=float(avg/np.average(-e0[shell], weights=weights)),
                                   fraction_reference_positive_energy=float(np.average(positive_energy, weights=weights)),
                                   mean_delta_e_reference_negative_energy=float(np.average(de[~positive_energy], weights=weights[~positive_energy])),
                                   positive_energy_share_of_net_increment=float(np.sum(weights[positive_energy]*de[positive_energy])/np.sum(weights*de)),
                                   max_abs_iso_mean_delta_e=float(np.max(np.abs(np.average(iso_changes[:, shell], weights=weights, axis=1)))),
                                   median_apocentre_pc=float(np.median(apocentre_pc[shell])),
                                   fraction_apocentre_beyond_2r=float(np.average(~local, weights=weights)),
                                   fraction_apocentre_beyond_rt=float(np.average(apocentre_pc[shell] > rt*1000, weights=weights)),
                                   late_mean_local_apocentre=float(np.average(late_particle[local], weights=weights[local])),
                                   late_mean_extended_apocentre=float(np.average(late_particle[~local], weights=weights[~local]))))
        result["models"][shape] = dict(
            pericentre=dict(elapsed_myr=(peri_time-track.time[0])*TIME_UNIT_MYR,
                            radius_kpc=radius, speed_kms=speed, force_equivalent_mass_msun=equivalent_mass),
            tidal_radius_pc=rt*1000, tidal_min_elapsed_myr=(rt_time-track.time[0])*TIME_UNIT_MYR,
            density_quadrature_checks=df_checks, max_retention_grid_difference=quadrature_difference,
            smooth_radius_pc=smooth_r.tolist(), matched_df_retained=retention.tolist(),
            bin_matched_df_retained=bins.tolist(), observed_density=observed, energy=energy_rows,
            energy_late_snapshots=len(late), end_elapsed_myr=(times[-1]-times[0])*TIME_UNIT_MYR)
        print(shape, "r_t", rt*1000, "DF errors", recovery, quadrature_difference, flush=True)
        print("20 pc", energy_rows[2], flush=True)

    plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5, "legend.fontsize": 7,
                         "axes.grid": True, "grid.alpha": .2, "lines.linewidth": 1.3,
                         "savefig.bbox": "tight"})
    figures = []

    def save(fig, name):
        path = output/f"dynamical_analytics_{name}.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        figures.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(6.65, 3.45), sharey=True, layout="constrained")
    for ax, (shape, colour) in zip(axes, COLOURS.items()):
        model = result["models"][shape]
        rr = np.asarray(model["smooth_radius_pc"])
        rt = model["tidal_radius_pc"]
        ax.axhline(1, color=".6", lw=.7)
        ax.plot(rr, model["matched_df_retained"], color="black", label="Energy cut: actual initial DF")
        ax.plot(rr, np.exp(-rr/rt), color=".45", ls="--", label=r"Earlier $e^{-r/r_t}$ placeholder")
        if shape == "cusp":
            ax.plot(rr, kepler_truncation_factor(rr/rt), color="#985ba5", ls=":", label=r"Kepler, pure $r^{-1}$ cusp")
        for engine, marker, offset in [("live", "o", .97), ("frozen", "s", 1.03)]:
            rows = model["observed_density"][engine]
            y = np.array([v["median"] for v in rows])
            ax.errorbar(RADII*offset, y, yerr=[y-[v["p16"] for v in rows], [v["p84"] for v in rows]-y],
                        fmt=marker, ms=4, color=colour, mfc=colour if engine == "live" else "white",
                        capsize=2, label=f"Simulation: {engine}")
        ax.set(xscale="log", xlabel="Radius [pc]", title=f"{shape.title()}: one 88-Myr passage", ylim=(0, 1.13))
        ax.set_xticks([3, 10, 20, 35, 70], ["3", "10", "20", "35", "70"])
        ax.legend(loc="lower left", fontsize=6.5)
    axes[0].set_ylabel("Retained shell density relative to isolation")
    save(fig, "density")

    fig, axes = plt.subplots(2, 2, figsize=(6.65, 5.1), sharex=True, layout="constrained")
    for col, (shape, colour) in enumerate(COLOURS.items()):
        model = result["models"][shape]
        # The 3-pc shell is inside the Plummer core: the point-mass frequency
        # in the legacy shock estimate is inapplicable there.
        rows = model["energy"][1:]
        energy_radii = RADII[1:]
        rr = np.geomspace(9, 77, 200)
        peri = model["pericentre"]
        shock = shock_heating(rr, peri["force_equivalent_mass_msun"], peri["radius_kpc"]*1000,
                              peri["speed_kms"])
        imp = .5*shock["dv"]**2
        ax = axes[0, col]
        ax.fill_between(rr, imp*(1+shock["omega_tau"]**2)**(-2.5), imp*shock["A"], color=".8",
                        label=r"Scalar pericentre, $\gamma=1.5$--2.5")
        ax.plot(rr, imp*shock["A"], color="black", lw=1)
        ax.plot(rr, imp, ls="--", color=".4", label="Scalar pericentre, no correction")
        ax.errorbar(energy_radii, [v["mean_delta_e"] for v in rows], yerr=[v["sampling_sem"] for v in rows],
                    fmt="o", color=colour, ms=4, capsize=2, label="Frozen full orbit: mean")
        ax.plot(energy_radii, [v["median_delta_e"] for v in rows], "s:", color=colour, ms=4,
                mfc="white", label="Frozen full orbit: median")
        ax.plot(energy_radii, [v["mean_delta_e_reference_negative_energy"] for v in rows], "^--",
                color=colour, ms=4, mfc="white", label=r"Mean of subset with $E_{\rm ref}<0$")
        ax.set(xscale="log", title=shape.title())
        ax.set_yscale("log")
        ax.set_ylim(1e-8, 100)
        ax.legend(loc="lower right", fontsize=6.3)
        ratio = np.array([v["measured_to_predicted"] for v in rows])
        error = [v["sampling_sem"]/v["predicted_delta_e"] for v in rows]
        axes[1, col].errorbar(energy_radii, ratio, yerr=error, fmt="o-", color=colour, capsize=2, ms=4)
        axes[1, col].axhspan(.5, 2., color=".92")
        axes[1, col].axhline(1, color=".5", lw=.7)
        axes[1, col].set(xlabel="Initial shell radius [pc]")
        axes[1, col].set_yscale("log")
        axes[1, col].set_ylim(.3, 1e6)
        axes[1, col].set_xticks([10, 20, 35, 50, 70], ["10", "20", "35", "50", "70"])
    axes[0, 0].set_ylabel(r"Specific energy change [km$^2$ s$^{-2}$]")
    axes[1, 0].set_ylabel(r"Mean / scalar prediction ($\gamma=1.5$)")
    save(fig, "heating")

    fig, axes = plt.subplots(1, 2, figsize=(6.65, 3.2), layout="constrained")
    for shape, colour in COLOURS.items():
        model = result["models"][shape]
        rows = model["energy"]
        axes[0].plot(RADII, [v["fraction_apocentre_beyond_2r"] for v in rows], "o-", color=colour,
                     label=f"{shape.title()}: apocentre > 2r")
        axes[0].plot(RADII, [v["fraction_apocentre_beyond_rt"] for v in rows], "s--", color=colour,
                     label=f"{shape.title()}: apocentre > tidal scale")
        for field, style in [("late_mean_local_apocentre", "o-"), ("late_mean_extended_apocentre", "s--")]:
            axes[1].plot(RADII, [v[field] for v in rows], style, color=colour,
                         label=f"{shape.title()}: {'apocentre < 2r' if 'local' in field else 'apocentre > 2r'}")
    axes[0].set(ylabel="Fraction of the initial shell cohort", ylim=(-.02, 1.02))
    axes[1].set(ylabel=r"Late mean energy change [km$^2$ s$^{-2}$]")
    axes[1].set_yscale("symlog", linthresh=.01)
    for ax in axes:
        ax.set(xscale="log", xlabel="Initial shell radius r [pc]")
        ax.set_xticks([3, 10, 20, 35, 70], ["3", "10", "20", "35", "70"])
        ax.legend(fontsize=6.5, loc="upper left")
    save(fig, "orbital_cohorts")
    record(ANALYSIS/"summary.json")
    source_files = ["src/ocen_dm/dynamics/models.py", "src/ocen_dm/dynamics/diagnostics.py",
                    "src/ocen_dm/dynamics/orbits.py", "src/ocen_dm/dynamics/experiment.py",
                    "src/ocen_dm/dynamics/__init__.py", "src/ocen_dm/tails/truncated_equilibrium.py",
                    "tests/test_dynamical_analytics.py"]
    result.update(source_sha256=digest(Path(__file__)), inputs_sha256=hashes,
                  supporting_code_sha256={p: digest(ROOT/p) for p in source_files},
                  figures_sha256={str(p.relative_to(ROOT)): digest(p) for p in figures})
    (data_output/"dynamical_analytics.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
