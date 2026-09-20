**The premise that the nine histories are ready to drive production simulations is contradicted by the code.** The fast integrators use the wrong acceleration conversion, and the class-3 calculation labels AGAMA time units as Gyr. These need resolving before interpreting differences between orbital classes. [progenitor_orbits.py:234](< /Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:234>), [progenitor_orbits.py:255](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:255>), [friction_test.py:43](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/friction_test.py:43>), [bar_migration.py:39](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:39>), [bar_migration.py:191](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:191>).

My verdict: **agree with gyrfalcON + AGAMA for pilots; conditionally accept a rigid nucleus; reject the proposed sampling, friction prescription and convergence programme as sufficient for production.**

1. **Orbit correctness comes first.**

   With \(C=1.02271\), velocities in km/s and time in Gyr require
   \[
   \dot x=Cv,\qquad \dot v=C\,a_{\rm[km^2\,s^{-2}\,kpc^{-1}]}.
   \]
   The fast code multiplies the drift by \(C\), but divides gravitational and friction accelerations by \(C\). Its kicks are therefore \(C^{-2}=0.9561\) of the correct values: **4.39% too weak**, independently of timestep refinement. [progenitor_orbits.py:230](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:230>), [progenitor_orbits.py:254](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:254>), [progenitor_orbits.py:265](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:265>).

   I queried the installed AGAMA: its configured time unit is **977.792 Myr**. Thus the class-3 integration from 8 to 0 spans **7.822 Gyr**, although the history declares Gyr. Correct the clock and bar-history convention together; merely relabelling one array would introduce another mismatch. [bar_migration.py:34](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:34>), [bar_migration.py:82](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:82>), [bar_migration.py:185](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:185>).

   There is also a specification mismatch: class 1 comprises **three different Galactic potentials**, whereas the IC proposal assigns McMillan17 to classes 1 and 2. Replaying a trajectory in the wrong potential does not preserve its tidal history. [progenitor_orbits.py:367](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:367>), [NBODY_IC_PROPOSAL.md:35](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:35>).

2. **Code choice: gyrfalcON first, but the runtime estimate is not established.**

   gyrfalcON supplies the required individual softenings and block steps. For nine independent histories, running several serial simulations concurrently is a reasonable use of ten P-cores. Start with measured concurrency, rather than assuming ten processes retain single-process performance. [gyrfalcON.1:105](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/man/man1/gyrfalcON.1:105>), [gyrfalcON.1:167](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/man/man1/gyrfalcON.1:167>).

   **pyfalcon** is useful for prototyping custom coupling, but its supplied example is a global-step leapfrog; production block stepping and velocity-dependent drag integration become your responsibility. At the assessment’s assumed 0.7 seconds per full force evaluation, \(10\,{\rm Gyr}/0.02\,{\rm Myr}\) costs **97 hours**, before host forces and output. Python overhead is not the principal issue. [pyfalcon documentation](https://github.com/GalacticDynamics-Oxford/pyfalcon), [NBODY_CODE_ASSESSMENT.md:15](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_CODE_ASSESSMENT.md:15>).

   **Gadget-4** merits an early representative benchmark if increased resolution is required. Its parallelism and hierarchical integration are relevant, but ten cores do not imply a tenfold speedup. Its spline softening is specified as *Plummer-equivalent*; numerical epsilon values cannot be transferred blindly between engines. [Gadget-4 documentation](https://wwwmpa.mpa-garching.mpg.de/gadget4/05_parameterfile/).

   **Restricted N-body** is suitable for screening and controls, not the final remnant measurement. However, describing AGAMA’s method as simply lacking self-gravity is imprecise: its example periodically updates the satellite potential from the particles using a spherical monopole. It misses nonspherical collective response. [AGAMA example](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_nbody_simulation.py).

   The proposed cost arithmetic gives **1.94 + 2.78 = 4.72 hours**. The unsupported input is the **0.02-second partial step**, including tree maintenance and external-force evaluation. The reported bar benchmark uses only \(10^4\) particles. Moreover, DM passing through the nucleus also needs short steps; this requirement survives replacing the live nucleus. Treat **5–8 hours as an unvalidated estimate**. [NBODY_CODE_ASSESSMENT.md:16](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_CODE_ASSESSMENT.md:16>), [NBODY_CODE_ASSESSMENT.md:23](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_CODE_ASSESSMENT.md:23>), [NBODY_CODE_ASSESSMENT.md:37](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_CODE_ASSESSMENT.md:37>).

3. **A rigid nucleus is defensible as an approximation; “live is wrong by construction” is not.**

   Re-evaluating the proposal’s formula, \(0.1N\,t_{\rm cross}/\ln N\), with \(t_{\rm cross}=0.3\) Myr gives:

   | Nucleus particles | Relaxation estimate |
   |---:|---:|
   | \(2\times10^4\) | 60.6 Myr |
   | \(10^5\) | 261 Myr |
   | \(10^6\) | 2.17 Gyr |

   Thus 60 Myr is correct; 3 Gyr is a rough overestimate. These are **unsoftened estimates**, not demonstrations of evaporation within a Gyr. Softening changes the minimum impact parameter and Coulomb logarithm. For illustration, \(\ln(7/0.5)=2.64\), substantially below \(\ln(20000)=9.90\). It can suppress relaxation appreciably, at the cost of force bias. [NBODY_IC_PROPOSAL.md:10](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:10>).

   For a true Plummer potential with scale \(a=7\) pc, the force differs from a point mass by **2.87% at 50 pc** and **0.90% at 90 pc**, but **45% at 10 pc**. Outer DM particles can have inner pericentres, so their present radius alone does not justify the approximation.

   Also, \(r_{\rm half,3D}=1.305a\): if 7 pc means the three-dimensional half-mass radius, use **\(a=5.36\) pc**. If it means projected half-light radius, \(a=7\) pc is consistent with a Plummer model. The proposal conflates these definitions. [NBODY_IC_PROPOSAL.md:23](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:23>).

   **The implementation trap is more serious:** falcON defaults to kernel P1, whereas P0 is Plummer; its individual pair softening is \((\epsilon_i+\epsilon_j)/2\). A nucleus assigned 7 pc interacting with a 50 pc particle therefore does not generate a 7 pc Plummer field. Prefer an explicitly calibrated rigid potential with consistent centre dynamics, or verify the actual particle-pair force before constructing the DF. [gyrfalcON.1:177](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/man/man1/gyrfalcON.1:177>), [forces_C.h:126](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/inc/forces_C.h:126>), [gravity.h:185](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/inc/public/gravity.h:185>), [kernel.h:111](</Users/vasilybelokurov/Work/src/nemo_satellite_experiment/usr/dehnen/falcON/inc/public/kernel.h:111>).

   The literature justification is partly wrong: **Boldrini et al. (2020) used live GCs with roughly \(10^4\) stellar particles, equal particle masses across components, and softening tests**, not single massive GC particles. Their result does not establish convergence here, but contradicts the citation’s claimed support. [NBODY_IC_PROPOSAL.md:25](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:25>), [Boldrini et al., sections 3 and Appendix A](https://academic.oup.com/mnras/article/492/3/3169/5700296).

4. **Reject radius-only thinning and universal mass-ratio rules.**

   Inverse-probability weighting preserves the expected initial mass distribution. It does **not** stop a heavy particle initially near apocentre from traversing the nucleus immediately. Subsequent tides change its orbit; numerical mass segregation can additionally drive heavy particles inward. The proposed radial thinning is therefore insufficient. [NBODY_IC_PROPOSAL.md:28](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:28>).

   Zemp’s method includes **orbital-pericentre refinement**, not just shells. The paper generally recommends neighbouring mass ratios around 2 and a protective refinement radius at least ten times the inner shell radius. For 0.3 kpc, that suggests testing a several-kpc protection region—not declaring 0.3 kpc sufficient. External tides still require contamination monitoring. [Zemp et al. (2008), section 2.3](https://academic.oup.com/mnras/article/386/3/1543/1061468).

   The quoted heating ratios are arithmetically correct **only for the quoted densities**: 0.00444 and 0.444. But \(2500\,M_\odot\,{\rm pc}^{-3}\) is approximately the **central** density of the proposed \(a=7\) pc Plummer nucleus. At 10 pc it is **155**, giving ratios **0.0715 and 7.15**, assuming the same DM density and equal velocity/Coulomb-log factors. With \(r_{\rm half,3D}=7\) pc, the corresponding ratios are **0.0845 and 8.45**. An NFW mass of \(10^{10}M_\odot\) alone cannot determine its density at 10 pc: concentration and mass/redshift convention are missing. [NBODY_IC_PROPOSAL.md:14](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:14>).

   For dwarf stars, increasing stellar particle mass to satisfy a ratio does not cure heating by massive DM particles. The relevant diffusion depends on DM particle mass, density, velocities and effective impact-parameter limits. Ludlow’s experiments specifically find weak dependence on stellar particle mass. [Ludlow et al. (2021)](https://arxiv.org/abs/2105.03561).

   Nor is “softening ≥ light-particle spacing” sufficient. For a \(10^6M_\odot\) particle at relative speed 20 km/s, \(b_{90}=Gm/v^2\approx10.8\) pc. Softening enough to suppress scattering must also preserve the required central force resolution. Exclude such particles from the scientific region rather than trying to make them harmless there.

5. **One million initial DM particles does not establish remnant resolution.**

   The proposal’s final count assumes every survivor has mass \(10^3M_\odot\). Under that assumption:

   | Bound DM | Survivors | Counting-noise floor |
   |---:|---:|---:|
   | \(10^5M_\odot\) | 100 | 10% |
   | \(10^6M_\odot\) | 1,000 | 3.2% |
   | \(10^7M_\odot\) | 10,000 | 1% |

   These are sampling errors, **not dynamical-convergence guarantees**. For unequal masses use
   \[
   N_{\rm eff}=(\sum m_i)^2/\sum m_i^2.
   \]
   One surviving \(10^5M_\odot\) particle can supply the entire smallest claimed remnant. [NBODY_IC_PROPOSAL.md:28](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:28>), [NBODY_IC_PROPOSAL.md:50](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:50>).

   I would target **at least \(10^3\) effective survivors for an integrated mass**, preferably \(10^4\) for a radial profile. Resolving a \(10^5M_\odot\) remnant then requires inner masses around **100 or 10 \(M_\odot\)** respectively. These are planning targets to validate, not universal thresholds. Artificial stripping depends on force resolution and discreteness as well as survivor count. [van den Bosch & Ogiya (2018)](https://arxiv.org/abs/1801.05427).

6. **A moving coordinate frame does not supply dynamical friction.**

   For inertial position \(x=X(t)+r\),
   \[
   \ddot r=a_{\rm self}+a_{\rm host}(X+r,t)+a_{\rm df,i}-\ddot X.
   \]
   The host centre is **\(-X(t)\)** in these translating coordinates. If the prescribed orbit satisfies \(\ddot X=a_{\rm host}(X)+a_{\rm df,COM}\), merely adding the fictitious acceleration \(-\ddot X\) leaves a residual \(-a_{\rm df,COM}\) at the origin unless physical drag is supplied too. The satellite drifts off the prescribed path. The AGAMA frame example is not, by itself, a friction implementation. [NBODY_IC_PROPOSAL.md:36](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:36>), [AGAMA reference.tex:1317](</Users/vasilybelokurov/Work/venvs/.venv/lib/python3.13/site-packages/agama/doc/reference.tex:1317>).

   A correctly constructed prescribed tidal history is a legitimate **conditional experiment**. It is not self-consistent orbital decay when the simulated bound mass differs from the mass used to obtain the orbit. The adopted exponential histories end at **\(6.74\times10^7\), \(1.01\times10^7\), and \(3.55\times10^6M_\odot\)** respectively; they do not all finish at nucleus mass. [progenitor_orbits.py:288](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:288>), [progenitor_orbits.py:302](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/progenitor_orbits.py:302>).

   Iterate the forward orbit with measured bound mass and size, or couple drag online to the remnant’s bulk motion. Applying Chandrasekhar drag independently to internal particle velocities would spuriously alter internal dynamics. Applying identical drag to escaped debris would also be wrong. Even coupled semi-analytic friction misses the responsive host wake, its memory, host recoil/deformation and background streaming; some debris backreaction can be captured if debris remains self-gravitating.

   Class 3 likewise cannot accept an arbitrary intact \(10^{10}\)–\(10^{11}M_\odot\) dwarf on a test-particle trajectory as a physical history. The orbit document itself conditions bar migration on substantial prior envelope loss. Its quoted threshold needs recomputation after the unit correction. [PROGENITOR_ORBITS.md:173](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/PROGENITOR_ORBITS.md:173>).

7. **AGAMA equilibrium construction: agree, with substantially stronger validation.**

   Construct each stellar/DM DF in the **combined potential**, including the nucleus from the beginning. Check recovered densities, velocity moments, DF positivity and convergence of the potential grid. For spherical models, component-wise Eddington/quasispherical construction may suffice; iterative action-based modelling is not mandatory. The sampled equilibrium must match the **softened force actually evolved**, especially near the nucleus.

   Ten dynamical times tests initial transients, not ten-Gyr secular heating. The proposal does not specify which dynamical time it means. Use short equilibrium checks plus long isolated controls, and define halo concentration, truncation, stellar size, anisotropy and the nucleus’s formation history. Equilibrium alone does not determine whether nuclear growth previously contracted the DM. [NBODY_IC_PROPOSAL.md:33](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:33>).

   **The bound-mass diagnostic needs independent validation.** The reusable scaffold currently substitutes a single Plummer potential based on total candidate bound mass; that is unsuitable for a nucleus + extended stellar envelope + DM remnant. [diagnostics.py:58](</Users/vasilybelokurov/Work/Code/satellite_experiment/src/diagnostics.py:58>).

   There is no conserved Jacobi integral in an eccentric orbit through a decelerating bar. Report iterative self-bound mass, aperture mass, and persistence of membership across subsequent passages separately. Use the actual tidal field to define any instantaneous tidal boundary. The existing radius routine uses \(R[M/(3M_{\rm gal})]^{1/3}\), not the general barred-field calculation. [composite.py:164](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/mass_models/composite.py:164>).

   Do not assume 50–90 pc for class 3. As an illustration—not a Hunter24 calculation—a circular flat-rotation model with \(v_c=200\) km/s gives **31 pc at 0.4 kpc**, versus 79 pc at 1.6 kpc. The proposed scale separation is therefore not established. [NBODY_IC_PROPOSAL.md:47](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:47>).

**Before production, I would require these tests.**

- Correct units; compare identical-potential orbits with an independent adaptive integrator. Verify translating-frame/inertial-frame equivalence, drag signs and bar time/angle interpolation.
- Direct-force checks for every mass/softening pairing, including nucleus–DM; convergence in timestep and force tolerance independently.
- Full-duration isolated controls measuring nucleus structure, dwarf-star heating, DM diffusion and heavy-particle contamination.
- Pilots covering a class-1 orbit, coupled class-2 decay and the strongest class-3 passage; include both early stripping and late remnant evolution.
- At least three mass resolutions, several seeds, independently varied numerical softenings, and rigid-versus-live nucleus comparisons through relevant pericentres. A change in rigid-nucleus scale changes the physical model; it is not purely a numerical convergence test.
- A declared acceptance target—for example, under 10% change in surviving DM mass across the two finest resolutions, with adequate \(N_{\rm eff}\)—and denser diagnostics around pericentre. Fifty-Myr snapshots cannot diagnose shocks lasting roughly a Myr. [NBODY_IC_PROPOSAL.md:39](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/NBODY_IC_PROPOSAL.md:39>), [bar_migration.py:203](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/tails/bar_migration.py:203>).

Existing tests do not provide these assurances: the frictionless comparison accepts 0.1-kpc extrema differences; the barred comparison accepts 5% median energy disagreement and 0.2 in success fraction. Those are screening tests, not validation of the numerical coupling. [test_progenitor_orbits.py:30](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/tests/test_progenitor_orbits.py:30>), [test_friction_test.py:10](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/tests/test_friction_test.py:10>).

For data integrity, preserve particle IDs, component labels, masses and softenings across snapshots/restarts: the scaffold’s current NPZ writer stores only mass, position, velocity and time. Its launcher executes configurable shell text, so configurations must be trusted executable inputs. [nemo_io.py:35](</Users/vasilybelokurov/Work/Code/satellite_experiment/src/nemo_io.py:35>), [run_gyrfalcon.py:35](</Users/vasilybelokurov/Work/Code/satellite_experiment/src/run_gyrfalcon.py:35>), [run_gyrfalcon.py:53](</Users/vasilybelokurov/Work/Code/satellite_experiment/src/run_gyrfalcon.py:53>).

**1. Verified versus assumed.** I read the cited implementations, tests, proposals, local falcON/AGAMA documentation and external papers; independently recalculated the numbers above; and queried AGAMA’s unit setting. I did not run simulations, benchmarks or the test suite, and modified nothing. The heating estimates assume the stated DM density and comparable velocity/logarithmic factors. The tidal-radius example assumes a circular flat-rotation model. The runtime assessment accepts the documented timings without independently reproducing them.

**2. Strongest argument against my answer.** A rigid nucleus may suppress precisely the stellar scattering and time-dependent central response that determine survival of tightly bound DM over ten Gyr. Consequently, even excellent collisionless convergence could yield a systematically excessive remnant. My conditional acceptance requires showing that the relevant DM orbits spend too little time in the nucleus for those effects to control the answer.

**3. What remains necessary to be sure.** Determine whether 7 pc is projected or three-dimensional; specify the initial halo/stellar profiles and nuclear growth history; measure survivor orbital pericentres and physical stellar-scattering times; recompute the orbital histories; demonstrate force/DF consistency and uncontaminated central sampling; validate the bound-mass estimator; and measure actual full-resolution throughput, memory and restart fidelity.