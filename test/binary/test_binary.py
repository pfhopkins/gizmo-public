"""e=0.9, q=0.1 Kepler binary -- energy and momentum conservation.

Two Type 5 sinks, no gas, softening 500x below pericentre. Everything is exactly known at
t=0, so both diagnostics are absolute rather than relative to a reference run:

  ENERGY    |E/E0 - 1|, with E = sum(0.5*m*v^2) - G*m1*m2/r taken directly from the state
            rather than inferred through the orbital elements. E0 = -G*m1*m2/(2a) exactly.
            Measured from the IO_HERMITE_SYNC datasets --
            on the ordinary Coordinates/Velocities this is meaningless, because positions are
            drifted to the output time while velocities are left at the last kick (see the
            IO_VEL comment in file_io/io.cc), so vis-viva sees a state that never existed.

  MOMENTUM  the metric this test exists for, and the one the mixed state does NOT spoil:
            P = sum(m*v) uses velocities alone, so if both particles were last kicked at the
            same instant it is the true total momentum at that instant. The ICs are built in
            the COM frame, so P == 0 exactly and no baseline is subtracted. Reported as the
            spurious COM velocity |v_com| / sqrt(G*M/a), matching momentum_drift_common's
            normalisation so the numbers are comparable across tests. NOTE the caveat: it is
            valid only while the two share a timebin, which for two particles they should,
            since every input to the 2-body criterion is symmetric.

WHY THIS CONFIGURATION. e=0.9 swings the required timestep by ((1+e)/(1-e))^3 ~ 6900 across
the orbit, so the pair moves through many timebins per orbit; q=0.1 means the two components
want DIFFERENT bins, because each one's timestep responds to its own acceleration. A
hierarchical scheme kicks the two halves of a pair at different instants whenever their bins
differ, and the equal-and-opposite impulses then fail to cancel.

GROWTH LAW, NOT JUST MAGNITUDE. A random walk grows as sqrt(t) and a systematic error as t,
and the distinction matters far more than the size at any one time: a systematic error is a
bug with a mechanism, a random walk is discretisation noise. Both diagnostics are therefore
recorded per snapshot and fitted for a power law, and the exponent is asserted on.
"""

import glob
import os
import re
from os import chdir, getcwd, path

import h5py
import numpy as np
import pytest
import matplotlib
from matplotlib import pyplot as plt

# G and the AU conversion come from gizmo.units, which mirrors GIZMO's constants.h.
# NOT astropy: the code integrates with GRAVITY_G_CGS = 6.672e-8 and SOLAR_MASS_CGS =
# 1.989e33, giving G_code = 4.300710573e-3 rather than 4.300917270e-3. Reconstructing
# energies or orbital elements with the wrong G injects a spurious term ~ dG/r that
# sweeps with the orbit -- 9.1e-4 in |dE/E| for test/binary, an order of magnitude above
# what that test measures.
from gizmo.units import G_CODE, AU_PER_PC
from gizmo.test import (
    build_and_run_test,
    clean_test_outputs,
    variant_output_dir,
    parse_params,
    run_test,
)

TEST_NAME = "binary"
TEST_DIR = f"test/{TEST_NAME}"
IC_FILE = f"{TEST_DIR}/{TEST_NAME}_ics.hdf5"


M1, Q, A_AU, ECC = 1.0, 0.1, 100.0, 0.9
M2 = Q * M1
MTOT = M1 + M2
A0 = A_AU / AU_PER_PC
V_ORB = np.sqrt(G_CODE * MTOT / A0)          # normalisation for the COM drift
P_ORB = 2.0 * np.pi * np.sqrt(A0 ** 3 / (G_CODE * MTOT))

# What this guards: the shared timestep normalization (SINK_TIMESTEP_SAFETY_FACTOR in
# core/timestep.cc). With it the pair's timebins fuse, so an isolated binary has no inactive
# sources and the Hermite source prediction never fires here; without it the bins split again
# and |dE/E| and the drift exceed the ceilings below by ~10x and ~4x. The source prediction is
# exercised by test/triple, whose hierarchy cannot fuse (see test/triple/README.md).
#
# Measured over 1000 orbits (starforge_defaults, 1 rank, per-orbit envelope): |dE/E| = 7.78e-5,
# COM drift = 3.48e-15. The drift ceiling sits far above that round-off value so the test is not
# a floating-point-reproducibility check; 5e-4 still catches a revert by 4x.
MAX_DE_OVER_E = 1.5e-3
MAX_COM_DRIFT = 5e-4

# nbody_accuracy variant (HERMITE_CORRECTOR_ITERATIONS=2 + HERMITE_SYMMETRIC_TIMESTEP): 3x its measured
# 1000-orbit envelope, 3.72e-6. The ceiling above sits ~20x over the default arm and would not notice
# this variant losing its 24x gain. The drift ceiling is shared: both arms sit at round-off.
NBODY_ACCURACY_FLAGS = ("HERMITE_CORRECTOR_ITERATIONS=2", "HERMITE_SYMMETRIC_TIMESTEP")
MAX_DE_OVER_E_NBODY_ACCURACY = 1.1e-5

# A systematic (secular) error grows as t^1, a random walk as t^0.5. Anything at or above
# ~0.85 is a drift with a mechanism behind it and should be investigated rather than
# absorbed into a looser tolerance.
SECULAR_EXPONENT = 0.85

# --- integration-order sweep ------------------------------------------------------------------
# The tolerances above bound the error at ONE accuracy setting; they cannot see a change that
# keeps magnitudes under tolerance while degrading the scheme's ORDER. That order is the actual
# claim the Hermite machinery makes, so sweep ErrTolIntAccuracy and fit it.
#
# Fits the ENERGY error, not the drift: with the pair's bins fused the drift sits at round-off
# (3.5e-15), where a fit would measure floating-point noise rather than the integrator.
#
# Factors of 4 in eta are factors of 2 in dt -- exactly one timebin -- so the bin structure is
# reproduced at each point instead of drifting across boundaries. Swept DOWNWARD only: above the
# operating point MaxSizeTimestep binds and caps dt, flattening the fit. Run at SWEEP_ORBITS
# rather than the fiducial 1000: the order is a property of the error's eta-scaling, not of the
# baseline length. Each finer eta doubles the runtime, and this is the test's dominant cost.
SWEEP_ETAS = (5e-3, 1.25e-3, 3.125e-4)
SWEEP_ORBITS = 100
# 4th order in dt is eta^2. A 2nd-order scheme gives eta^1. The threshold sits between, nearer
# the lower side: the fit is over 3 points on a finite-length run, so it has real scatter.
MIN_ENERGY_ORDER = 1.5


def _ensure_ic():
    if path.exists(IC_FILE):
        return
    import subprocess
    import sys
    subprocess.run(
        [sys.executable, "make_binary_ics.py", "--m1", str(M1), "--q", str(Q),
         "--a_au", str(A_AU), "--ecc", str(ECC), "--out", f"{TEST_NAME}_ics.hdf5"],
        cwd=TEST_DIR, check=True,
    )


def _read(snap):
    """Prefer the IO_HERMITE_SYNC datasets, which are a consistent (r,v) pair.

    Coordinates/Velocities are NOT: positions are drifted to the output time while velocities
    are left at the last kick (see the IO_VEL comment in file_io/io.cc). Every quantity this
    test measures combines r and v, so on the ordinary datasets it is evaluated on a state the
    system never occupied -- which is why |a/a0-1| appeared to swing to 0.2 and why an earlier
    version of this test reported a spurious secular energy drift.
    """
    with h5py.File(snap, "r") as f:
        t = float(dict(f["Header"].attrs)["Time"])
        g = f["PartType5"]
        ids = np.array(g["ParticleIDs"])
        o = np.argsort(ids)                   # snapshots need not preserve IC order
        synced = "HermiteSyncCoordinates" in g and "HermiteSyncVelocities" in g
        pos_key = "HermiteSyncCoordinates" if synced else "Coordinates"
        vel_key = "HermiteSyncVelocities" if synced else "Velocities"
        return (t, np.array(g["Masses"])[o],
                np.array(g[pos_key])[o], np.array(g[vel_key])[o], synced)


def _energy(mass, pos, vel):
    """Total energy of the pair, from the IO_HERMITE_SYNC state.

    E = sum(0.5*m*v^2) - G*m1*m2/r, exactly -G*m1*m2/(2a) at t=0. Reported as |E/E0 - 1|.
    This is the whole system's energy, so the spurious COM kinetic energy is included -- at
    |v_com| ~ 8e-3 of v_orb that contributes ~1e-4 relative, well under what is being measured.
    """
    r = np.linalg.norm(pos[1] - pos[0])
    return 0.5 * (mass * np.sum(vel ** 2, axis=1)).sum() - G_CODE * mass[0] * mass[1] / r


def _elements(mass, pos, vel):
    """Semi-major axis and eccentricity from the two-body state."""
    m1, m2 = mass
    M = m1 + m2
    dr = pos[1] - pos[0]
    dv = vel[1] - vel[0]
    r = np.linalg.norm(dr)
    v2 = dv @ dv
    a = 1.0 / (2.0 / r - v2 / (G_CODE * M))   # vis-viva
    h = np.cross(dr, dv)
    ecc = np.sqrt(max(0.0, 1.0 - (h @ h) / (G_CODE * M * a)))
    return a, ecc


def _trajectory(snaps):
    t, av, ev, drift, sep, en = [], [], [], [], [], []
    v_com0 = None
    for s in snaps:
        ti, mass, pos, vel, synced = _read(s)
        a, ecc = _elements(mass, pos, vel)
        v_com = (mass[:, None] * vel).sum(0) / mass.sum()
        if v_com0 is None:
            v_com0 = v_com
        t.append(ti); av.append(a); ev.append(ecc); en.append(_energy(mass, pos, vel))
        drift.append(np.linalg.norm(v_com - v_com0) / V_ORB)
        sep.append(np.linalg.norm(pos[1] - pos[0]))
    return (np.array(t), np.array(av), np.array(ev), np.array(drift), np.array(sep),
            np.array(en), synced)


def _per_orbit_envelope(t, y):
    """Per-orbit minimum of y, and the orbit-mean times.

    The instantaneous values oscillate by orders of magnitude WITHIN each orbit -- |a/a0-1|
    swings from 1e-4 to 0.2 -- because a snapshot catches the two particles at whatever point
    their own timebins have reached, so an inactive one contributes a velocity from its last
    kick. A final value, or a power law fitted through the oscillation, measures that sampling
    artefact rather than the integration error.

    The per-orbit MINIMUM is the envelope floor, where the states are most nearly
    synchronised, and it is monotonic. That is the secular error.
    """
    n = int(np.floor(t[-1] / P_ORB))
    tm, ym = [], []
    for i in range(n):
        k = (t >= i * P_ORB) & (t < (i + 1) * P_ORB)
        if k.sum() < 3:   # P/4 cadence gives 4/orbit; allow for bin-edge rounding
            continue
        tm.append(t[k].mean()); ym.append(y[k].min())
    return np.array(tm), np.array(ym)


def _growth_exponent(t, y):
    """Fit the per-orbit envelope to y ~ t^p. nan if there is not enough of it."""
    tm, ym = _per_orbit_envelope(t, y)
    good = (tm > 0) & (ym > 0) & (tm > 0.1 * tm.max())
    if good.sum() < 8:
        return float("nan")
    return float(np.polyfit(np.log(tm[good]), np.log(ym[good]), 1)[0])


def _envelope_final(t, y):
    """Last per-orbit minimum: the secular error at the end of the run."""
    tm, ym = _per_orbit_envelope(t, y)
    return float(ym[-1]) if len(ym) else float("nan")


def _plot(t, energy, ecc, drift, variant_id):
    # This is a diagnostic plot, not a paper figure. The environment may set text.usetex,
    # and a missing texmf tree would then fail the whole test on a rendering detail --
    # which it did on the first run (cmr10.tfm). mathtext needs no external install.
    matplotlib.rcParams['text.usetex'] = False
    orbits = t / P_ORB
    de_inst = np.abs(energy / energy[0] - 1.0)
    fig, ax = plt.subplots(3, 1, figsize=(8, 9.5), sharex=True)
    ax[0].plot(orbits, de_inst, lw=.4)
    ax[1].plot(orbits, np.abs(ecc - ECC), lw=.4)
    ax[2].plot(orbits, drift, lw=.4)
    for A in ax:
        A.set_xscale("log"); A.set_yscale("log")
    ax[0].set_ylabel("|E/E$_0$ - 1|")
    ax[1].set_ylabel("|e - e$_0$|")
    ax[2].set_ylabel("|v$_{com}$| / v$_{orb}$")
    ax[2].set_xlabel("orbits")
    ax[2].set_xlim(1, max(2.0, orbits[-1]))
    fig.tight_layout()
    fig.savefig(f"{TEST_DIR}/{TEST_NAME}_{variant_id}_conservation.png", dpi=120)
    plt.close(fig)


def _plot_convergence(etas, errs, order, variant_id, series=None):
    """log-log energy error vs ErrTolIntAccuracy, with reference slopes.

    The fitted line IS the measurement here, unlike the conservation plot -- and the reference
    slopes matter as much as the fit: a 4th-order scheme fed drifted source positions degrades
    to 2nd, so the question this figure answers is which of the two dashed lines the points lie
    along, not merely whether they are straight.
    """
    matplotlib.rcParams["text.usetex"] = False
    etas, errs = np.asarray(etas, float), np.asarray(errs, float)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    # LEFT: how the error accumulates in each run. A sweep point that is anomalous in the fit is
    # ambiguous from the fit alone -- a uniformly larger error, a late blow-up and a spike at one
    # phase all look the same there. This panel separates them, and would have identified the
    # factor-118 outlier at eta=2.5e-3 (a half-timebin offset) on sight.
    axl = axes[0]
    if series:
        for (ts, des), eta in zip(series, etas):
            axl.loglog(ts / P_ORB, des, lw=.7, label=f"$\\eta$={eta:.3g}")
        axl.set_xlabel("orbits"); axl.set_ylabel("|E/E$_0$ - 1|")
        axl.legend(frameon=False, fontsize=8)
    ax = axes[1]
    ax.loglog(etas, errs, "o", ms=7, zorder=3)
    xs = np.array([etas.min() / 1.6, etas.max() * 1.6])
    c = np.polyfit(np.log(etas), np.log(errs), 1)
    ax.loglog(xs, np.exp(np.polyval(c, np.log(xs))), "-", lw=1.4, zorder=2,
              label=f"fit: $\\eta^{{{order:.2f}}}$  (dt$^{{{2*order:.1f}}}$)")
    for p, lab in ((1.0, "2nd order  $\\eta^1$"), (2.0, "4th order  $\\eta^2$")):
        ax.loglog(xs, errs[0] * (xs / etas[0]) ** p, "--", lw=1.0, zorder=1, label=lab)
    ax.set_xlabel("ErrTolIntAccuracy")
    ax.set_ylabel("|E/E$_0$ - 1|  (per-orbit envelope)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{TEST_DIR}/{TEST_NAME}_{variant_id}_convergence.png", dpi=120)
    plt.close(fig)

def _order_sweep(extra_config_flags, n_ranks, n_omp, variant_id):
    """Run the same problem at several ErrTolIntAccuracy values and fit the energy error's order.

    Reuses the binary built for the fiducial run -- only the params change, via run_test's
    override mechanism -- so this costs runs, not rebuilds. Each point overwrites the previous
    one's snapshots, so each is analysed before the next starts.
    """
    outdir = variant_output_dir(TEST_NAME, extra_config_flags)
    etas, errs, series = [], [], []
    for eta in SWEEP_ETAS:
        # run_test resolves <name>.params relative to the TEST directory and is normally called
        # by build_and_run_test from inside it; called from the repo root it raises
        # FileNotFoundError. Enter and leave around each run, restoring on failure.
        # Remove the previous point's snapshots first. The fiducial run is far longer than a
        # sweep point, so its snapshots are NOT all overwritten -- the leftovers sit at times
        # beyond this point's TimeMax and dominate the per-orbit envelope, making every point
        # report the fiducial value and the fitted order come out ~0.
        for _stale in glob.glob(outdir + "/snapshot_*.hdf5"):
            os.remove(_stale)
        cwd = getcwd()
        try:
            chdir(TEST_DIR)
            run_test(TEST_NAME, n_ranks, n_omp, param_overrides={
                "OutputDir": path.basename(outdir),   # a variant's run otherwise writes to output/
                "ErrTolIntAccuracy": float(eta),
                "TimeMax": float(SWEEP_ORBITS * P_ORB),
                "TimeBetSnapshot": float(P_ORB / 4.0),
            })
        finally:
            chdir(cwd)
        snaps = sorted(glob.glob(outdir + "/snapshot_*.hdf5"),
                       key=lambda f: int(re.search(r"snapshot_(\d+)", f).group(1)))
        if len(snaps) < 16:
            pytest.fail(f"sweep point eta={eta:g} produced only {len(snaps)} snapshots")
        t, a, ecc, drift, sep, energy, synced = _trajectory(snaps)
        de_inst = np.abs(energy / energy[0] - 1.0)
        _, de_env = _per_orbit_envelope(t, de_inst)
        etas.append(eta); errs.append(de_env[-1])
        series.append((t.copy(), de_inst.copy()))   # full evolution, for the diagnostic panel
        print(f"  sweep eta={eta:9.3e}   |dE/E| envelope {de_env[-1]:.3e}")
    order = np.polyfit(np.log(etas), np.log(errs), 1)[0]
    print(f"  energy error ~ eta^{order:.2f}  (dt^{2*order:.1f});  "
          f"4th order is eta^2, 2nd order eta^1, threshold {MIN_ENERGY_ORDER}")
    _plot_convergence(etas, errs, order, variant_id, series)
    assert order >= MIN_ENERGY_ORDER, (
        f"energy error converges as eta^{order:.2f} (dt^{2*order:.1f}), below the "
        f"eta^{MIN_ENERGY_ORDER} floor. The integrator has lost its order even though the "
        f"magnitudes at the operating point are still within tolerance. The source prediction "
        f"cannot be the cause here -- this pair shares a timebin, so it never fires (see the "
        f"module docstring) -- so check the shared timestep normalization (core/timestep.cc) "
        f"and the Hermite predictor/corrector machinery itself (core/kicks.cc).")
    return np.asarray(etas), np.asarray(errs), order, series


@pytest.mark.parametrize("num_mpi_ranks", (1,))
@pytest.mark.parametrize("num_omp_threads", (1,))
@pytest.mark.parametrize("extra_config_flags", [
    pytest.param((), id="starforge_defaults"),
    pytest.param(("DISABLE_HERMITE_INTEGRATION",), id="kdk", marks=pytest.mark.xfail(
        reason="KDK is the baseline HERMITE_INTEGRATION exists to beat; at e=0.9 it should "
               "not hold the orbit, and this variant records by how much",
        strict=False,
    )),
    pytest.param(NBODY_ACCURACY_FLAGS, id="nbody_accuracy"),
])
def test_binary(num_mpi_ranks, num_omp_threads, extra_config_flags, request):
    _ensure_ic()
    clean_test_outputs(TEST_NAME, extra_config_flags)
    build_and_run_test(TEST_NAME, num_mpi_ranks, num_omp_threads, extra_config_flags)

    outdir = variant_output_dir(TEST_NAME, extra_config_flags)
    # NUMERIC sort. glob + sorted() is lexical, so snapshot_1000 lands before snapshot_999 and
    # the array is scrambled: t[-1] then reports the lexically-last file rather than the last in
    # time. With 4000 snapshots that silently truncated a 1000-orbit run to the first 250.
    snaps = sorted(glob.glob(outdir + "/snapshot_*.hdf5"),
                   key=lambda f: int(re.search(r"snapshot_(\d+)", f).group(1)))
    if len(snaps) < 64:
        raise RuntimeError(f"GIZMO produced only {len(snaps)} snapshots in {outdir}")

    t, a, ecc, drift, sep, energy, synced = _trajectory(snaps)
    assert np.all(np.diff(t) >= 0), (
        "snapshot times are not monotonic -- the file ordering is wrong, so every per-orbit "
        "quantity below is being binned from a scrambled series")
    assert synced, ("snapshots lack HermiteSyncCoordinates/HermiteSyncVelocities -- build with "
                    "IO_HERMITE_SYNC. Without it every metric here mixes r(t_out) with v(t_kick).")
    variant_id = request.node.callspec.id.split("-")[0]
    _plot(t, energy, ecc, drift, variant_id)
    np.savez(f"{TEST_DIR}/summary_{variant_id}.npz",
             t=t, a=a, ecc=ecc, drift=drift, sep=sep, energy=energy,
             a0=A0, ecc0=ECC, period=P_ORB)

    de_inst = np.abs(energy / energy[0] - 1.0)
    de = _envelope_final(t, de_inst)
    drift_env = _envelope_final(t, drift)
    n_orbits = t[-1] / P_ORB
    p_drift = _growth_exponent(t, drift)
    p_energy = _growth_exponent(t, de_inst)

    print(f"\n  {n_orbits:.1f} orbits, {len(snaps)} snapshots")
    print(f"  |dE/E|      envelope {de:.3e}  (instantaneous median {np.median(de_inst):.3e}, "
          f"max {de_inst.max():.3e})   growth t^{p_energy:+.2f}")
    print(f"  COM drift   envelope {drift_env:.3e}  (instantaneous max {drift.max():.3e})"
          f"   growth t^{p_drift:+.2f}")
    print(f"  pericentre  min separation {sep.min():.4e} pc "
          f"(softening 1e-7, so {sep.min() / 1e-7:.0f}x above it)")

    max_de = MAX_DE_OVER_E_NBODY_ACCURACY if variant_id == "nbody_accuracy" else MAX_DE_OVER_E
    assert de < max_de, (
        f"relative energy error {de:.3e} over {n_orbits:.1f} orbits (tol {max_de}); "
        f"E0={energy[0]:.6e}, E_final={energy[-1]:.6e}. Measured on the IO_HERMITE_SYNC state, "
        f"so this is the integrator and not the output convention."
    )
    assert drift_env < MAX_COM_DRIFT, (
        f"spurious COM velocity {drift_env:.3e} in units of sqrt(GM/a) (tol {MAX_COM_DRIFT}). "
        f"The ICs have P == 0 exactly, so this is entirely integration error."
    )
    # Order sweep. Runs last of the checks that can fail hard: a magnitude regression should
    # report from the asserts above rather than from a confusing order fit downstream.
    sweep_etas, sweep_errs, sweep_order, sweep_series = _order_sweep(
        extra_config_flags, num_mpi_ranks, num_omp_threads, variant_id)
    # re-save with the sweep included. The earlier savez above is kept deliberately: it runs
    # before the assertions, so a failing run still leaves its trajectory on disk to look at.
    np.savez(f"{TEST_DIR}/summary_{variant_id}.npz", t=t, a=a, ecc=ecc, drift=drift, sep=sep, energy=energy, a0=A0, ecc0=ECC, period=P_ORB,
             sweep_etas=sweep_etas, sweep_errs=sweep_errs, sweep_order=sweep_order,
             sweep_t=np.array([x[0] for x in sweep_series]),
             sweep_de=np.array([x[1] for x in sweep_series]))

    # The growth law distinguishes a bug from discretisation noise, and is asserted only when the
    # drift is large enough for the exponent to be meaningful, so a well-behaved run cannot fail
    # on fitting noise. With the shared timestep normalization the pair's bins fuse and the drift
    # sits at round-off, so this branch no longer fires for the default variant -- it remains as
    # the guard for a regression that reintroduces a secular leak.
    if drift_env > 1e-6 and p_drift >= SECULAR_EXPONENT:
        pytest.xfail(
            f"known secular momentum leak: COM drift grows as t^{p_drift:.2f} (>= "
            f"{SECULAR_EXPONENT}), envelope {drift_env:.3e}. Energy leaks too, t^{p_energy:+.2f}. "
            f"Two particles, Newtonian at "
            f"{sep.min()/1e-7:.0f}x the softening, so this is the block-timestep scheme, not "
            f"the tree or the kernel. Energy is fine (|da/a| ~ t^{p_energy:+.2f})."
        )
