"""Plummer cluster of equal-mass circular binaries — integration / preservation test.

Each Plummer particle is replaced by a binary (2 Type-5 sinks, 100 AU separation,
circular orbit). Uses the gravity- and integration-relevant settings from
SINGLE_STAR_STARFORGE_DEFAULTS: HERMITE_INTEGRATION, GRAVITY_ACCURATE_FEWBODY_INTEGRATION,
SINGLE_STAR_TIMESTEPPING, etc.

The Lagrange-radii / density-profile preservation is measured on the raw stars
(consecutive ParticleID pairs 2k-1, 2k define binary k in the IC).

Two variants: the default tree gravity, and SINGLE_STAR_DIRECT_GRAVITY, which replaces every
star-star pair with an exact brute-force sum. With no gas in this problem the second is a pure
direct N-body integration, so the two together measure what the tree approximation costs here.
"""

from os import path
import glob

import h5py
import numpy as np
import pytest
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
    default_mpi_ranks,
    default_omp_threads,
    variant_output_dir,
    assert_final_time,
    parse_params,
)

TEST_NAME = "plummer_binaries"
TEST_DIR = f"test/{TEST_NAME}"
IC_FILE = f"{TEST_DIR}/{TEST_NAME}_ics.hdf5"


def _physical_cpu_count():
    """Number of physical CPU cores on this host (counts hyperthreads as 1).

    On Linux we parse /proc/cpuinfo for unique (physical id, core id) pairs;
    elsewhere we fall back to os.cpu_count() // 2 assuming 2-way SMT.
    """
    import os
    try:
        with open("/proc/cpuinfo") as f:
            text = f.read()
        cores, cur = set(), {}
        for line in text.splitlines():
            if not line.strip():
                if "physical id" in cur and "core id" in cur:
                    cores.add((cur["physical id"], cur["core id"]))
                cur = {}
            elif ":" in line:
                k, v = line.split(":", 1)
                cur[k.strip()] = v.strip()
        if cores:
            return len(cores)
    except OSError:
        pass
    return max(1, (os.cpu_count() or 4) // 2)


# Benchmarked optimum on a 16-physical-core node: 2 MPI ranks x 8 OMP threads. Capped in total
# rather than scaled with the node: there are only 2*N_BINARIES particles and the deep timestep
# hierarchy leaves a handful of them active per step, so wider parallelism adds synchronisation
# without adding work.
PB_MAX_CORES = 8
PB_NUM_MPI_RANKS = 2
PB_NUM_OMP_THREADS = max(1, min(PB_MAX_CORES, _physical_cpu_count()) // PB_NUM_MPI_RANKS)

# KDK resolves the same hard binaries without the Hermite integrator, taking ~2.4x as many steps
# per unit time, and is an xfail on energy conservation either way. A tenth of the run is enough
# to measure how fast it loses the binaries.
KDK_TIME_FRACTION = 0.1

# BOX_PERIODIC is not in this test's Config.sh, and the cluster (~1 pc) is tiny against the
# 300 pc box regardless, so pair separations never wrap.
PERIODIC = False

# Cluster parameters (in code units: pc - km/s - Msun)
SCALE_RADIUS = 1.0           # pc
M_STAR = 1.0                 # Msun (per star)
N_BINARIES = 256
M_CLUSTER = 2 * N_BINARIES * M_STAR
BINARY_SEPARATION_AU = 1000.0
BOXSIZE = 300.0

LAGRANGE_FRACTIONS = (0.5,)
LAGRANGE_TOL = (0.15,)
LAGRANGE_PERCENTILES_DENSE = np.arange(1, 100, dtype=float)


def _ensure_ic():
    if path.isfile(IC_FILE):
        with h5py.File(IC_FILE, "r") as F:
            n_ok = int(F["Header"].attrs["NumPart_Total"][5]) == 2 * N_BINARIES
            # Check the SEPARATION too, not just the count. Running the generator with its own
            # defaults produces a file with the right particle count and a different problem --
            # a 10x tighter binary is ~30x shorter in period and turns a 47 min run into hours,
            # while looking valid to a count-only check.
            sep_ok = False
            if n_ok:
                x = F["PartType5/Coordinates"][:]
                sep = np.median(np.linalg.norm(x[0::2] - x[1::2], axis=1)) * AU_PER_PC
                sep_ok = abs(sep / BINARY_SEPARATION_AU - 1.0) < 1e-3
            if n_ok and sep_ok:
                return
    import importlib.util

    test_dir = path.dirname(path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "make_plummer_binaries_ics", path.join(test_dir, "make_plummer_binaries_ics.py")
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.make_plummer_binaries_ics(
        N_binaries=N_BINARIES, m_star=M_STAR, a_cluster=SCALE_RADIUS,
        binary_separation_au=BINARY_SEPARATION_AU, boxsize=BOXSIZE,
        seed=42, outfile=IC_FILE,
    )


def _load_snapshot(snap):
    """Prefer the IO_HERMITE_SYNC datasets: a consistent (r,v) pair at the output time.

    Coordinates/Velocities are NOT such a pair -- positions are drifted to the output time while
    velocities are left at the last kick (see the IO_VEL comment in file_io/io.cc) -- and in a
    cluster on deep individual timesteps each star was last kicked at a different instant. Any
    quantity combining r and v, or summing v across stars, is then evaluated on a configuration
    the system never occupied.
    """
    with h5py.File(snap, "r") as F:
        g = F["PartType5"]
        synced = "HermiteSyncCoordinates" in g and "HermiteSyncVelocities" in g
        pos = g["HermiteSyncCoordinates"][:] if synced else g["Coordinates"][:]
        vel = g["HermiteSyncVelocities"][:] if synced else g["Velocities"][:]
        mass = g["Masses"][:]
        pot = g["Potential"][:]
        ids = g["ParticleIDs"][:]
        boxsize = float(F["Header"].attrs["BoxSize"])
    # Sort by ID so binary pairs are contiguous (IDs 1,2 are binary 0; 3,4 are binary 1; ...).
    order = np.argsort(ids)
    return pos[order], vel[order], mass[order], pot[order], boxsize


def _wrap(dx, box):
    return dx - box * np.round(dx / box)


def _radii_from_center(pos, boxsize, periodic):
    dx = pos
    if periodic:
        dx = _wrap(dx, boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def _radii_from_com(pos, mass):
    """Radii in the cluster's mass-weighted center-of-mass frame.

    The cluster's COM drifts noticeably over the run (~0.7 pc by t=32.4) due to
    asymmetric ejections and numerical drift. Measuring radii from the fixed
    origin would conflate that bulk-translation with real cluster evolution
    and would inflate Lagrange-radii drift even when the cluster is intact.
    """
    com = np.average(pos, axis=0, weights=mass)
    return np.sqrt(np.sum((pos - com) ** 2, axis=1))


def _potential_energy(pos, mass, boxsize, periodic, softening=0.0):
    """Direct pairwise potential energy from the positions given.

    Not the snapshot's Potential field: that is evaluated by the tree at the DRIFTED positions,
    so pairing it with last-kick velocities gives a kinetic and a potential term belonging to
    different times -- the two halves of the energy disagree, and the resulting "conservation"
    number carries a sampling artifact set by where in each star's timestep the output landed.
    At N=512 the exact O(N^2) sum is milliseconds, so there is no reason to approximate it.

    Absolute values differ slightly from the tree's by its opening error (ErrTolTheta); only
    energy DRIFT is comparable across that change, which is what the assertions use.
    """
    d = pos[:, None, :] - pos[None, :, :]
    if periodic:
        d -= boxsize * np.round(d / boxsize)
    r2 = np.sum(d * d, axis=-1) + softening * softening
    np.fill_diagonal(r2, np.inf)                       # drop self-pairs
    mm = mass[:, None] * mass[None, :]
    return -0.5 * G_CODE * np.sum(mm / np.sqrt(r2))    # 0.5 for double counting


def _total_energy(vel, mass, pot=None, pos=None, boxsize=None, periodic=False):
    """Total energy. With pos given, the potential is recomputed from those positions so that it
    shares a clock with vel; otherwise it falls back to the snapshot's tree-evaluated field."""
    ke = 0.5 * np.sum(mass * np.sum(vel**2, axis=1))
    pe = (_potential_energy(pos, mass, boxsize, periodic) if pos is not None
          else 0.5 * np.sum(mass * pot))
    return ke + pe, ke, pe


def _conservation_trajectories(snap_paths):
    """Energy and total-momentum histories over a run's snapshots.

    Momentum is here because it is the one conserved quantity that resolves the difference between
    tree and direct star-star gravity. Total momentum is conserved iff the pairwise forces are
    antisymmetric, which direct summation makes exact and the tree only approximates -- measured on
    these ICs as sum(m*a) of 6.5e-17 vs 2.1e-6 relative at t=0. Energy cannot see it: both variants
    sit on the same floor (median |dE|/KE0 of 4.8e-4 vs 4.4e-4), set by the integrator and by
    snapshot velocities being half-kicked, not by the force error.

    Returned as a drift in the cluster's bulk velocity, |P(t) - P(0)| / M_total, so it carries units
    of velocity (km/s here) and reads as "spurious bulk motion the cluster acquired".
    """
    times, energies, momenta, ke0 = [], [], [], None
    for s in snap_paths:
        pos, vel, mass, pot, boxsize = _load_snapshot(s)
        with h5py.File(s, "r") as F:
            t = float(F["Header"].attrs["Time"])
        e, ke, _ = _total_energy(vel, mass, pos=pos, boxsize=boxsize, periodic=PERIODIC)
        times.append(t)
        energies.append(e)
        momenta.append(np.sum(mass[:, None] * vel, axis=0))
        if ke0 is None:
            ke0 = ke
            m_total = mass.sum()
    momenta = np.array(momenta)
    com_vel_drift = np.linalg.norm(momenta - momenta[0], axis=1) / m_total
    return np.array(times), np.array(energies), ke0, com_vel_drift


def _radial_density_profile(r, mass, rbins):
    counts_mass, _ = np.histogram(r, bins=rbins, weights=mass)
    vol = (4.0 / 3.0) * np.pi * (rbins[1:] ** 3 - rbins[:-1] ** 3)
    return counts_mass / vol


def _lagrange_radii(r, mass, fractions=LAGRANGE_FRACTIONS):
    order = np.argsort(r)
    cum = np.cumsum(mass[order]) / mass.sum()
    r_sorted = r[order]
    return np.array([r_sorted[np.searchsorted(cum, f)] for f in fractions])


def _summary_npz_path(variant_id):
    return f"{TEST_DIR}/summary_{variant_id}.npz"


def _plot_signed_log(ax, x, y, color=None, label=None):
    """Plot |y| with positive y solid, negative y dashed. Each contiguous same-sign
    run is plotted as its own Line2D so a transition from positive to negative never
    visually connects through the zero-crossing."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = (x > 0) & (y != 0)
    if not valid.any():
        return
    state = np.where(valid, (y > 0).astype(int), -1)
    bounds = np.flatnonzero(np.diff(state) != 0) + 1
    starts = np.concatenate(([0], bounds))
    ends = np.concatenate((bounds, [len(x)]))
    line = None
    labeled = False
    for s, e in zip(starts, ends):
        if state[s] < 0:
            continue
        seg_x = x[s:e]
        is_pos = bool(state[s])
        seg_y = y[s:e] if is_pos else -y[s:e]
        kwargs = dict(color=color if line is None else line.get_color())
        if not labeled and label is not None:
            kwargs["label"] = label
            labeled = True
        if len(seg_x) == 1:
            ln, = ax.plot(seg_x, seg_y, marker=("o" if is_pos else "x"),
                          linestyle="", **kwargs)
        else:
            ln, = ax.plot(seg_x, seg_y, "-" if is_pos else "--", **kwargs)
        if line is None:
            line = ln


def _plot_summary():
    files = sorted(glob.glob(f"{TEST_DIR}/summary_*.npz"))
    if not files:
        return
    fig, (ax_rho, ax_lag, ax_e, ax_p) = plt.subplots(4, 1, figsize=(7, 16))
    init_done = False
    init_lag_dense = None
    # com_vel_drift is deliberately absent from 'required': summaries written before it existed
    # should still contribute to the other three panels rather than being dropped entirely.
    required = {"variant_id", "rc", "rho_initial", "rho_final",
                "r_lag_initial_dense", "r_lag_final_dense",
                "times", "energies", "ke0"}
    for f in files:
        data = np.load(f)
        if not required.issubset(set(data.files)):
            continue
        vid = str(data["variant_id"])
        rc = data["rc"]
        if not init_done:
            ax_rho.loglog(rc, data["rho_initial"], "k-", lw=2, label="initial")
            init_lag_dense = data["r_lag_initial_dense"]
            init_done = True
        line, = ax_rho.loglog(rc, data["rho_final"], "--", label=vid)
        rel_lag = (data["r_lag_final_dense"] - init_lag_dense) / init_lag_dense
        _plot_signed_log(ax_lag, LAGRANGE_PERCENTILES_DENSE, rel_lag,
                         color=line.get_color(), label=vid)
        times = data["times"]
        energies = data["energies"]
        ke0 = float(data["ke0"])
        rel_e = (energies - energies[0]) / abs(ke0)
        _plot_signed_log(ax_e, times, rel_e, color=line.get_color(), label=vid)
        if "com_vel_drift" in data.files:
            dP = data["com_vel_drift"]
            good = (times > 0) & (dP > 0)
            if good.any():
                ax_p.loglog(times[good], dP[good], "-", color=line.get_color(), label=vid)
    ax_rho.set_xlabel("r [pc]")
    ax_rho.set_ylabel(r"$\rho(r)$ [Msun/pc$^3$]")
    ax_rho.set_title(f"{TEST_NAME}: final density profile (initial in black)")
    ax_rho.legend(fontsize=8)

    ax_lag.set_yscale("log")
    ax_lag.set_xlabel("Lagrange percentile")
    ax_lag.set_ylabel(r"$|r_{\rm final} - r_{\rm initial}| / r_{\rm initial}$  (solid: $>0$, dashed: $<0$)")
    ax_lag.legend(fontsize=8)

    ax_e.set_xscale("log")
    ax_e.set_yscale("log")
    ax_e.set_xlabel("t [code time units]")
    ax_e.set_ylabel(r"$|E(t) - E(0)| / |KE_0|$  (solid: $>0$, dashed: $<0$)")
    ax_e.legend(fontsize=8)

    # Momentum is conserved iff the pairwise forces are antisymmetric, so this panel -- unlike the
    # energy one above, where both variants sit on the same integrator/IO floor -- is where exact
    # star-star summation separates from the tree.
    ax_p.set_xscale("log")
    ax_p.set_yscale("log")
    ax_p.set_xlabel("t [code time units]")
    ax_p.set_ylabel(r"$|P(t) - P(0)| / M_{\rm tot}$  [km/s]")
    ax_p.set_title("spurious bulk velocity (force antisymmetry)", fontsize=9)
    ax_p.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{TEST_DIR}/{TEST_NAME}_summary.png", dpi=120)
    plt.close()


def _plot_variant_density_evolution(variant_id, snaps):
    """One panel per variant: density profile at each snapshot, colored by time."""
    rbins = np.geomspace(SCALE_RADIUS / 4, 5 * SCALE_RADIUS, 30)
    rc = np.sqrt(rbins[:-1] * rbins[1:])
    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.get_cmap("viridis")
    n = len(snaps)
    for i, s in enumerate(snaps):
        pos, _, mass, _, _ = _load_snapshot(s)
        with h5py.File(s, "r") as F:
            t = float(F["Header"].attrs["Time"])
        r = _radii_from_com(pos, mass)
        rho = _radial_density_profile(r, mass, rbins)
        ax.loglog(rc, rho, "-", color=cmap(i / max(n - 1, 1)), label=f"t = {t:.2f}")
    ax.set_xlabel("r [pc]")
    ax.set_ylabel(r"$\rho(r)$ [Msun/pc$^3$]")
    ax.set_title(f"{TEST_NAME} [{variant_id}]: density profile per snapshot")
    ax.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(f"{TEST_DIR}/{TEST_NAME}_{variant_id}_density.png", dpi=120)
    plt.close()


@pytest.mark.parametrize("num_mpi_ranks", (PB_NUM_MPI_RANKS,))
@pytest.mark.parametrize("num_omp_threads", (PB_NUM_OMP_THREADS,))
@pytest.mark.parametrize("extra_config_flags", [
    pytest.param((), id="starforge_defaults"),
    # Every particle here is a star, so this variant routes 100% of the gravity through the
    # brute-force pass and none through the tree -- the tree's opening error is removed entirely
    # rather than reduced. That makes it the sharpest available check on SINGLE_STAR_DIRECT_GRAVITY,
    # and the pair with starforge_defaults isolates what the tree approximation costs a cluster of
    # hard binaries. Compatible only because SINGLE_STAR_STARFORGE_DEFAULTS sets
    # SINGLE_STAR_TIMESTEPPING=0, which leaves SINGLE_STAR_FIND_BINARIES off; the two are a
    # compile-time #error together (see precompiler_logic.h).
    pytest.param(("SINGLE_STAR_DIRECT_GRAVITY",), id="direct_gravity"),
    # KETJU/MSTAR integrates each binary as a regularized chain; the tree still sees every star.
    pytest.param(("KETJU_REGULARIZATION",), id="ketju"),
])
def test_plummer_binaries(num_mpi_ranks, num_omp_threads, extra_config_flags, request):
    _ensure_ic()

    # KETJU variant: one crossing time (3.24) rather than ten, for run time. The coupling
    # conserves energy as well as Hermite here (|dE|/KE0 ~3e-5 at t=3.24, 5.7e-4 at 32.4 vs
    # Hermite's 1.4e-3); the "~1%/crossing residual" the variant was first written around was
    # the coupling bugs fixed on this branch, not a property of the method.
    is_ketju = "KETJU_REGULARIZATION" in extra_config_flags
    ketju_timemax_override = 3.24 if is_ketju else None

    clean_test_outputs(TEST_NAME, extra_config_flags)
    time_max = float(parse_params(f"{TEST_DIR}/{TEST_NAME}.params")["TimeMax"])
    overrides = None
    if is_ketju:
        # KETJU keeps every chain member active every step (~512 force updates per step), so the
        # default TreeDomainUpdateFrequency=0.1 would decompose every step (~80% of wall time).
        time_max = ketju_timemax_override
        overrides = {"TimeMax": ketju_timemax_override, "TreeDomainUpdateFrequency": 10}
    build_and_run_test(TEST_NAME, num_mpi_ranks, num_omp_threads, extra_config_flags,
                       param_overrides=overrides)

    outputdir = variant_output_dir(TEST_NAME, extra_config_flags)
    snaps = sorted(glob.glob(outputdir + "/snapshot_*.hdf5"))
    if len(snaps) < 2:
        raise RuntimeError(f"GIZMO did not produce enough snapshots in {outputdir}")
    assert_final_time(snaps[-1], TEST_NAME, time_max=time_max)

    pos0, vel0, mass0, pot0, boxsize = _load_snapshot(snaps[0])
    posf, velf, massf, potf, _ = _load_snapshot(snaps[-1])
    # Use COM-frame radii so the panels reflect cluster evolution, not
    # COM drift through the box-fixed coordinate frame.
    r0 = _radii_from_com(pos0, mass0)
    rf = _radii_from_com(posf, massf)

    rbins = np.geomspace(SCALE_RADIUS / 4, 5 * SCALE_RADIUS, 20)
    rc = np.sqrt(rbins[:-1] * rbins[1:])
    rho_initial = _radial_density_profile(r0, mass0, rbins)
    rho_final = _radial_density_profile(rf, massf, rbins)
    r_lag_initial = _lagrange_radii(r0, mass0)
    r_lag_final = _lagrange_radii(rf, massf)
    dense_fracs = LAGRANGE_PERCENTILES_DENSE / 100.0
    r_lag_initial_dense = _lagrange_radii(r0, mass0, fractions=dense_fracs)
    r_lag_final_dense = _lagrange_radii(rf, massf, fractions=dense_fracs)
    times, energies, ke0_traj, com_vel_drift = _conservation_trajectories(snaps)

    variant_id = request.node.callspec.id.split("-")[0]
    np.savez(
        _summary_npz_path(variant_id),
        variant_id=variant_id,
        rc=rc,
        rho_initial=rho_initial,
        rho_final=rho_final,
        r_lag_initial=r_lag_initial,
        r_lag_final=r_lag_final,
        r_lag_initial_dense=r_lag_initial_dense,
        r_lag_final_dense=r_lag_final_dense,
        fractions=np.array(LAGRANGE_FRACTIONS),
        fractions_dense=dense_fracs,
        times=times,
        energies=energies,
        ke0=ke0_traj,
        com_vel_drift=com_vel_drift,
    )
    _plot_summary()

    e0, ke0, pe0 = _total_energy(vel0, mass0, pot0)
    ef, _, _ = _total_energy(velf, massf, potf)
    rel_e_err = abs(ef - e0) / abs(ke0)
    # 2% for the KETJU variant is the threshold it was written with; the measured error is
    # ~3e-5, so this catches only gross breakage. Tighten once the chain-member snapshot state
    # (chord positions/velocities mid-step) is written exactly and the metric can use it.
    energy_tol = 0.02 if is_ketju else 0.01
    assert rel_e_err < energy_tol, (
        f"Energy not conserved: |dE|/KE_0 = {rel_e_err:.4f} (>{energy_tol*100:.0f}%)  "
        f"(E0={e0:.4g}, Ef={ef:.4g}, KE0={ke0:.4g}, PE0={pe0:.4g})"
    )

    rel = np.abs(r_lag_final - r_lag_initial) / r_lag_initial
    failures = [
        f"r_{int(f * 100)}: |dr/r|={rel[i]:.3f} (tol {LAGRANGE_TOL[i]}); "
        f"r0={r_lag_initial[i]:.4g}, rf={r_lag_final[i]:.4g}"
        for i, f in enumerate(LAGRANGE_FRACTIONS) if rel[i] > LAGRANGE_TOL[i]
    ]
    assert not failures, "Lagrange radii drifted:\n  " + "\n  ".join(failures)
