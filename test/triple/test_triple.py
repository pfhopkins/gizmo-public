"""Hierarchical triple -- momentum and energy conservation across a deep timebin hierarchy.

A 0.8+0.2 Msun eccentric 4 AU inner binary + 1 Msun tertiary at 100 AU, mutually inclined and
randomly oriented: period ratio 88, so the tertiary runs ~6 timebins coarser than the inner
stars at every sync point. This is the configuration test/binary cannot reach -- a bound pair
splits by at most ~1 bin -- and the configuration of the production momentum violation (a
hardening triple). Forces on the fine stars are evaluated with the coarse tertiary seen
mid-step, up to 2^6 of their own steps from its last sync.

Two separate splits are in play, deliberately. The 6-bin inner<->outer gap is the deep one. The
unequal inner masses add a 1-bin split INSIDE the pair, so the inner binary -- which carries 86%
of the total energy -- is itself sensitive to the source prediction. With the previous equal
masses it shared a bin by symmetry and was blind to it, which made the energy diagnostic a
6:1 mixture dominated by a component that could not respond to what was under test. See
make_triple_ics.py for the full reasoning and for what that trade costs.

What is asserted: an energy ceiling that bounds gross breakage only. This configuration is the
one the Old*-based source prediction (gravity/forcetree.cc) exists for -- the hierarchy cannot
fuse -- but the committed assertion does not discriminate that fix, and nothing in the committed
suite does; see the README for what did and what restoring a guard would take.

Measured from the IO_HERMITE_SYNC datasets, for the same reason as test/binary: the plain
Coordinates/Velocities are a mixed state (positions drifted, velocities at the last kick).
"""

import glob
import re
from os import path

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
)

TEST_NAME = "triple"
TEST_DIR = f"test/{TEST_NAME}"
IC_FILE = f"{TEST_DIR}/{TEST_NAME}_ics.hdf5"


M1_IN, M2_IN, M3 = 0.8, 0.2, 1.0    # unequal: splits the pair 1 bin, see make_triple_ics.py
A_IN = 4.0 / AU_PER_PC
A_OUT = 100.0 / AU_PER_PC
# Eccentric, inclined, randomly oriented, unequal-mass: a fully general hierarchical triple.
# Nothing here is a special case -- no shared bin by symmetry, no constant separation, no
# coplanar or axis-aligned degeneracy. See make_triple_ics.py for the two constraints that bound
# it (Mardling-Aarseth hierarchy, and staying clear of the Kozai window).
E_IN, E_OUT, I_MUT = 0.5, 0.3, 25.0
ORIENT_SEED = 20250825
MTOT = M1_IN + M2_IN + M3
P_OUT = 2.0 * np.pi * np.sqrt(A_OUT ** 3 / (G_CODE * MTOT))
V_OUT = np.sqrt(G_CODE * MTOT / A_OUT)

# 3x the per-outer-orbit median |dE/E| envelope measured on this configuration: 5.66e-5 over 50
# orbits (0.8+0.2 Msun eccentric inner binary, inclined tertiary, eta=1.25e-3). This bounds gross
# breakage only: the source-prediction defect passed on magnitude and failed only on the COM
# drift's growth exponent, a check removed pending recalibration (see README.md).
MAX_DE_OVER_E = 1.7e-4

# nbody_accuracy variant (HERMITE_CORRECTOR_ITERATIONS=2 + HERMITE_SYMMETRIC_TIMESTEP): 3x its measured
# envelope, 4.10e-5 over the same 50 orbits.
NBODY_ACCURACY_FLAGS = ("HERMITE_CORRECTOR_ITERATIONS=2", "HERMITE_SYMMETRIC_TIMESTEP")
MAX_DE_OVER_E_NBODY_ACCURACY = 1.2e-4

def _ic_matches():
    """Does the IC on disk have the configuration this test is calibrated for?

    Checking existence alone is not enough. The generator's parameters have changed before, and
    a stale file silently runs the OLD configuration under the new tolerances -- the same failure
    that put test/plummer_binaries on a 100 AU IC while its test described 1000 AU. Recover the
    elements and compare, so a parameter change regenerates instead of being ignored.
    """
    if not path.exists(IC_FILE):
        return False
    try:
        with h5py.File(IC_FILE, "r") as f:
            g = f["PartType5"]
            o = np.argsort(np.array(g["ParticleIDs"]))
            m = np.array(g["Masses"])[o]
            x = np.array(g["Coordinates"])[o]
            v = np.array(g["Velocities"])[o]
        if len(m) != 3:
            return False
        m_in_tot = m[0] + m[1]
        r, dv = x[0] - x[1], v[0] - v[1]
        rn = np.linalg.norm(r)
        a = 1.0 / (2.0 / rn - dv @ dv / (G_CODE * m_in_tot))
        h = np.cross(r, dv)
        e = np.linalg.norm(np.cross(dv, h) / (G_CODE * m_in_tot) - r / rn)
        return abs(a / A_IN - 1.0) < 1e-3 and abs(e - E_IN) < 1e-3
    except (OSError, KeyError):
        return False


def _ensure_ic():
    if _ic_matches():
        return
    import subprocess
    import sys
    subprocess.run(
        [sys.executable, "make_triple_ics.py",
         "--m1", str(M1_IN), "--m2", str(M2_IN), "--m3", str(M3),
         "--a_in_au", "4.0", "--a_out_au", "100.0",
         "--e_in", str(E_IN), "--e_out", str(E_OUT), "--i_mut", str(I_MUT),
         "--seed", str(ORIENT_SEED), "--out", f"{TEST_NAME}_ics.hdf5"],
        cwd=TEST_DIR, check=True,
    )


def _read(snap):
    with h5py.File(snap, "r") as f:
        t = float(dict(f["Header"].attrs)["Time"])
        g = f["PartType5"]
        o = np.argsort(np.array(g["ParticleIDs"]))
        synced = "HermiteSyncCoordinates" in g and "HermiteSyncVelocities" in g
        assert synced, "IO_HERMITE_SYNC datasets missing -- the metrics are meaningless on the mixed state"
        return (t, np.array(g["Masses"])[o],
                np.array(g["HermiteSyncCoordinates"])[o],
                np.array(g["HermiteSyncVelocities"])[o])


def _trajectory(snaps):
    t, de, dr, a_in = [], [], [], []
    e0 = v0 = None
    for s in snaps:
        ti, m, x, v = _read(s)
        E = 0.5 * (m * np.sum(v ** 2, axis=1)).sum()
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                E -= G_CODE * m[i] * m[j] / np.linalg.norm(x[j] - x[i])
        vc = (m[:, None] * v).sum(0) / m.sum()
        # inner semi-major axis as an orbit-health check (vis-viva on the inner pair)
        rr = np.linalg.norm(x[1] - x[0])
        vv = np.sum((v[1] - v[0]) ** 2)
        ain = 1.0 / (2.0 / rr - vv / (G_CODE * (M1_IN + M2_IN)))
        if e0 is None:
            e0, v0 = E, vc
        t.append(ti)
        de.append(abs(E / e0 - 1.0))
        dr.append(np.linalg.norm(vc - v0) / V_OUT)
        a_in.append(ain)
    return np.array(t), np.array(de), np.array(dr), np.array(a_in)


def _envelope(t, y):
    """Per-OUTER-orbit MEDIAN.

    Not the minimum. With IO_HERMITE_SYNC the (r,v) pair is already consistent, so the
    within-orbit spread is the inner binary's phase rather than a sampling artifact; taking the
    minimum tracks the dips of an oscillation instead of the error, and reports a falling trend
    where the error is in fact growing.
    """
    orb = (t / P_OUT).astype(int)
    xs, ys = [], []
    for i in range(orb.max() + 1):
        k = orb == i
        if k.sum() < 3:
            continue
        xs.append(t[k].mean())
        ys.append(np.median(y[k]))
    return np.array(xs), np.array(ys)


def _growth_exponent(t, y):
    xs, ys = _envelope(t, y)
    k = (xs > 0) & (ys > 0)
    if k.sum() < 8:
        return np.nan
    h = k & (xs >= 0.5 * xs[k].max())
    return np.polyfit(np.log(xs[h]), np.log(ys[h]), 1)[0]


def _plot(t, de, dr, variant_id):
    matplotlib.rcParams['text.usetex'] = False
    orbits = t / P_OUT
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    ax[0].plot(orbits, de, lw=.4)
    ax[1].plot(orbits, dr, lw=.4)
    for A in ax:
        A.set_xscale("log")
        A.set_yscale("log")
    ax[0].set_ylabel("|E/E$_0$ - 1|")
    ax[1].set_ylabel("|v$_{com}$| / v$_{out}$")
    ax[1].set_xlabel("outer orbits")
    ax[1].set_xlim(1, max(2.0, orbits[-1]))
    fig.tight_layout()
    fig.savefig(f"{TEST_DIR}/{TEST_NAME}_{variant_id}_conservation.png", dpi=120)
    plt.close(fig)


@pytest.mark.parametrize("num_mpi_ranks", (1,))
@pytest.mark.parametrize("num_omp_threads", (1,))
# Hermite only. A KDK variant here spent a 50-orbit run producing numbers nothing asserted on:
# every ceiling below is calibrated for Hermite. The non-Hermite build path of the
# IO_HERMITE_SYNC output block is covered by test/binary's kdk variant.
@pytest.mark.parametrize("extra_config_flags", [
    pytest.param((), id="starforge_defaults"),
    pytest.param(NBODY_ACCURACY_FLAGS, id="nbody_accuracy"),
])
def test_triple(num_mpi_ranks, num_omp_threads, extra_config_flags, request):
    _ensure_ic()
    clean_test_outputs(TEST_NAME, extra_config_flags)
    build_and_run_test(TEST_NAME, num_mpi_ranks, num_omp_threads, extra_config_flags)

    outdir = variant_output_dir(TEST_NAME, extra_config_flags)
    snaps = sorted(glob.glob(outdir + "/snapshot_*.hdf5"),
                   key=lambda f: int(re.search(r"snapshot_(\d+)", f).group(1)))
    assert len(snaps) >= 64, f"only {len(snaps)} snapshots -- the run died early"

    t, de, dr, a_in = _trajectory(snaps)
    assert np.all(np.diff(t) >= 0), "snapshot times are not monotonic -- file ordering is wrong"

    n_orbits = t[-1] / P_OUT
    _, de_env = _envelope(t, de)
    _, dr_env = _envelope(t, dr)
    p_de = _growth_exponent(t, de)
    p_dr = _growth_exponent(t, dr)

    variant_id = request.node.callspec.id.split("-")[0]
    _plot(t, de, dr, variant_id)
    np.savez(f"{TEST_DIR}/summary_{variant_id}.npz",
             t=t, de=de, drift=dr, a_in=a_in, period_out=P_OUT)

    print(f"  {n_orbits:.1f} outer orbits, {len(snaps)} snapshots")
    print(f"  |dE/E|      envelope {de_env[-1]:.3e}   growth t^{p_de:+.2f}")
    print(f"  COM drift   band max {dr_env.max():.3e}   growth t^{p_dr:+.2f}")
    print(f"  inner orbit a/a0 - 1 = {a_in[-1]/a_in[0]-1:+.3e}")

    # One assertion. The COM drift ceiling and growth exponents are reported above but not
    # asserted: they were calibrated on a configuration this IC no longer produces. The
    # convergence sweep was withdrawn; README.md records why and what a valid one would need.
    max_de = MAX_DE_OVER_E_NBODY_ACCURACY if variant_id == "nbody_accuracy" else MAX_DE_OVER_E
    assert de_env[-1] < max_de, (
        f"relative energy error {de_env[-1]:.3e} over {n_orbits:.0f} outer orbits "
        f"(tol {max_de})")
