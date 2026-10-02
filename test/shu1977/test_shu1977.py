"""Shu (1977) singular isothermal sphere collapse: checks that exactly one sink forms, compares
spurious COM drift with and without RANDOMIZE_GRAVTREE, and runs the collapse once under
TREE_INTEGRITY_AUDITS.

The sphere is symmetric and at rest, so its net momentum should stay zero by symmetry and any
drift is spurious. Non-periodic setup, so this exercises the move/enlarge-root-node path.

The momentum accounting must include the sink (PartType5), not just the gas: this run accretes
gas into a sink, and gas alone would report that physical transfer as a huge false drift.
"""

import os
import sys

import pytest
import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
from os import path
from meshoid import Meshoid
from gizmo.test import (build_and_run_test, flush_colorbar, assert_final_time,
                        get_final_snapshot, default_omp_threads, clean_test_outputs,
                        variant_output_dir)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from momentum_drift_common import (  # noqa: E402
    DRIFT_SANITY_CEILING, assert_randomized_drift, mass_bookkeeping_error, measure_and_record,
    plot_momentum_drift, report_momentum_drift,
)

# gas + sink: both carry momentum, and mass moves between them during the run
MOMENTUM_PARTTYPES = ("PartType0", "PartType5")


def plot_shu1977_density_slice(coords, rho, boxsize, output_dir="."):
    """Plot a density slice through the Shu 1977 collapse center."""
    center = np.average(coords, axis=0)
    size = boxsize * 0.1
    M = Meshoid(coords)
    rho_slice = M.Slice(np.log10(rho), res=1024, plane="z", center=center, size=size, order=0)
    fig, ax = plt.subplots(figsize=(6, 6))
    half = size / 2
    im = ax.imshow(rho_slice.T, origin="lower", cmap="inferno", extent=[-half, half, -half, half])
    flush_colorbar(im, ax=ax, label="log10(Density)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Shu 1977 Collapse - Density Slice")
    fig.savefig(output_dir + "/Density_2D.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


@pytest.mark.parametrize("num_mpi_ranks", (8,))
@pytest.mark.parametrize("num_omp_threads", (default_omp_threads(),))
@pytest.mark.parametrize(
    "extra_config_flags",
    # tree_audits: fatal tree audits through sink formation and accretion. Every swallow rearranges
    # the particle list under Hermite sinks with no spawning, the case MAINTAIN_TREE_IN_REARRANGE
    # must cover; without it the first post-swallow moments refresh aborts the run.
    [(), ("RANDOMIZE_GRAVTREE",), ("TREE_INTEGRITY_AUDITS",)],
    ids=["baseline", "randomize", "tree_audits"],
)
def test_shu1977(num_mpi_ranks, num_omp_threads, extra_config_flags, request):
    test_name = "shu1977"
    test_dir = f"test/{test_name}"
    variant_id = request.node.callspec.id.split("-")[0]

    clean_test_outputs(test_name, extra_config_flags)
    build_and_run_test(test_name, num_mpi_ranks, num_omp_threads, extra_config_flags)

    final_snap = get_final_snapshot(test_name, extra_config_flags)
    assert_final_time(final_snap, test_name)

    with h5py.File(final_snap, "r") as f:
        num_sinks = f["Header"].attrs["NumPart_ThisFile"][5]

    assert num_sinks == 1, f"[{variant_id}] Expected exactly 1 PartType5 particle, got {num_sinks}"

    # --- spurious COM drift from correlated tree-force errors (RANDOMIZE_GRAVTREE) ---
    traj = measure_and_record(test_dir, variant_id,
                              variant_output_dir(test_name, extra_config_flags),
                              parttype=MOMENTUM_PARTTYPES)
    if traj is not None:
        # If mass is not conserved across the types we summed, the drift is measuring
        # bookkeeping rather than force error and the comparison below is meaningless.
        mass_err = mass_bookkeeping_error(traj)
        assert mass_err < 1e-6, (
            f"[{variant_id}] total mass over {MOMENTUM_PARTTYPES} varies by {mass_err:.3e} -- "
            f"some momentum-carrying particles are unaccounted for, so the COM drift is not a "
            f"clean measure of force error"
        )
        final_drift = report_momentum_drift(test_dir, test_name, variant_id, traj)
        plot_momentum_drift(test_dir, test_name)
        assert final_drift < DRIFT_SANITY_CEILING, (
            f"[{variant_id}] spurious COM velocity reached {final_drift:.3e} of the internal "
            f"velocity dispersion (sanity ceiling {DRIFT_SANITY_CEILING})"
        )
        assert_randomized_drift(test_name, variant_id, final_drift)

    with h5py.File(final_snap, "r") as F:
        coords = F["PartType0/Coordinates"][:]
        rho = F["PartType0/Density"][:]
        boxsize = F["Header"].attrs["BoxSize"]
    plot_shu1977_density_slice(coords, rho, boxsize, output_dir=test_dir)
