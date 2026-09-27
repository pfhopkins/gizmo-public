"""fewbody: a suite of cold-collapse 3-10 body problems, as an energy-conservation stress test.

Each problem is an isolated system of 3-10 Salpeter-sampled stars dropped from rest and allowed to
collapse. Cold collapse is the cheapest way to manufacture the thing that actually breaks gravity
solvers: a deep near-radial infall ending in close encounters that span orders of magnitude in
timestep. Repeating it across a suite makes the result a statement about the method rather than
about one lucky orbit.

Each problem is its own GIZMO run, and the runs are executed concurrently, sized to the cores
available. That costs a little orchestration but buys the thing that matters: because a run
contains exactly one system and nothing else, its energy budget is exact and attributable. A
failure names the problem, its N and its mass spectrum, instead of being averaged away among
sixty others sharing a box.

Each problem is followed for the same number of ITS OWN free-fall times, so a 3-body system of
half-solar stars and a 10-body system containing a 40 Msun star are compared at equal dynamical
age rather than equal wall time.

Four variants, a 2x2 over the two things that can spoil conservation, so each is isolated rather
than inferred:

  tree             production solver, individual timesteps -- tree + integration + unequal-dt
  tree_equaldt     FORCE_EQUAL_TIMESTEPS                   -- tree + integration
  direct_gravity   SINGLE_STAR_DIRECT_GRAVITY              -- integration + unequal-dt
  direct_equaldt   both                                    -- integration alone

SINGLE_STAR_DIRECT_GRAVITY sums every star-star pair exactly; with no gas here that removes the
tree from the force calculation entirely. FORCE_EQUAL_TIMESTEPS puts every particle on one
universal step. So a difference down a column is the tree, a difference across a row is the
timestep hierarchy, and the last row measures the integrator with nothing else in the way -- no
single number would tell you which of the three you were looking at.

Three KETJU_REGULARIZATION variants hand the close encounters to MSTAR's regularised chain
integrator instead of Hermite. Every star is a chain centre (KetjuMinStarMass is below the lightest
star) and centres closer than 2R share a chain, so R sets how much of each problem is regularised:

  ketju          R = 0.005 pc, tree        -- most of each system is one chain after the collapse
  ketju_r200au   R = 0.001 pc, tree        -- only close subsystems are chains
  ketju_direct   R = 0.005 pc, direct sum  -- as ketju, without tree force error in the coupling

Chains use the host's star-star softening, so the softened energy reported below is also the
Hamiltonian MSTAR conserves. That energy is valid for chain members because it is only reported at
full synchronization, where every chain has just ended its step and been given its true velocity.

Energy is read from the in-code synced diagnostic (ENERGY_BUDGET_DIAGNOSTIC), not from snapshots.
Snapshots write Velocities at kick-time against drift-time positions, an O(dt/2t_dyn) error per
particle; measured on this suite that is 3-10% with no secular trend, an order of magnitude above
the 1% under test, and recomputing the potential exactly removes only half of it because the
residual lives in the velocities. plummer_binaries gets away with the snapshot estimate because
at N=512 those per-particle errors average down to 5e-4; at N=3-10, where one star can hold a
third of the binding energy, they do not.
"""

import glob
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from os import path

import h5py
import numpy as np
import pytest
from matplotlib import pyplot as plt

from gizmo.test import (
    build_gizmo_for_test,
    clean_test_outputs,
    variant_output_dir,
    variant_suffix,
    parse_params,
)

TEST_NAME = "fewbody"
TEST_DIR = f"test/{TEST_NAME}"

# --- the suite -------------------------------------------------------------------------------
N_PROBLEMS = 48
N_MIN, N_MAX = 3, 10
M_MIN, M_MAX = 0.5, 50.0     # Msun, Salpeter
R_SYS = 0.01                 # pc
BOXSIZE = 1.0                # pc; only has to be big enough to contain one collapsing system
SEED = 42
N_TFF = 20.0                 # how many free-fall times to follow each problem for
N_SNAPS = 20                 # snapshots per problem

# --- parallelism ----------------------------------------------------------------------------
# 2 ranks per problem, pure MPI. Two rather than one because the parts of the solver most likely
# to be wrong are the parallel ones -- the direct-summation star table is assembled with an
# Allgatherv over per-task counts/offsets, and node sink moments are exchanged across top-level
# nodes; on a single rank both degenerate and stop being exercised. A 3-body problem split over
# two ranks also puts one or two particles on a rank, which is a worthwhile edge case in itself.
# No OpenMP (threads=0 leaves OPENMP out of Config.sh rather than compiling it in and running one
# thread): a handful of particles leaves almost nothing active per step, so threads would buy
# synchronisation rather than throughput, and dropping them removes one of the two sources of
# reduction-order nondeterminism.
RANKS_PER_PROBLEM = 2
OMP_THREADS = 0

# Per-problem wall-clock ceiling. A pathological problem should fail loudly on its own rather
# than hang the suite; unbounded, one stuck run would stall the whole pytest session.
PROBLEM_TIMEOUT_S = float(os.environ.get("FEWBODY_PROBLEM_TIMEOUT", "1800"))

# 10%. The measured worst case over all 48 problems is 0.034 with the source-prediction fix and
# 0.038 without it, so this bounds gross breakage only and does not discriminate the fix (the
# per-problem median does, by 2-2.4x, but no stored-baseline comparison is implemented). The
# number of problems over threshold is not asserted: it counts chaotic outcomes.
ENERGY_TOL = 0.10  # >10% is a fail, for both variants

# KETJU parameters, only written for the KETJU variants (KetjuRegionRadius comes from the variant).
# The minimum centre mass is below the suite's 0.5 Msun floor, so every star can anchor a chain.
KETJU_PARAMS = {"KetjuMinStarMass": 0.1, "KetjuMinBHMass": 0, "KetjuUseStarStarSoftening": 1}

PLOT_PATH = f"{TEST_DIR}/{TEST_NAME}_energy.png"
CURVES_PATH = f"{TEST_DIR}/{TEST_NAME}_energy_curves.png"
PAIRS_PATH = f"{TEST_DIR}/{TEST_NAME}_energy_pairs.png"


def _pct(x):
    """Format a fraction as a percentage, safe under text.usetex.

    With usetex on -- which this repo's matplotlibrc sets -- a bare '%' opens a LaTeX comment and
    silently eats the rest of the label, so "tolerance (1%)" renders as "tolerance (1".
    """
    txt = f"{x:.0%}"
    return txt.replace("%", r"\%") if plt.rcParams.get("text.usetex", False) else txt


# legend.frameon is off in this repo's matplotlibrc, which leaves the legend's sample markers
# floating among the data points they are meant to explain. Draw a box for these plots.
LEGEND_KW = dict(fontsize=8, frameon=True, framealpha=0.95, edgecolor="0.4")


def _available_cores():
    """Cores we may actually use: the affinity mask, so a Slurm cpuset is respected."""
    n = int(os.environ.get("FEWBODY_JOBS", "0"))
    if n > 0:
        return n
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def _n_concurrent():
    return max(1, _available_cores() // RANKS_PER_PROBLEM)


def _make_suite():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "make_fewbody_ics", path.join(path.dirname(path.abspath(__file__)), "make_fewbody_ics.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.make_fewbody_suite(
        n_problems=N_PROBLEMS, n_min=N_MIN, n_max=N_MAX, m_min=M_MIN, m_max=M_MAX,
        r_sys=R_SYS, boxsize=BOXSIZE, seed=SEED,
        outdir=path.join(TEST_DIR, "ics"), prefix=TEST_NAME)


def _write_problem_params(base_text, prob, out_rel, extra=None):
    """One params file per problem: same base settings, only the IC-dependent values differ."""
    tmax = N_TFF * prob["t_ff"]
    over = {
        "InitCondFile": path.join("ics", prob["name"] + "_ics"),
        "OutputDir": out_rel,
        "TimeMax": f"{tmax:.8g}",
        "TimeBetSnapshot": f"{tmax / N_SNAPS:.8g}",
        "TimeBetStatistics": f"{tmax / N_SNAPS:.8g}",
        # cap a single step well inside a free-fall time; the integrator's own criteria set the
        # real step, this only stops the very first step from stepping over the collapse
        "MaxSizeTimestep": f"{prob['t_ff'] / 100.0:.8g}",
        "BoxSize": f"{BOXSIZE:.8g}",
    }
    over.update({k: f"{v:.8g}" if isinstance(v, float) else str(v) for k, v in (extra or {}).items()})
    lines, seen = [], set()
    for line in base_text.split("\n"):
        k = line.split()
        if k and k[0] in over:
            lines.append(f"{k[0]:35s}{over[k[0]]}")
            seen.add(k[0])
        else:
            lines.append(line)
    lines += [f"{k:35s}{v}" for k, v in over.items() if k not in seen]
    p = path.join(TEST_DIR, "params", prob["name"] + ".params")
    os.makedirs(path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(f"% generated by test_fewbody.py for problem {prob['index']} "
                f"(N={prob['n']}, M={prob['m_total']:.3g} Msun, t_ff={prob['t_ff']:.4g})\n")
        f.write("\n".join(lines))
    return p


def _run_problem(prob, params_rel, out_abs):
    """Run one problem. Returns (index, returncode, message, wall seconds)."""
    os.makedirs(out_abs, exist_ok=True)
    log = path.join(out_abs, "run.log")
    # mpirun, not srun, even inside a Slurm allocation: concurrent srun steps contend for the
    # step allocation and serialise, which would defeat the whole point of running these in
    # parallel. --bind-to none is essential -- without it every concurrent mpirun binds its ranks
    # to the same low-numbered cores and they fight over them.
    cmd = ["mpirun", "-np", str(RANKS_PER_PROBLEM), "--bind-to", "none",
           "--oversubscribe", "./GIZMO", params_rel, "0"]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    t0 = time.monotonic()
    try:
        with open(log, "w") as fh:
            r = subprocess.run(cmd, cwd=TEST_DIR, stdout=fh, stderr=subprocess.STDOUT,
                               timeout=PROBLEM_TIMEOUT_S, check=False, env=env)
        return (prob["index"], r.returncode, ("ok" if r.returncode == 0 else f"exit {r.returncode}"),
                time.monotonic() - t0)
    except subprocess.TimeoutExpired:
        return prob["index"], -1, f"timeout after {PROBLEM_TIMEOUT_S:.0f}s", time.monotonic() - t0


SYNC_RE = re.compile(r"Energy \(synced,grav\) t=(\S+) E_kin=(\S+) E_pot=(\S+) E_tot=(\S+)")


def _energy_curve(outdir):
    """(times, relative energy error, worst) for one problem, from the in-code synced diagnostic.

    Read from run.log rather than the snapshots. ENERGY_BUDGET_DIAGNOSTIC reports only at full
    synchronization, the one point in the step where every particle's velocity and the potential
    share a clock, so the sum is a genuine conserved quantity. A snapshot never satisfies that:
    io.cc writes Velocities at kick-time against drift-time positions, which is an O(dt/2t_dyn)
    error per particle -- percent-level here. Measured on this suite, the snapshot estimate gave
    3-10% scatter with no secular trend, i.e. pure instrument noise an order of magnitude above
    the 1% being tested; recomputing the potential with an exact brute-force sum removed only half
    of it, because the residual is in the velocities and no choice of potential can reach it.

    Normalised by |E(0)|, not the initial kinetic energy, since these start at rest and KE(0) = 0.
    """
    log = path.join(outdir, "run.log")
    if not path.isfile(log):
        return np.array([]), np.array([]), np.nan
    t, e = [], []
    with open(log, errors="replace") as fh:
        for line in fh:
            mm = SYNC_RE.search(line)
            if mm:
                t.append(float(mm.group(1)))
                e.append(float(mm.group(4)))
    if len(e) < 2:
        return np.array(t), np.array([]), np.nan
    t, e = np.array(t), np.array(e)
    rel = (e - e[0]) / abs(e[0])
    return t, rel, float(np.max(np.abs(rel)))


def _results_path(variant_id):
    return f"{TEST_DIR}/results_{variant_id}.json"


def _plot_summary():
    files = sorted(glob.glob(f"{TEST_DIR}/results_*.json"))
    if not files:
        return
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    marks = {"tree": ("o", "#2a78d6"), "tree_equaldt": ("v", "#1baf7a"),
             "direct_gravity": ("s", "#eb6834"), "direct_equaldt": ("D", "#eda100"),
             "ketju": ("P", "#8e44ad"), "ketju_r200au": ("X", "#c0392b"),
             "ketju_direct": ("*", "#555555")}
    for f in files:
        with open(f) as fh:
            d = json.load(fh)
        vid = d["variant_id"]
        mk, col = marks.get(vid, ("^", "#1baf7a"))
        n = np.array([p["n"] for p in d["problems"]], dtype=float)
        w = np.array([p["worst"] for p in d["problems"]], dtype=float)
        ok = np.isfinite(w) & (w > 0)
        # small horizontal offset per variant so overlapping integer N stay legible
        dx = {"tree": -0.18, "tree_equaldt": -0.06, "direct_gravity": 0.06, "direct_equaldt": 0.18,
              "ketju": -0.3, "ketju_r200au": -0.24, "ketju_direct": 0.3}
        ax.semilogy(n[ok] + dx.get(vid, 0.0), w[ok],
                    mk, ms=5, mew=0, color=col, alpha=0.8, label=vid)
    ax.axhline(ENERGY_TOL, color="k", ls="--", lw=1, label=f"tolerance ({_pct(ENERGY_TOL)})")
    ax.set_xlabel("N bodies in problem")
    ax.set_ylabel(r"worst $|E(t)-E(0)| / |E(0)|$ over the run")
    ax.set_title(f"{TEST_NAME}: per-problem energy conservation\n"
                 "tree = tree force error + integration error; direct = integration error alone",
                 fontsize=9)
    ax.grid(True, axis="y", alpha=0.25, lw=0.6)
    ax.legend(**LEGEND_KW)
    fig.tight_layout()
    fig.savefig(PLOT_PATH, dpi=120)
    plt.close(fig)


def _plot_per_problem():
    """One curve per problem: |dE|/|E0| against dynamical age, log-log, one panel per variant.

    Time is in units of each problem's own free-fall time rather than code units, because the
    problems span 6x in t_ff -- plotted against absolute time the curves would be sheared apart
    by nothing more interesting than how massive each system happens to be. Colour is N, since
    the interesting question is whether error grows with the number of bodies.
    """
    files = sorted(glob.glob(f"{TEST_DIR}/results_*.json"))
    if not files:
        return
    data = []
    for f in files:
        with open(f) as fh:
            data.append(json.load(fh))
    fig, axes = plt.subplots(1, len(data), figsize=(4.6 * len(data), 4.4),
                             sharex=True, sharey=True, squeeze=False)
    ns = [q["n"] for d in data for q in d["problems"]]
    norm = plt.Normalize(min(ns), max(ns))
    cmap = plt.get_cmap("viridis")
    for ax, d in zip(axes[0], data):
        for q in d["problems"]:
            t = np.asarray(q.get("times", []), dtype=float)
            rel = np.abs(np.asarray(q.get("rel", []), dtype=float))
            if t.size < 2 or not q.get("t_ff"):
                continue
            x = t / q["t_ff"]
            good = (x > 0) & (rel > 0)
            if good.any():
                ax.loglog(x[good], rel[good], "-", lw=1.0, alpha=0.75,
                          color=cmap(norm(q["n"])))
        ax.axhline(ENERGY_TOL, color="k", ls="--", lw=1.2)
        # start at a tenth of a free-fall time: before that the systems are still in free fall,
        # nothing has interacted, and the curves only show the integrator ticking over at machine
        # level -- decades of empty axis that squash the part where the error actually develops
        ax.set_xlim(left=0.1)
        ax.set_xlabel(r"$t / t_{\rm ff}$")
        ax.set_title(d["variant_id"], fontsize=10)
        ax.grid(True, which="major", alpha=0.25, lw=0.6)
    axes[0][0].set_ylabel(r"$|E(t)-E(0)| / |E(0)|$")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes[0], pad=0.01); cb.set_label("N bodies")
    axes[0][-1].plot([], [], "k--", lw=1.2, label=f"tolerance ({_pct(ENERGY_TOL)})")
    axes[0][-1].legend(loc="lower right", **LEGEND_KW)
    fig.suptitle(f"{TEST_NAME}: energy error per problem "
                 "(tree = tree + integration error; direct = integration error alone)",
                 fontsize=10)
    fig.savefig(CURVES_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig)


# The four single-factor comparisons the 2x2 design affords. Each holds one thing fixed and
# changes the other, so a panel where the points sit on the diagonal says that factor did not
# matter for these problems, and vertical displacement measures how much it did.
PAIR_PANELS = [
    ("tree",           "direct_gravity", "effect of the tree (individual dt)"),
    ("tree_equaldt",   "direct_equaldt", "effect of the tree (equal dt)"),
    ("tree",           "tree_equaldt",   "effect of unequal dt (tree)"),
    ("direct_gravity", "direct_equaldt", "effect of unequal dt (direct)"),
    ("tree",           "ketju",          "effect of KETJU (tree)"),
    ("direct_gravity", "ketju_direct",   "effect of KETJU (direct)"),
]


def _plot_pairwise():
    """Per-problem worst error of one variant against another, one panel per isolated factor.

    Paired by problem index, so each point is the same initial condition integrated two ways and
    the comparison is not confounded by which problems happen to be hard. Points on the 1:1 line
    mean the factor changed nothing.
    """
    res = {}
    for f in glob.glob(f"{TEST_DIR}/results_*.json"):
        with open(f) as fh:
            d = json.load(fh)
        res[d["variant_id"]] = {q["index"]: q for q in d["problems"]}
    if not res:
        return
    ns = [q["n"] for v in res.values() for q in v.values()]
    norm = plt.Normalize(min(ns), max(ns)); cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(3, 2, figsize=(9.2, 12.9))
    # limits from the data rather than fixed, so the panels are not mostly empty decades; shared
    # across all panels so the 1:1 diagonal means the same thing everywhere and panels stay
    # comparable. The tolerance is forced into range so its line is always visible.
    allw = np.array([q["worst"] for v in res.values() for q in v.values()], dtype=float)
    allw = allw[np.isfinite(allw) & (allw > 0)]
    if allw.size:
        lo = min(allw.min() / 3.0, ENERGY_TOL / 3.0)
        hi = max(allw.max() * 3.0, ENERGY_TOL * 3.0)
    else:
        lo, hi = 1e-6, 1.0
    for ax, (xa, ya, what) in zip(axes.ravel(), PAIR_PANELS):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.plot([lo, hi], [lo, hi], "-", color="0.6", lw=1, zorder=1)
        ax.axhline(ENERGY_TOL, color="k", ls=":", lw=0.8, zorder=1)
        ax.axvline(ENERGY_TOL, color="k", ls=":", lw=0.8, zorder=1)
        ax.set_xlabel(xa.replace("_", r"\_")); ax.set_ylabel(ya.replace("_", r"\_"))
        ax.set_title(what, fontsize=9)
        if xa not in res or ya not in res:
            ax.text(0.5, 0.5, "pending", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="0.45")
            continue
        common = sorted(set(res[xa]) & set(res[ya]))
        x = np.array([res[xa][i]["worst"] for i in common], dtype=float)
        y = np.array([res[ya][i]["worst"] for i in common], dtype=float)
        c = np.array([res[xa][i]["n"] for i in common], dtype=float)
        good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        ax.scatter(x[good], y[good], s=26, c=cmap(norm(c[good])), edgecolors="none",
                   alpha=0.85, zorder=3)
        if good.any():  # >1 means the y-axis variant is the worse one
            r = np.median(y[good] / x[good])
            ax.text(0.04, 0.94, f"median ratio {r:.2f}", transform=ax.transAxes,
                    fontsize=8, va="top")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=axes, pad=0.02, shrink=0.85); cb.set_label("N bodies")
    fig.suptitle(f"{TEST_NAME}: worst $|E(t)-E(0)|/|E(0)|$ per problem, variant vs variant\n"
                 f"grey line is 1:1, dotted lines mark the {_pct(ENERGY_TOL)} tolerance",
                 fontsize=10)
    fig.savefig(PAIRS_PATH, dpi=130, bbox_inches="tight")
    plt.close(fig)


# A 2x2 cross over the two things that can spoil conservation, so each is isolated rather than
# inferred. Reading down a column removes the tree; reading across a row removes the unequal
# timesteps; the bottom-right corner is the integrator on its own with nothing else in the way:
#
#                          individual dt              FORCE_EQUAL_TIMESTEPS
#   tree                   tree + integ + dt          tree + integ
#   SINGLE_STAR_DIRECT_..  integ + dt                 integ alone
#
# FORCE_EQUAL_TIMESTEPS puts every particle on one universal step (the minimum over all of them),
# which is affordable at N<=10 and has a useful side effect: every step is then a full
# synchronization, so the synced energy diagnostic reports on all of them instead of only at the
# rare moments the hierarchy happens to line up.
#
# The KETJU variants (see the module docstring) use individual timesteps only: KETJU_REGULARIZATION
# and FORCE_EQUAL_TIMESTEPS are a compile-time #error together.
@pytest.mark.parametrize("extra_config_flags,ketju_radius", [
    pytest.param((), None, id="tree"),
    pytest.param(("FORCE_EQUAL_TIMESTEPS",), None, id="tree_equaldt"),
    pytest.param(("SINGLE_STAR_DIRECT_GRAVITY",), None, id="direct_gravity"),
    pytest.param(("SINGLE_STAR_DIRECT_GRAVITY", "FORCE_EQUAL_TIMESTEPS"), None, id="direct_equaldt"),
    pytest.param(("KETJU_REGULARIZATION",), 0.005, id="ketju"),
    pytest.param(("KETJU_REGULARIZATION",), 0.001, id="ketju_r200au"),
    pytest.param(("KETJU_REGULARIZATION", "SINGLE_STAR_DIRECT_GRAVITY"), 0.005, id="ketju_direct"),
])
def test_fewbody(extra_config_flags, ketju_radius, request):
    variant_id = request.node.callspec.id.split("-")[0]
    problems = _make_suite()
    base_text = open(f"{TEST_DIR}/{TEST_NAME}.params").read()
    parse_params(f"{TEST_DIR}/{TEST_NAME}.params")  # fail early on a malformed base params file

    extra_params = None
    dir_key = extra_config_flags  # output directory identity
    if ketju_radius is not None:
        extra_params = {**KETJU_PARAMS, "KetjuRegionRadius": ketju_radius}
        # KETJU variants that differ only in R share build flags: key their directories on R too
        dir_key = extra_config_flags + (f"R{ketju_radius:g}",)
    out_root_rel = "output" + variant_suffix(dir_key)   # relative to TEST_DIR
    out_root = variant_output_dir(TEST_NAME, dir_key)   # relative to repo root

    skip_run = bool(os.environ.get("GIZMO_TEST_SKIP_BUILD_RUN"))
    if not skip_run:
        clean_test_outputs(TEST_NAME, dir_key)
        # one build serves every problem in this variant
        build_gizmo_for_test(TEST_NAME, OMP_THREADS, extra_config_flags)

    jobs = []
    for p in problems:
        rel = path.join(out_root_rel, p["name"])
        jobs.append((p, _write_problem_params(base_text, p, rel, extra_params),
                     path.join(out_root, p["name"])))

    failures, walls = [], {}
    if not skip_run:
        n_par = _n_concurrent()
        print(f"\n[{TEST_NAME}/{variant_id}] {len(jobs)} problems, {RANKS_PER_PROBLEM} rank(s) each, "
              f"{n_par} at a time on {_available_cores()} cores")
        with ThreadPoolExecutor(max_workers=n_par) as ex:
            for idx, rc, msg, wall in ex.map(lambda j: _run_problem(j[0], path.relpath(j[1], TEST_DIR), j[2]), jobs):
                walls[idx] = wall
                if rc != 0:
                    failures.append(f"problem {idx}: GIZMO {msg} (see {jobs[idx][2]}/run.log)")

    # collect, whether or not everything ran, so a partial suite still reports what it has
    records, worst_overall, offenders = [], 0.0, []
    for p, _, out_abs in jobs:
        t, rel, worst = _energy_curve(out_abs)
        expect_tmax = N_TFF * p["t_ff"]
        reached = bool(len(t) and abs(t[-1] - expect_tmax) < 1e-6 * expect_tmax)
        records.append({**p, "worst": worst, "n_snaps": int(len(t)),
                        "t_final": (float(t[-1]) if len(t) else None), "reached_timemax": reached,
                        "wall_s": walls.get(p["index"]),
                        "times": [float(x) for x in t], "rel": [float(x) for x in rel]})
        if not np.isfinite(worst):
            failures.append(f"problem {p['index']}: too few snapshots to measure energy "
                            f"({len(t)} in {out_abs})")
            continue
        if not reached:
            failures.append(f"problem {p['index']}: stopped at t={t[-1]:.4g}, expected "
                            f"{expect_tmax:.4g} -- run truncated, energy error not comparable")
        worst_overall = max(worst_overall, worst)
        if worst >= ENERGY_TOL:
            offenders.append(f"  problem {p['index']:3d}: N={p['n']:2d} "
                             f"M={p['m_total']:6.2f} Msun (m {p['m_min']:.2f}-{p['m_max']:.2f}) "
                             f"|dE|/|E0| = {worst:.4f}")

    with open(_results_path(variant_id), "w") as fh:
        json.dump({"variant_id": variant_id, "energy_tol": ENERGY_TOL,
                   "ranks_per_problem": RANKS_PER_PROBLEM, "n_tff": N_TFF,
                   "ketju_params": extra_params,
                   "problems": records}, fh, indent=1)
    _plot_summary()
    _plot_per_problem()
    _plot_pairwise()

    assert not failures, (
        f"{variant_id}: {len(failures)} problem(s) did not produce a usable run:\n  "
        + "\n  ".join(failures[:20]))
    assert not offenders, (
        f"{variant_id}: {len(offenders)}/{len(jobs)} problems exceed {ENERGY_TOL:.0%} energy "
        f"error (worst {worst_overall:.4f}):\n" + "\n".join(offenders[:20])
        + "\nCompare the two variants: if direct_gravity fails too this is the integrator; if "
          "only tree fails it is tree force error.")
