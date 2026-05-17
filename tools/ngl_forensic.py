#!/usr/bin/env python3
"""Forensic analysis of NGL dump mismatches.

For the first call where dumps A and B disagree, compute for each missing/extra pair:
  r² (using each side's positions), cutoff² (per search_mode), margin = (r² - cutoff²)/cutoff².

Phil's rule: never call something "FP error" without showing it at the actual FP
floor (~1e-15 relative for double, ~1e-7 for float). compact_xyzh is float32 so
the relevant floor on h-related quantities is ~1e-7. Any margin much larger than
1e-7 is NOT FP — real divergence.

Usage:
  ngl_forensic.py <dump_dir_A> <dump_dir_B> [call_seq_to_diff]
"""

import os, sys, struct, glob, math

# v3 = basic, v4 = basic+detail (since the 64-bit-CSR migration; offsets are
# int64). v1/v2 (32-bit offsets) are no longer produced and not supported.
MAGIC_V3 = b'NGLDMPv3'
MAGIC_V4 = b'NGLDMPv4'

# Search modes (must match neighbor_list.h)
NGB_SEARCH_ONEWAY = 0
NGB_SEARCH_SYMMETRIC = 1

def parse_dump(path):
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    magic = data[pos:pos+8]; pos += 8
    if magic == MAGIC_V3: detail = False
    elif magic == MAGIC_V4: detail = True
    else: raise ValueError(f"{path}: bad magic {magic}")
    rank = struct.unpack_from('<i', data, pos)[0]; pos += 4
    seq = struct.unpack_from('<i', data, pos)[0]; pos += 4
    caller = data[pos:pos+32].rstrip(b'\x00').decode(); pos += 32
    num_active = struct.unpack_from('<i', data, pos)[0]; pos += 4
    num_total = struct.unpack_from('<i', data, pos)[0]; pos += 4
    total_pairs = struct.unpack_from('<q', data, pos)[0]; pos += 8
    active = list(struct.unpack_from(f'<{num_active}i', data, pos)); pos += 4 * num_active
    offsets = list(struct.unpack_from(f'<{num_active+1}q', data, pos)); pos += 8 * (num_active + 1)
    if total_pairs > 0:
        neighbors = list(struct.unpack_from(f'<{total_pairs}i', data, pos))
        pos += 4 * total_pairs
    else:
        neighbors = []

    detail_dict = None
    if detail:
        search_mode = struct.unpack_from('<i', data, pos)[0]; pos += 4
        srf = struct.unpack_from('<d', data, pos)[0]; pos += 8
        h_inflate = struct.unpack_from('<d', data, pos)[0]; pos += 8
        pflags = list(struct.unpack_from('<3i', data, pos)); pos += 12
        box_sizes = list(struct.unpack_from('<3d', data, pos)); pos += 24
        # Per-active: 3 doubles pos + 1 double h_q + 1 double h_kernel
        active_pos, active_hq, active_hk = [], [], []
        for a in range(num_active):
            xyz = struct.unpack_from('<3d', data, pos); pos += 24
            hq = struct.unpack_from('<d', data, pos)[0]; pos += 8
            hk = struct.unpack_from('<d', data, pos)[0]; pos += 8
            active_pos.append(xyz); active_hq.append(hq); active_hk.append(hk)
        # Per-neighbor: 3 doubles pos + 1 double hk + 1 float compact_h
        ngb_pos, ngb_hk, ngb_ch = [], [], []
        for p in range(total_pairs):
            xyz = struct.unpack_from('<3d', data, pos); pos += 24
            hk = struct.unpack_from('<d', data, pos)[0]; pos += 8
            ch = struct.unpack_from('<f', data, pos)[0]; pos += 4
            ngb_pos.append(xyz); ngb_hk.append(hk); ngb_ch.append(ch)
        detail_dict = dict(search_mode=search_mode, srf=srf, h_inflate=h_inflate,
                           pflags=pflags, box_sizes=box_sizes,
                           active_pos=active_pos, active_hq=active_hq, active_hk=active_hk,
                           ngb_pos=ngb_pos, ngb_hk=ngb_hk, ngb_ch=ngb_ch)
    return dict(rank=rank, seq=seq, caller=caller, num_active=num_active, num_total=num_total,
                total_pairs=total_pairs, active=active, offsets=offsets, neighbors=neighbors,
                detail=detail_dict, path=path)

def periodic_dr(dx, box, periodic):
    if not periodic or box <= 0: return dx
    if dx > 0.5 * box: dx -= box
    elif dx < -0.5 * box: dx += box
    return dx

def compute_r2(pos_a, pos_b, pflags, box):
    dx = periodic_dr(pos_a[0] - pos_b[0], box[0], pflags[0])
    dy = periodic_dr(pos_a[1] - pos_b[1], box[1], pflags[1])
    dz = periodic_dr(pos_a[2] - pos_b[2], box[2], pflags[2])
    return dx*dx + dy*dy + dz*dz

def parse_header_only(path):
    """Parse just the call metadata + neighbor list (no detail) — cheap."""
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    magic = data[pos:pos+8]; pos += 8
    if magic not in (MAGIC_V3, MAGIC_V4):
        raise ValueError(f"{path}: bad magic")
    rank = struct.unpack_from('<i', data, pos)[0]; pos += 4
    seq = struct.unpack_from('<i', data, pos)[0]; pos += 4
    caller = data[pos:pos+32].rstrip(b'\x00').decode(); pos += 32
    num_active = struct.unpack_from('<i', data, pos)[0]; pos += 4
    num_total = struct.unpack_from('<i', data, pos)[0]; pos += 4
    total_pairs = struct.unpack_from('<q', data, pos)[0]; pos += 8
    return num_active, num_total, total_pairs, caller

def parse_neighbor_section(path):
    """Parse only header + active + offsets + neighbors. Skip the detail tail.
    Cheap — no full unpacking of position arrays."""
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    magic = data[pos:pos+8]; pos += 8
    if magic not in (MAGIC_V3, MAGIC_V4):
        raise ValueError(f"{path}: bad magic")
    rank = struct.unpack_from('<i', data, pos)[0]; pos += 4
    seq = struct.unpack_from('<i', data, pos)[0]; pos += 4
    caller = data[pos:pos+32].rstrip(b'\x00').decode(); pos += 32
    num_active = struct.unpack_from('<i', data, pos)[0]; pos += 4
    num_total = struct.unpack_from('<i', data, pos)[0]; pos += 4
    total_pairs = struct.unpack_from('<q', data, pos)[0]; pos += 8
    # Slice rather than unpack — keep raw bytes; tuple-compare on raw is fast.
    active_bytes = data[pos:pos + 4 * num_active]; pos += 4 * num_active
    offsets_bytes = data[pos:pos + 8 * (num_active + 1)]; pos += 8 * (num_active + 1)
    neighbors_bytes = data[pos:pos + 4 * total_pairs] if total_pairs > 0 else b''
    return (caller, num_active, num_total, total_pairs,
            active_bytes, offsets_bytes, neighbors_bytes)

def find_first_neighbor_set_diff(dir_a, dir_b):
    """Find the first call where the neighbor SET differs (rows have different
    members). Detail-only differences (positions, compact_h) don't trigger.
    Fast: byte-compares the neighbor-list sections, only full-parses on the
    first hit."""
    files_a = sorted(glob.glob(os.path.join(dir_a, "ngl_rank*_call*.bin")))
    common = sorted(set(os.path.basename(p) for p in files_a)
                  & set(os.path.basename(p) for p in glob.glob(os.path.join(dir_b, "ngl_rank*_call*.bin"))))
    for name in common:
        pa = os.path.join(dir_a, name)
        pb = os.path.join(dir_b, name)
        sa = parse_neighbor_section(pa)
        sb = parse_neighbor_section(pb)
        if sa != sb:
            print(f"  first neighbor-set diff at {name}: "
                  f"caller={sa[0]} A_tp={sa[3]} B_tp={sb[3]}", file=sys.stderr)
            return name, parse_dump(pa), parse_dump(pb)
    return None, None, None

def find_first_diff_call(dir_a, dir_b):
    return find_first_neighbor_set_diff(dir_a, dir_b)

def analyze_pair(name_a_or_b, side_label, a_active_pos, a_active_hq, a_active_hk,
                 j_pos, j_hk, search_mode, srf, h_inflate, pflags, box):
    """Compute r², cutoff² for active vs j on side `side_label`."""
    r2 = compute_r2(a_active_pos, j_pos, pflags, box)
    h_q_eff = a_active_hq * srf  # search_radius_factor inflation
    if search_mode == NGB_SEARCH_ONEWAY:
        cutoff2 = h_q_eff * h_q_eff
    else:  # SYMMETRIC
        # Walker uses max(h_q, h_j) where h_j comes from compact_xyzh
        # which is (P[j].KernelRadius * h_inflate). We compute both forms.
        h_j_compact = j_hk * h_inflate  # what the walker actually uses for h_j
        h_eff = max(h_q_eff, h_j_compact)
        cutoff2 = h_eff * h_eff
    return r2, cutoff2

def main():
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} <dir_A> <dir_B>", file=sys.stderr)
        sys.exit(2)
    dir_a, dir_b = sys.argv[1], sys.argv[2]

    name, a, b = find_first_diff_call(dir_a, dir_b)
    if name is None:
        print("No mismatching call found in common set.")
        sys.exit(0)

    print(f"=== first divergent call: {name} ===")
    print(f"A: caller={a['caller']} num_active={a['num_active']} num_total={a['num_total']} tp={a['total_pairs']}")
    print(f"B: caller={b['caller']} num_active={b['num_active']} num_total={b['num_total']} tp={b['total_pairs']}")

    if a['detail'] is None or b['detail'] is None:
        print("ERROR: detail dump missing — re-run with GIZMO_NGL_DUMP_DETAIL=1")
        sys.exit(2)
    if a['active'] != b['active']:
        print("ERROR: active index lists differ — cannot align by active position")
        sys.exit(2)
    if a['num_active'] != b['num_active']:
        print(f"ERROR: num_active differs: A={a['num_active']} B={b['num_active']}")
        sys.exit(2)

    da, db = a['detail'], b['detail']
    sm = da['search_mode']
    srf = da['srf']
    hinf = da['h_inflate']
    pflags = da['pflags']
    box = da['box_sizes']
    sm_label = {0: "ONEWAY", 1: "SYMMETRIC"}.get(sm, f"sm={sm}")
    print(f"search_mode={sm_label} srf={srf} h_inflate={hinf}")
    print(f"box={box} periodic_flags={pflags}")

    # Per-active row diff. For each mismatched pair, compute the margin using
    # BOTH A's and B's recorded positions/h. The walker's pruning happens on
    # whichever side's compact_xyzh values; the question is whether the side
    # that EXCLUDED the pair has r² > cutoff² IN ITS OWN VIEW.
    rows_with_diff = 0
    margin_log = []
    for ai in range(a['num_active']):
        beg_a, end_a = a['offsets'][ai], a['offsets'][ai+1]
        beg_b, end_b = b['offsets'][ai], b['offsets'][ai+1]
        row_a = set(a['neighbors'][beg_a:end_a])
        row_b = set(b['neighbors'][beg_b:end_b])
        if row_a == row_b: continue
        rows_with_diff += 1
        active_idx = a['active'][ai]
        a_pos_q = da['active_pos'][ai]; a_hq = da['active_hq'][ai]
        b_pos_q = db['active_pos'][ai]; b_hq = db['active_hq'][ai]

        # Lookups by particle j → its recorded (pos, hk, compact_h) on each side.
        idx_to_a, idx_to_b = {}, {}
        for k in range(beg_a, end_a):
            idx_to_a[a['neighbors'][k]] = (da['ngb_pos'][k], da['ngb_hk'][k], da['ngb_ch'][k])
        for k in range(beg_b, end_b):
            idx_to_b[b['neighbors'][k]] = (db['ngb_pos'][k], db['ngb_hk'][k], db['ngb_ch'][k])

        missing = row_a - row_b   # A includes, B does not
        extras  = row_b - row_a   # B includes, A does not

        # For "missing": A included it (pass), B excluded it (would-be-fail).
        # We need B's view of this pair (B's positions/h for active and j).
        # B's row doesn't contain j (it was excluded), so we don't have j's
        # pos/h from B's dump. Use the pos/h relationship: positions across
        # the two runs differ by at most the drift accumulated by call 19.
        # The most informative data we have: A's view (pair was kept) gives
        # the lower-bound margin. If A's margin is far from zero (1%), it's
        # not a near-boundary FP case in A — definitive.
        for j in sorted(missing):
            jp_a, jhk_a, jch_a = idx_to_a[j]
            r2_a, cut2_a = analyze_pair(a['path'], 'A', a_pos_q, a_hq, da['active_hk'][ai],
                                         jp_a, jhk_a, sm, srf, hinf, pflags, box)
            margin_a = (r2_a - cut2_a) / cut2_a if cut2_a > 0 else float('inf')
            # Compute B's view using B's active position; B's j position is
            # unknown (j absent from B's row), so use A's j as best estimate.
            # Position drift between runs is the question; if A's r² is far
            # below cutoff, the pair is robustly an A-include and the issue
            # is B's exclusion, not a near-boundary flip.
            r2_b_est, cut2_b = analyze_pair(b['path'], 'B(est)', b_pos_q, b_hq, db['active_hk'][ai],
                                             jp_a, jhk_a, sm, srf, hinf, pflags, box)
            margin_b_est = (r2_b_est - cut2_b) / cut2_b if cut2_b > 0 else float('inf')
            # Position drift between A and B for the active query
            d2 = sum((a_pos_q[k] - b_pos_q[k])**2 for k in range(3))
            pos_drift = math.sqrt(d2)
            margin_log.append((min(abs(margin_a), abs(margin_b_est)),
                              'A_includes_B_excludes', active_idx, j,
                              math.sqrt(r2_a), math.sqrt(cut2_a), margin_a,
                              margin_b_est, pos_drift))

        for j in sorted(extras):
            jp_b, jhk_b, jch_b = idx_to_b[j]
            r2_b, cut2_b = analyze_pair(b['path'], 'B', b_pos_q, b_hq, db['active_hk'][ai],
                                         jp_b, jhk_b, sm, srf, hinf, pflags, box)
            margin_b = (r2_b - cut2_b) / cut2_b if cut2_b > 0 else float('inf')
            r2_a_est, cut2_a = analyze_pair(a['path'], 'A(est)', a_pos_q, a_hq, da['active_hk'][ai],
                                             jp_b, jhk_b, sm, srf, hinf, pflags, box)
            margin_a_est = (r2_a_est - cut2_a) / cut2_a if cut2_a > 0 else float('inf')
            d2 = sum((a_pos_q[k] - b_pos_q[k])**2 for k in range(3))
            pos_drift = math.sqrt(d2)
            margin_log.append((min(abs(margin_b), abs(margin_a_est)),
                              'B_includes_A_excludes', active_idx, j,
                              math.sqrt(r2_b), math.sqrt(cut2_b), margin_b,
                              margin_a_est, pos_drift))

    print(f"\nrows with diff: {rows_with_diff} / {a['num_active']}")
    print(f"total mismatched pair entries: {len(margin_log)}")

    # Sort by min |margin| (i.e. the side that's CLOSEST to flipping)
    margin_log.sort()
    print(f"\n=== mismatched pair details ===")
    print(f"{'kind':<24} {'active':<8} {'j':<8} {'r_inc':<12} {'cut_inc':<12} {'m_inc':<11} {'m_exc(est)':<12} {'pos_drift':<10}")
    for m_min, kind, ai, j, r_inc, cut_inc, m_inc, m_exc_est, pdrift in margin_log[:50]:
        print(f"{kind:<24} {ai:<8} {j:<8} {r_inc:<12.4e} {cut_inc:<12.4e} {m_inc:<+11.3e} {m_exc_est:<+12.3e} {pdrift:<10.3e}")

    # Classification per Phil's rule
    # compact_xyzh is float32 → relevant floor on h-related quantities is ~1e-7
    # For double-precision positions, position-comparison floor is ~1e-15
    # The cutoff involves h (float32 precision), so use 1e-7 as threshold
    FP_FLOOR_FLOAT = 1e-7
    # Use min(|margin_inc|, |margin_exc_est|) — the side closer to the
    # boundary is the relevant FP-floor question.
    near_boundary = [tup for tup in margin_log if tup[0] < FP_FLOOR_FLOAT]
    real_drift   = [tup for tup in margin_log if tup[0] >= FP_FLOOR_FLOAT]

    print(f"\n=== classification (FP-float floor = {FP_FLOOR_FLOAT}) ===")
    print(f"near-boundary (|margin| < 1e-7): {len(near_boundary)}")
    print(f"real-drift    (|margin| >= 1e-7): {len(real_drift)}")
    if real_drift:
        print(f"\n*** {len(real_drift)} mismatched pair(s) have margin >> FP floor — NOT FP nondeterminism ***")
        sys.exit(1)
    else:
        print(f"\n=== all mismatches consistent with float-precision near-boundary ===")
        sys.exit(0)

if __name__ == '__main__':
    main()
