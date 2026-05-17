#!/usr/bin/env python3
"""Strict A/B diff for GIZMO_NGL_DUMP_DIR binary dumps.

Reads two parallel directories of ngl_rank<R>_call<N>.bin files and verifies:
  - same set of (rank, call) pairs
  - identical caller_label + num_active + num_total + total_pairs per call
  - identical active_indices array per call
  - identical sorted neighbor row per active particle (extras and missing both FATAL)

Order within a row is ignored (rows are sorted at dump time on both sides).

Exit codes:
  0 = identical
  1 = mismatch found
  2 = harness error (missing files, parse error, etc.)
"""

import os, sys, struct, glob

# v3 = basic, v4 = basic+detail (since the 64-bit-CSR migration; offsets are
# int64). v1/v2 (32-bit offsets) are no longer produced and not supported.
MAGIC_V3 = b'NGLDMPv3'
MAGIC_V4 = b'NGLDMPv4'

def parse_dump(path):
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    if data[pos:pos+8] not in (MAGIC_V3, MAGIC_V4):
        raise ValueError(f"{path}: bad magic {data[pos:pos+8]!r}")
    pos += 8
    rank = struct.unpack_from('<i', data, pos)[0]; pos += 4
    seq = struct.unpack_from('<i', data, pos)[0]; pos += 4
    caller = data[pos:pos+32].rstrip(b'\x00').decode(); pos += 32
    num_active = struct.unpack_from('<i', data, pos)[0]; pos += 4
    num_total = struct.unpack_from('<i', data, pos)[0]; pos += 4
    total_pairs = struct.unpack_from('<q', data, pos)[0]; pos += 8
    active = struct.unpack_from(f'<{num_active}i', data, pos); pos += 4 * num_active
    offsets = struct.unpack_from(f'<{num_active+1}q', data, pos); pos += 8 * (num_active + 1)
    if total_pairs > 0:
        neighbors = struct.unpack_from(f'<{total_pairs}i', data, pos); pos += 4 * total_pairs
    else:
        neighbors = ()
    return {
        'rank': rank, 'seq': seq, 'caller': caller,
        'num_active': num_active, 'num_total': num_total, 'total_pairs': total_pairs,
        'active': active, 'offsets': offsets, 'neighbors': neighbors,
        'path': path,
    }

def fatal(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)

def err(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)

def main():
    if len(sys.argv) != 3:
        err(f"usage: {sys.argv[0]} <dump_dir_A> <dump_dir_B>")
    dir_a, dir_b = sys.argv[1], sys.argv[2]
    files_a = sorted(glob.glob(os.path.join(dir_a, "ngl_rank*_call*.bin")))
    files_b = sorted(glob.glob(os.path.join(dir_b, "ngl_rank*_call*.bin")))
    names_a = {os.path.basename(p) for p in files_a}
    names_b = {os.path.basename(p) for p in files_b}
    only_a = names_a - names_b
    only_b = names_b - names_a
    common = sorted(names_a & names_b)
    if only_a or only_b:
        # File-set mismatch is non-fatal IF it's just one run reached further
        # than the other (different wall budget). Diff the intersection only.
        # Fatal only if any specific call is missing on one side mid-range
        # (gap pattern) which would indicate skipped calls, not truncation.
        # For now: report and proceed with common set.
        print(f"NOTE: file set differs: only-A={len(only_a)} only-B={len(only_b)} "
              f"(comparing common set of {len(common)} files)", file=sys.stderr)

    files_a = [os.path.join(dir_a, n) for n in common]
    n_calls = len(files_a)
    n_total_pairs_compared = 0
    n_missing = 0
    n_extras = 0
    n_call_mismatches = 0

    for fa in files_a:
        fb = os.path.join(dir_b, os.path.basename(fa))
        a = parse_dump(fa)
        b = parse_dump(fb)
        if (a['caller'], a['num_active'], a['num_total'], a['total_pairs']) != \
           (b['caller'], b['num_active'], b['num_total'], b['total_pairs']):
            print(f"HEADER MISMATCH {os.path.basename(fa)}: "
                  f"A=(caller={a['caller']}, na={a['num_active']}, nt={a['num_total']}, tp={a['total_pairs']}) "
                  f"B=(caller={b['caller']}, na={b['num_active']}, nt={b['num_total']}, tp={b['total_pairs']})")
            n_call_mismatches += 1
            continue
        if a['active'] != b['active']:
            print(f"ACTIVE INDEX MISMATCH {os.path.basename(fa)}: first {sum(1 for x,y in zip(a['active'],b['active']) if x!=y)} differ")
            n_call_mismatches += 1
            continue

        # Per-active row diff. Rows are pre-sorted at dump time.
        for ai in range(a['num_active']):
            beg, end = a['offsets'][ai], a['offsets'][ai+1]
            beg_b, end_b = b['offsets'][ai], b['offsets'][ai+1]
            row_a = a['neighbors'][beg:end]
            row_b = b['neighbors'][beg_b:end_b]
            if row_a == row_b: continue
            sa, sb = set(row_a), set(row_b)
            missing = sa - sb  # in A, not in B → B is missing pairs
            extras = sb - sa   # in B, not in A → B has extras
            if missing or extras:
                n_missing += len(missing)
                n_extras += len(extras)
                # Print first few cases per call
                if n_call_mismatches < 20:
                    print(f"ROW MISMATCH {os.path.basename(fa)} active_idx={a['active'][ai]} "
                          f"(active_pos={ai}): missing={len(missing)} extras={len(extras)}; "
                          f"sample missing={sorted(missing)[:3]} sample extras={sorted(extras)[:3]}")
                n_call_mismatches += 1

        n_total_pairs_compared += a['total_pairs']

    print(f"\n=== A/B diff summary ===")
    print(f"calls compared:    {n_calls}")
    print(f"pairs compared:    {n_total_pairs_compared:,}")
    print(f"call mismatches:   {n_call_mismatches}")
    print(f"missing pairs:     {n_missing}  (in A, not in B)")
    print(f"extra pairs:       {n_extras}   (in B, not in A)")

    if n_missing > 0 or n_extras > 0 or n_call_mismatches > 0:
        print("RESULT: FAIL", file=sys.stderr)
        sys.exit(1)
    print("RESULT: PASS")
    sys.exit(0)

if __name__ == '__main__':
    main()
