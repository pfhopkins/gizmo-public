TODO
====

Features
--------

Optimization
------------
- The divide and conquer code has a lot of dynamic allocation, which might be
  avoidable. May not be significant for performance, though.

- `compute_T` might benefit from parallelization, but not clear whether comms
  overhead too big for this to be useful. If implemented, needs careful benchmarking. 

- `get_gbs_rows_for_groups` could account for uneven group sizes for more even
  computational load. But generally probably a fairly small effect.

? Currently the timers do not seem to capture enough things, for large task numbers some
  30% of the time is missed from the subcategories

Other
-----

+ Automated tests
    + Basics setup and some tests done
    ? More comprehensive testing

- Automated benchmarks
