# Fluid Point-Source Validation

This focused functional suite checks the public fluid injection path:

- 2D `points_y` remains a legacy alias for `points_z`, while specifying both
  or neither is rejected;
- non-finite coordinates, rates, amounts, and invalid or overflowing schedules
  fail as configuration errors;
- box coordinates use the same relative boundary tolerance as the runtime
  simplex lookup;
- every enabled point is located in the actual current mesh, even before its
  schedule starts;
- shared-boundary ownership is deterministic, and single- and multi-point
  source totals are conserved by nodal assembly.
- nonzero sources on fixed-pressure nodes or non-finite/non-positive hydraulic
  mass are rejected instead of being silently discarded;
- PT-driven moving-mesh remeshing fails closed until its pre-existing hydraulic
  remesh lifecycle is safe.

Run the 2D suite from the repository root:

```bash
python tests/functional/point_source_validation/run_point_source_validation.py --clean
```

Pass a separately built 3D executable to add the 3D coordinate and shared-face
cases:

```bash
python tests/functional/point_source_validation/run_point_source_validation.py \
  --clean --exe-3d /path/to/dynearthsol3d
```
