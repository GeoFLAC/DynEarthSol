# MMG remeshing design notes

The freeze policy in [remeshing.cxx](../remeshing.cxx) limits interpolation of
fields such as localized plastic strain by preserving quiet regions.

## Rejected alternatives

The freeze-policy development recorded the following outcomes. [Original design record](https://github.com/GeoFLAC/DynEarthSol/blob/8ded239cd1541d37c9cb43c849db753411171635/remeshing.cxx#L1359-L1367)

| Alternative | Recorded outcome |
| --- | --- |
| An unbounded flood fill to free undersized elements, or using `init_elem_size_n` directly in R5 | Freed MMG's own gradation transition again at every remesh. |
| Expanding the region freed by R1 around damage from boundary flattening on the direct remesh path | Produced NaNs. |
| Dissolving small required islands | Increased bad outputs and retries. |
| Using the metric floor as R1's tiny-element cutoff | Freed nodes pinned by R2 and diffused the localized strain band. |
| Adding an extra ring of hinge elements | Had no effect. |
| Using `SPLIT` for the R1 quality check | Distortion that needed node motion reached the quality gate's broader unfreeze. This increased remeshing without improving field preservation. |
