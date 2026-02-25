# GoSPL Coupling Improvements: Linear Rate Extrapolation & Adaptive Frequency

## Background

The previous coupling scheme applied the full elevation change `Δh` as a lump sum every N DES
steps. This creates impulsive perturbations that grow with N. The new scheme distributes erosion
continuously using linear rate extrapolation and adapts N based on how rapidly erosion rates change.

---

## Diagram 1: The Problem — Old Lump-Sum vs. New Gradual Application

```
OLD BEHAVIOR (lump-sum, N=5 steps):
─────────────────────────────────────────────────────────────────────
 DES steps:   1      2      3      4      5 (GoSPL)  6      7  ...

 Topography:  ─────────────────────────────▼──────────────────────
                                            │
                                         JUMP! Full Δh applied
                                         all at once → impulsive
                                         perturbation, grows with N

NEW BEHAVIOR (gradual + extrapolated, N=5 steps):
─────────────────────────────────────────────────────────────────────
 DES steps:   1      2      3      4      5 (GoSPL)  6      7  ...

 Topography:  ╲      ╲      ╲      ╲      ╲──────── ╲      ╲
               ╲      ╲      ╲      ╲                ╲      ╲
                ─      ─      ─      ─                ─      ─
                ▲ each step applies rate_i × dt      ▲ new rate
                  (no sudden jump, smooth evolution)    computed
```

---

## Diagram 2: Full Coupling Cycle Flowchart

```
Every DES time step:
┌─────────────────────────────────────────────────────┐
│  step_counter++                                      │
│  accumulated_dt += dt                                │
└──────────────────────────┬──────────────────────────┘
                           │
           ┌───────────────▼──────────────────┐
           │  step_counter < adaptive_N ?      │
           └───────────────┬──────────────────┘
                    YES    │    NO
          ┌────────────────┘    └──────────────────────────┐
          │                                                 │
          ▼                                                 ▼
 ┌─────────────────────┐                     ┌─────────────────────────┐
 │  NON-COUPLING STEP  │                     │     COUPLING STEP        │
 │                     │                     │                          │
 │  Apply pending rate │                     │  1. Push DES topo → GoSPL│
 │  with linear extrap │                     │  2. Query elev (before)  │
 │  (see Diagram 3)    │                     │  3. Run GoSPL erosion    │
 │                     │                     │  4. Query elev (after)   │
 │  coord[z] +=        │                     │  5. Δh = after - before  │
 │    rate_i × dt      │                     │  6. rate = Δh / total_dt │
 └─────────────────────┘                     │  7. Apply 1st step: rate×dt│
                                             │  8. Save rate → pending  │
                                             │  9. Adapt N (Diagram 4)  │
                                             │ 10. Reset step_counter=0 │
                                             └─────────────────────────┘
```

---

## Diagram 3: Linear Rate Extrapolation Detail

A **coupling interval** = the N DES steps between two consecutive GoSPL calls.
`prev_rate` = rate computed at the end of the previous interval.
`pending_rate` = rate computed at the end of the most recent interval (used during the next one).

```
Idea: the erosion rate is probably changing smoothly, so project forward.

  Rate
  (m/s)            ← previous         ← current          ← next
   │                 interval           interval            interval
   │                                                       (being filled)
   │   prev_rate ●                    ● pending_rate
   │              ╲                  /  ╲
   │               ╲  (known slope) /    ╲  extrapolated
   │                ╲              /      ╲  trend
   │                 ╲            /        *──────────────*
   │        GoSPL ────●──────────●──────── ▲   steps 1…N  ▲ GoSPL
   │        call                 GoSPL     │               │ call
   │                             call      progress: 0→1  (clipped at 2×)
   │
   └─────────────────────────────────────────────────────── time

   rate_applied(step i) = pending_rate + (pending_rate - prev_rate) × (i / N)
                          ╰─────────────╯ ╰──────────────────────────────────╯
                            base rate           extrapolated adjustment

   Limiter: |rate_applied| ≤ 2 × |pending_rate|   (prevents runaway)
```

---

## Diagram 4: Adaptive Coupling Frequency Logic

```
After each GoSPL call, compute:

        √( Σ (rate_new - rate_old)² )
metric = ─────────────────────────────
              √( Σ rate_new² )

  = relative L2 change in erosion rates between cycles


           metric > tolerance (0.3)?
          ┌─────────────────────────┐
    YES   │                         │   NO
  ┌───────┘                         └────────────────────┐
  │                                                       │
  ▼                                                       ▼
Rates changing FAST                           metric < tolerance / 4 ?
→ Need more coupling                          ┌──────────────────────┐
                                        YES   │                      │ NO
  N = max(1, N / 2)              ┌────────────┘                      │
  (couple more often)            │                                    │
                                 ▼                                    ▼
                       Rates very stable                    Keep N unchanged
                       → Save compute
                       N = min(4×N_base, N × 2)
                       (couple less often)


  Bounds: 1 ≤ adaptive_N ≤ 4 × base_N
  Set tolerance = 0 → disables adaptive behavior, fixed N
```

---

## Diagram 5: Memory Layout Per Node (State Machine)

```
After each GoSPL call, the driver holds:

  ┌──────────────────────────────────────────────────────┐
  │  GoSPLDriver state (per surface node i)               │
  │                                                       │
  │  prev_erosion_rate[i]    ←── rate from the           │
  │                               previous interval       │
  │                                                       │
  │  pending_erosion_rate[i] ←── rate from the most      │
  │                               recent interval (m/s)   │
  │                                                       │
  │  has_prev_rate = true/false                           │
  │    (false during 1st interval: no extrapolation yet)  │
  └──────────────────────────────────────────────────────┘

  Each non-coupling step reads both arrays to
  compute the extrapolated rate for that node.

  At the start of the next GoSPL call:
    prev    ← pending        (shift history by one interval)
    pending ← new rates from GoSPL
```

---

## Configuration Parameters

| Parameter                    | Default | Description                                              |
|------------------------------|---------|----------------------------------------------------------|
| `gospl_coupling_frequency`   | `1`     | Base N: run GoSPL every N DES steps                      |
| `gospl_rate_change_tolerance`| `0.3`   | L2 rate-change threshold for adaptive N (0 = disabled)   |

### Behavioral summary by N

| `gospl_coupling_frequency` | `gospl_rate_change_tolerance` | Behavior                                        |
|----------------------------|-------------------------------|-------------------------------------------------|
| `1`                        | any                           | Identical to original (no extrapolation needed) |
| `N > 1`                    | `0`                           | Fixed N, gradual application with extrapolation |
| `N > 1`                    | `0.3` (default)               | Adaptive N + gradual application                |

---

## See Also

- [bc.cxx](../bc.cxx) — `use_gospl()` implementation (lines ~1589–1750)
- [gospl-driver.hpp](gospl-driver.hpp) — Driver state members
- [GOSPL_COUPLING.md](../GOSPL_COUPLING.md) — Full coupling documentation
- [implementation_plan.md](implementation_plan.md) — Original design document
