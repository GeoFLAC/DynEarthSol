---
name: Bug report
about: A build failure, crash, wrong result or hang
title: ''
labels: bug
---

<!-- Title with the symptom, not a diagnosis: "SIGFPE at startup when
     quality_check_step_interval = 0", not "missing zero guard". Paste output
     as text, not screenshots, so it can be searched and quoted. -->

## Reproduce

<!-- Enough for a maintainer to see the failure without a round-trip. -->

- Version: <!-- `git log -1 --format='%h %s'`, and the branch -->
- Build: <!-- the exact make line, e.g. `make ndims=2 openacc=1`; snapshot.diff has the flags -->
- Run: <!-- the exact command, and attach the cfg -- the smallest one that still fails -->
- Platform: <!-- OS, compiler and version; GPU driver and HPC SDK for openacc=1 -->

## Observed

<!-- Verbatim error text and the exit code. For a wrong result, the numbers as
     printed and the frame they came from. -->

```
```

## Expected

<!-- What should have happened and why you expected it: defaults.cfg, the
     docs, or an earlier version that behaved differently. -->

## Notes

<!-- Optional. What you already ruled out. A diagnosis is welcome, labelled as
     one ("looks like X, unverified"). If it is a regression, name the last
     version that worked. -->
