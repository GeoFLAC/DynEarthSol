# Development Notes

The main repository is on:
  https://github.com/GeoFLAC/DynEarthSol

To get the most timely progress, using Github:
  git clone https://github.com/GeoFLAC/DynEarthSol.git

## To-do list:

High priority:
* Simple benchmarks for rheology
* Local remeshing
* Different init_marker_spacing in regions

Low priority:
* Save output as vtk format directly
* Stress BC, esp no-stress sidewall
* Heatflux BC
* Frictional heating
* Adiabatic cooling (adiabatic temperature profile)
* Internal heating
* Volume changed induced stress

## Design notes:

* Avoid C++ stream for bulk output, as stream is slower than C-style IO.
* Avoid creating/destroying objects in inner for-loops.
* Avoid static variables and global variables.
* Meaning of error codes. First digit is the category, second the specific
  cause; 1x is the user's to fix, 2x the environment, 3x-6x ours. Exit via
  `die(EXIT_...)` in `utils.hpp`, which prints the code and its category, so
  never `exit()` a bare number.
   0: Normal exit
  10: Config error          11: bad value/unknown option  12: malformed .poly/.exo
  20: Cannot open file      21: read/write failed (HDF5)  22: restart mismatch
  30: Unsupported in this NDIMS                           31: library not built in
  40: Triangle/TetGen       41: MMG                       42: mesh quality/topology
  50: NaN/non-finite        51: marker/geometry lookup    52: resource exhausted
  60: Assertion violated    61: unreachable branch
  These are renumbered: earlier builds used 1 input, 2 IO, 10 triangulation,
  11 runtime, 12 assertion, so 10-12 mean something different in older logs.

## Benchmarks and regression testing

Regression tests live in `benchmarks-cores/`. All targets run from inside that
directory; see the comment block at the top of `benchmarks-cores/Makefile` for
the full variable and target reference.

```bash
make set                # run once and store as reference (creates orig-<CASE>/)
make cmp                # re-run and diff against reference  (most common)
make fresh-restart-cmp  # fresh run → checkpoint → restart → diff (determinism check)
```

`compare.py` is called automatically by `cmp`, `restart`, and
`fresh-restart-cmp`. It prints the relative difference (max, stddev) for every
field and exits 1 if any field exceeds 1e-8 or contains NaN/Inf. Run
`python compare.py --help` (or read its module docstring) for manual usage and
restart-troubleshooting instructions.

## Development and release workflow

1. **Prepare feature branch**: Develop your feature and ensure all local tests pass.
2. **Update Changelog**: Document the new features, bug fixes, and breaking changes in `CHANGELOG.md`.
3. **Acquire a draft DOI**: Go to [Zenodo](https://zenodo.org) and acquire a draft DOI for the new version (currently maintained by @chaseshyu).
4. **Update Versioning**:
   - Update `version`, `doi`, and `date-released` in `CITATION.cff`.
   - Update `description` in `.zenodo.json` for the release description on Zenodo.
5. **Create and Review PR**: Create a Pull Request against the master branch. Ensure CI/CD tests pass and request a code review from at least one other developer.
6. **Merge PR**: Once approved, merge the feature branch into master.
7. **Push the tag**: Tag the merge commit with the version number and push the tag (`git tag v2.0.2 && git push origin v2.0.2`). A GitHub Actions workflow packages the source code—specifically including all submodules—and attaches the tarball to a **draft** release. Both halves of that workflow—the packaging here and the Zenodo upload in step 8—run the file as it exists in the tagged commit, so tag a commit that already carries it. Dispatching the workflow by hand packages a tag pushed earlier, but publishing that release still runs the old file from its tag, so a tag predating this workflow is better deleted and re-cut.
8. **Publish the Release**: Open the draft on GitHub, set the title, summarize the changelog in the release notes, and publish. Published releases are immutable, so the tarball has to be attached before this point—assets cannot be added afterwards. Publishing uploads that same tarball to the Zenodo draft whose DOI is in `CITATION.cff`, for review and manual publication there.

Steps 7 and 8 from the command line, for tag `v2.0.2`:

```bash
git checkout master && git pull            # 7. tag the merge commit, which
git tag v2.0.2 && git push origin v2.0.2   #    packages the draft release
gh run list --workflow=release-with-submodules-and-zenodo.yml --limit 1
gh release view v2.0.2 --json isDraft,assets   # draft, with the tarball on it?

# 8. publish, which sends that same tarball to the Zenodo draft
gh release edit v2.0.2 --draft=false \
  --title "DynEarthSol version 2.0.2" --notes-file release-notes-v2.0.2.md

# package a tag pushed before the workflow existed, or re-package one
gh workflow run release-with-submodules-and-zenodo.yml --ref master -f tag=v2.0.2

# start over -- possible while it is a draft, impossible once published
gh release delete v2.0.2 --cleanup-tag --yes
```

A `benchmarks-v*` tag, for a benchmark dataset release, is packaged into a draft
release the same way; only step 8's Zenodo upload is skipped, because those tags
reserve no DOI in `CITATION.cff`.

A tag name that has ever carried a published release is spent: deleting the
release frees the tag but not the name, and neither does deleting the repository
([immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases)).
A release that has to be redone therefore takes the next version number, which
is why 2.0.1 does not exist.

## Git Submodule Workflow

Since this repository relies on Git submodules (e.g., `knn-bvh`), your development workflow needs to account for them. Submodules are essentially pointers to specific commits in other repositories.

**1. Pulling the latest changes:**
When you pull updates from the main repository, Git does not automatically update the submodule contents by default. To fetch and update everything in one command, use:
```bash
git pull --recurse-submodules
```
*(If you simply ran `git pull` and notice your submodules are out of sync, you can fix it by running `git submodule update`.)*

**2. Modifying a submodule:**
By default, submodules are in a "detached HEAD" state. If you need to modify the code inside a submodule:
1. `cd` into the submodule directory.
2. Checkout the appropriate branch (e.g., `git checkout main`).
3. Make your changes, `git add`, `git commit`, and `git push` from **within** the submodule directory.
4. `cd` back to the main repository root. The main repo will now see that the submodule pointer has changed. 
5. `git add` the submodule directory and commit this pointer update to the main repository.
