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
* The same sources build `dynearthsol2d` and `dynearthsol3d`, so every change
  must compile with both `make ndims=2` and `make ndims=3`; use `NDIMS`,
  `NODES_PER_ELEM` and `NSTR`, never a literal 2, 3, 4 or 6.
* Parallel loops use `#pragma omp parallel for default(none)` and declare every
  variable in `shared`/`private`/`firstprivate`, so nothing is shared by
  accident. What a clause list may contain changed between g++ 8 and 9, so CI
  compiles with g++ 8 to 15 to keep the pragmas portable. Keep the paired
  `#pragma acc parallel loop` clauses in step.
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

## Commit messages

The message records what the code cannot: what was wrong before, why this fix
and not another, and what it was measured to do. Conventions follow the Git and
Linux kernel patch guides and Chris Beams' seven rules; the split with code
comments follows Hutterer's "On commit messages".

```
type(scope): imperative summary, 72 characters or fewer, no period

The problem, in the present tense, and its user-visible effect. Whether this is
a one-line fix or a new feature, there is a problem that motivated it -- state
it first, so a reader knows whether to read on.

The approach, only where the diff does not make it obvious, and any alternative
considered and rejected. Numbers for every claim: a mesh that changed, a field
that moved, a timing -- give the case and the value.

Co-Authored-By: Name <email>
```

* `type` is one of `fix`, `feat`, `refactor`, `test`, `docs`, `ci`, `build`,
  `cleanup`; `scope` is optional and names the subsystem (`rsf`, `mmg`, `tpf`).
  Lowercase after the colon.
* **Comment versus commit.** A code comment states the invariant the code must
  hold *now* and what breaks if it is violated; it is bound to the code and rots
  when it carries history. The commit body holds the history: the wrong
  behaviour it replaced, the reasoning, the measurements. Write each fact in one
  place. When a fact must live in the code -- an ordering contract, a guard's
  reason -- put it in the comment or an assertion, and let the body say so
  ("asserted at the consumer") rather than repeat it.
* Refer to earlier work by its subject line or pull request number, never by
  hash alone: branches here are rebased and replayed, so hashes go stale.
* Trailers: `Co-Authored-By` for every AI tool or person who wrote part of the
  change; `Fixes #N` when the commit closes an issue.
* Length is set by the problem, not by habit. Most fixes fit three short
  paragraphs; if the body restates the diff's comments, cut one of the two.

## Issues and pull requests

Templates under [`.github/`](.github/) fill in when you open an issue or a pull
request; each says in its comments what it asks for and why, so that is not
repeated here.

* **Bugs**: one issue per symptom, with the version, the exact build and run
  lines, the smallest failing cfg, and the output as text. A diagnosis is
  welcome, labelled as one.
* **Enhancements**: the problem as it stands today, a concrete proposal, the
  alternatives you rejected, and the state that closes the issue.
* **Pull requests**: one topic per PR, commits that each build and pass
  `benchmarks-cores`, and a title that stands alone in `git log` -- the area
  and what changed, never "various fixes". The template's one DES-specific
  question is whether the change altered results -- CI cannot tell -- so run
  the comparison above and report it. Open as a draft until CI is green, then
  request review from at least one other developer.

## AI-assisted contributions

AI-assisted contributions are welcome. You are responsible for the output of
the AI tools you use, as you are for anything else you submit: review it before
you ask anyone else to, and be able to explain every line. Have your assistant
read this file and the templates under [`.github/`](.github/) before it drafts
an issue or a pull request; [`AGENTS.md`](AGENTS.md) at the repository root
tells it to. Disclose the tools used and what they produced in the pull
request's AI disclosure section, and name them in `Co-Authored-By` trailers as
the commit guideline asks. Review is a conversation between people: you may use
AI to improve your wording, but the reply to a review comment is yours, and a
pull request opened by an agent without a person behind it will not be
reviewed.

## Development and release workflow

1. **Prepare feature branch**: Develop your feature and ensure all local tests pass.
2. **Update Changelog**: Document the new features, bug fixes, and breaking changes in `CHANGELOG.md`.
3. **Acquire a draft DOI**: Go to [Zenodo](https://zenodo.org) and acquire a draft DOI for the new version (currently maintained by @chaseshyu).
4. **Update Versioning**:
   - Update `version`, `doi`, and `date-released` in `CITATION.cff`.
   - Update `description` in `.zenodo.json` for the release description on Zenodo.
5. **Create and Review PR**: Create a Pull Request against the master branch. Ensure CI/CD tests pass and request a code review from at least one other developer.
6. **Merge PR**: Once approved, merge the feature branch into master.
7. **Push the tag**: this is the only command-line step. Tag the merge commit and push it; a GitHub Actions workflow, triggered by `v*` and `benchmarks-v*` tags, then packages the source code—specifically including all submodules and a shallow `.git`—and attaches the tarball to a **draft** release. Both halves of that workflow, the packaging here and the Zenodo upload in step 9, run the file as it exists in the tagged commit, so tag a commit that already carries it.

   ```bash
   git tag v2.0.2 && git push origin v2.0.2

   # to start over: delete the draft release on the website, then the tag here
   git push origin --delete v2.0.2
   ```

8. **Compose the release on GitHub**: once the workflow finishes, the draft release it created is waiting under [Releases](https://github.com/GeoFLAC/DynEarthSol/releases), tarball already attached. Edit it there—set the title, summarize the changelog in the notes—and publish. Skip *Generate release notes* if the notes already list the merged pull requests, or the list appears twice. Publishing is the point of no return: a published release is immutable, so confirm the tarball is on the draft first, because no asset can be added afterwards.
9. **Publish the Zenodo deposit**: the tarball lands in the deposit named by `CITATION.cff`, linked from the workflow log. Review it, add *Repository URL* under Software by hand—the deposit API cannot—then publish.

A `benchmarks-v*` tag releases the benchmark dataset from its own branch, which
carries its own `CITATION.cff` DOI and `.zenodo.json`.

A tag pushed before this workflow existed can still be packaged: run *Release
with submodules and Zenodo* from the Actions tab with the tag as its input.
Publishing such a release runs the old file from its own tag, though, so
re-cutting the tag is usually the better move.

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
