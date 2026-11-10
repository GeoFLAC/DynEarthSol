# Contributing to DynEarthSol

The main repository is https://github.com/GeoFLAC/DynEarthSol. Clone it with
its submodules:

```bash
git clone --recurse-submodules https://github.com/GeoFLAC/DynEarthSol.git
```

## Design notes

* Avoid C++ stream for bulk output, as stream is slower than C-style IO.
* Avoid creating/destroying objects in inner for-loops.
* Avoid static variables and global variables.
* The same sources build `dynearthsol2d` and `dynearthsol3d`, so every change
  must compile with both `make ndims=2` and `make ndims=3`; use `NDIMS`,
  `NODES_PER_ELEM` and `NSTR`, never a literal 2, 3, 4 or 6.
* The code is C++17: the Makefile passes `-std=c++17` and `constants.hpp`
  stops an older standard with a `static_assert`. C++17 features are fine, but
  the code must also compile as C++20, so use nothing C++20 removed or
  deprecates (the `cxx20-compat` workflow checks this with g++ and clang) and
  add no compiler extensions; the variable-length arrays in `bc.cxx`, `ic.cxx`
  and `remeshing.cxx` are the known exceptions. The oldest supported g++ is 8,
  which needs `-lstdc++fs` to link `std::filesystem`; add that to the g++
  branch of the Makefile with the first use of `<filesystem>`.
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
  53: GoSPL step failed
  54: mechanical solve did not converge
  60: Assertion violated    61: unreachable branch

## Source files and classes

The code is procedural: plain data structs such as `Param` and `Variables` go
to free functions, and each subsystem has one `.cxx`/`.hpp` pair. Flat data and
free functions also port to the OpenACC path, where device code cannot make
virtual calls. Keep this shape. Each new file or class is one more
thing that a reader must find and keep in mind. The class rules follow the
Google C++ style guide and Ousterhout's case against "classitis" in *A
Philosophy of Software Design*.

* Put new code in the file of the subsystem that it belongs to. A helper that
  one file uses stays in that file, as a `static` function or in an anonymous
  namespace. A helper that several files use goes next to the data that it
  works on.
* Add a new `.cxx`/`.hpp` pair only for a new subsystem: its own data, its own
  entry points called from the time loop or from several files, and some
  hundreds of lines of code. A few dozen lines go into an existing file.
* Split a file only when it is very large -- several thousand lines -- and only
  along a seam that it already has: a group of functions with its own data and
  its own callers. The split is its own commit that moves the code unchanged,
  so the results do not change; edit the moved code in a later commit.
* Prefer a struct and free functions to a class. Add a class only when it owns
  a resource or keeps an invariant across calls: an open file, a search tree,
  an array layout that it checks. Do not add a class to group functions (a file
  or a namespace does that), to hold the state of one call, or as an interface
  with one implementation. Select among variants with a `switch` on the option
  code, as for `rheol_type`, not with a class hierarchy.

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

## Code comments

A comment is for the reader of the code as it is now. It gives what the code
cannot say by itself: the physics or the method, with its reference; units and
sign conventions; an invariant and what breaks without it; a known limit; a
TODO. The rules follow the Linux kernel coding style (say what the code does
and why, not how) and Wilson et al., *Best Practices for Scientific Computing*
(2014): document purpose, not mechanics.

* **No history.** Do not write what the code was before, when or why it
  changed, which bug or run found a problem, or who wrote it: no "previously",
  "used to", "no longer", "the old code", "since v2.0", dates, hashes, pull
  request numbers or author names. `git log` and `git blame` keep this for
  every line; write it in the commit message when a reader needs it. A history
  note stays in the file after the next change makes it false.
* **Docs too.** `CONTRIBUTING.md`, `README.md`, the `--help` text and
  `examples/defaults.cfg` say how the code works now; what changed goes in
  `CHANGELOG.md` and the commit message.
* **Keep the reason, drop the story.** What a reader needs to change the code
  safely stays in the code; Ousterhout puts it as "comments belong in the code,
  not the commit log". A guard says which case it handles -- "a degenerate
  facet puts the marker outside every element" -- not the run that found it.
* **No test campaigns.** No tables of measured cases, no verification lists, no
  timings. A number in the code gets its source in one clause: a paper, a
  derivation, a benchmark case.
* **One fact, one place.** Write a rule once, at the code that implements it;
  other places give the name of that function. An overview block lists names,
  not their rules.
* Do not restate the code. Do not keep code in comments: delete it, and the
  history keeps it.
* A function header is one to five lines; an inline note is one or two.

## Commit messages

The message records what the code cannot: what was wrong before, why this fix
and not another, and what it was measured to do, only as much as a reader needs
to judge or revert the change. The order in which things were tried, debug
sessions and verification matrices go in the pull request, not in the message.
Conventions follow the Git and Linux kernel patch guides and Chris Beams' seven
rules; the split with code comments follows Hutterer's "On commit messages".

```
type: imperative summary, 72 characters or fewer, no period

The problem, in the present tense, and its user-visible effect. Whether this is
a one-line fix or a new feature, there is a problem that motivated it -- state
it first, so a reader knows whether to read on.

The approach, only where the diff does not make it obvious, and any alternative
considered and rejected. Numbers for every claim: a mesh that changed, a field
that moved, a timing -- give the case and the value.

Co-Authored-By: Name <email>
```

* `type` is one of `fix`, `feat`, `refactor`, `test`, `docs`, `ci`, `build`,
  `cleanup`. A subsystem may follow in parentheses, `fix(rsf):`, when a series
  is confined to one; most commits do not need it. Lowercase after the colon.
* **Comment versus commit.** "Code comments" says what stays in the code. When
  the change's reason lives there -- an ordering contract, a guard's reason --
  the body says so ("asserted at the consumer") rather than repeat it.
* Refer to earlier work by its subject line or pull request number, never by
  hash alone: branches here are rebased and replayed, so hashes go stale.
* Trailers: `Co-Authored-By` for every AI tool or person who wrote part of the
  change; `Fixes #N` when the commit closes an issue.
* Length is set by the problem, not by habit. Most fixes fit three short
  paragraphs.

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

To package an existing tag again, run *Release with submodules and Zenodo* from
the Actions tab with the tag as its input.

A tag name that has ever carried a published release is spent: deleting the
release frees the tag but not the name, and neither does deleting the repository
([immutable releases](https://docs.github.com/en/code-security/concepts/supply-chain-security/immutable-releases)).
A release that has to be redone therefore takes the next version number.

## Git Submodule Workflow

This repository uses three Git submodules: `knn-bvh`, `mmg` and `nanoflann`.
A submodule is a pointer to one commit of another repository, so your workflow
must account for them.

**1. Pulling the latest changes:**
When you pull updates from the main repository, Git does not automatically update the submodule contents by default. To fetch and update everything in one command, use:
```bash
git pull --recurse-submodules
```
*(If you simply ran `git pull` and notice your submodules are out of sync, you can fix it by running `git submodule update`.)*

**2. Modifying a submodule:**
By default, submodules are in a "detached HEAD" state. Only `knn-bvh` is a
GeoFLAC repository; a change to `mmg` or `nanoflann` goes to its upstream
project first, and this repository then moves the pointer (step 5). To modify
`knn-bvh`:
1. `cd` into the submodule directory.
2. Checkout the appropriate branch (e.g., `git checkout main`).
3. Make your changes, `git add`, `git commit`, and `git push` from **within** the submodule directory.
4. `cd` back to the main repository root. The main repo will now see that the submodule pointer has changed. 
5. `git add` the submodule directory and commit this pointer update to the main repository.
