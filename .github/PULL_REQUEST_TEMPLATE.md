<!-- The sections follow what is known to help reviewers: an 80K-PR study found
     that stating the feedback you want and explaining the code changes are the
     two description elements that raise merge rate and reviewer engagement
     (https://arxiv.org/abs/2602.14611), and Google's reviewer guide has the
     reviewer read the description, then the main file, then the rest
     (https://google.github.io/eng-practices/review/reviewer/navigate.html).
     The rest is DES-specific. Keep every answer to a line or two and the whole
     body under about 400 words: a reviewer should finish it in three minutes.
     Detail that wants more space belongs in the commit messages, the docs or a
     design issue, linked from here. Delete what does not apply. -->

## Why

<!-- The problem and how you noticed it. For physics, the equation or the
     reference; for a fix, the wrong behaviour it replaces. -->

Fixes #

## What

<!-- One bullet per change, one or two lines each: what a user notices, not
     how. Prefix one a user must react to -- a renamed parameter, a new exit
     code, a new output file -- with **Breaking.** as the changelog does. If
     tests, cfgs or docs dominate the diff, say how large the production-code
     change is. -->

-

## For the reviewer

- Feedback wanted: <!-- quick look | physics check | full review; plus any open question -->
- Start at: <!-- the file and function that carry the change; the rest follows from it -->
- Check closely: <!-- the riskiest line or assumption: a sign, a threshold, a loop bound -->

## Numerical effect

<!-- CI builds every flag combination and checks restart determinism, but cannot
     tell whether THIS change altered results. If it can, compare a
     benchmarks-cores case before and after (`make set` on the parent commit,
     then `make cmp`) -- see "Benchmarks and regression testing" in
     CONTRIBUTING.md. -->

- Effect: <!-- none (docs/build only) | bit-exact | round-off, below compare.py's 1e-8 | intended, explained above -->
- Reproduce: <!-- cfg and command, and the number a reviewer should see, e.g. "mean T +5.1 K after 20 kyr" -->

## New or retired names

<!-- Only names someone outside the diff must act on: config parameters, make
     knobs, environment variables read, files written, console lines a script
     might scrape, exit codes. A rename is old -> new plus what a reader of the
     old name does. Keys inside a record and internal identifiers belong in the
     docs, not here. "none" is the common answer and a useful one. -->

- Added:
- Renamed or retired:

## Checklist

<!-- Answer with yes or a short note. No task-list checkboxes: GitHub lets
     anyone with write access toggle those in the rendered page. -->

- Follows [CONTRIBUTING.md](../blob/master/CONTRIBUTING.md); builds for `ndims=2` and `ndims=3`:
- Entry under `## [Unreleased]` in [CHANGELOG.md](../blob/master/CHANGELOG.md):
- New parameters registered in `input.cxx` and documented in `examples/defaults.cfg`:
- Test under `tests/functional/` or `benchmarks-cores/`:

## AI disclosure

<!-- Name any AI tools used and what they produced, or write "No AI tools used". -->
