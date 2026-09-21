# Agent instructions for DynEarthSol

Conventions are in [CONTRIBUTING.md](CONTRIBUTING.md); its "AI-assisted
contributions" section is the policy you work under. This file holds only what
you cannot infer from the tree.

- Parameter names come from `examples/defaults.cfg` or `./dynearthsol2d --help`.
  Never guess one.
- The same sources build `dynearthsol2d` and `dynearthsol3d`; every change must
  compile with both `make ndims=2` and `make ndims=3`. Compiling is not
  verification: anything that can alter results gets `make set` on the parent
  commit and `make cmp` in `benchmarks-cores/` before you report it done.
  `benchmarks-cores/test-tiny.cfg` is a one-second smoke test.
- Issues and pull requests use the templates under `.github/` verbatim; never
  bypass them with `gh pr create --fill` or a hand-written body, and never open
  either or reply to a review without the user's explicit instruction.
- Commit messages follow the "Commit messages" section, with a `Co-Authored-By`
  trailer naming you.
