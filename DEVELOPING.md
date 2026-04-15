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
* Meaning of error codes:
   1: User input error
   2: IO error
  10: Triangulation/tetrahedralization error
  11: Runtime error
  12: Assertion error (due to programming)

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
