The main repository is on:
  http://www.bitbucket.org/tan2/dynearthsol3d

To get the most timely progress, using Mercurial:
  hg clone http://www.bitbucket.org/tan2/dynearthsol3d

=============================================================================

To-do list:

High priority:
* Simple benchmarks for rheology
* Local remeshing
* different init_marker_spacing in regions

Low priority:
* Save output as vtk format directly
* Stress BC, esp no-stress sidewall
* Heatflux BC
* Frictional heating
* Adiabatic cooling (adiabatic temperature profile)
* Internal heating
* Volume changed induced stress

=============================================================================

Design notes:

* Avoid C++ stream for bulk ouput, as stream is slower than C-style IO.
* Avoid creating/destroying objects in inner for-loops.
* Avoid static variables and global variables.
* Meaning of error codes:
   1: User input error
   2: IO error
  10: Triangulation/tetrahedralization error
  11: Runtime error
  12: Assertion error (due to programming)

