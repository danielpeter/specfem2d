# Axisymmetric (2.5D) simulation - homogeneous acoustic halfspace

This example runs an axisymmetric simulation in a homogeneous acoustic halfspace.
The source is a explosion CMT (equivalent to a force point acting as a pressure source) located on the symmetry axis (x==0).

The simulation runs in parallel using 2 MPI processes.

The output seismograms in folder `OUTPUT_FILES/` can be compared against analytical or fully 3D simulations.
In particular, the example `EXAMPLES/applications/homogeneous_acoustic` in SPECFEM3D Cartesian uses the same setup and provides such.

Reference outputs are provided in folder `REF_SEIS/`.

