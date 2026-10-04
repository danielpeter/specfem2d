# Axisymmetric (2.5D) simulation - homogeneous elastic halfspace

This example runs Lamb's problem as an axisymmetric simulation in a homogeneous elastic halfspace.
The source is a single force pointing downwards located at the surface and on the symmetry axis (x==0).

The simulation runs in parallel using 2 MPI processes.

The output seismograms in folder `OUTPUT_FILES/` can be compared against analytical or fully 3D simulations.
In particular, the example `EXAMPLES/applications/homogeneous_halfspace_HEX8_elastic_absorbing_Stacey_5sides`
in SPECFEM3D Cartesian can be used to provide such (after modifying the source accordingly).

Reference outputs are provided in folder `REF_SEIS/`.

