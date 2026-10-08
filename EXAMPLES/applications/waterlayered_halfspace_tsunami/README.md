# Waterlayered halfspace - tsunami simulation

This example runs a simulation for a CMT source with an integrated Kuepper wavelet as source time function
in a model with an oceanic layer on top of a homogeneous elastic halfspace.
The simulation includes **gravity** effects in both acoustic and elastic domains, which produce a tsunami wave in the ocean layer.

The setup follows the validation example shown in Fig. 2 of Hou et al. (2026),
"Spectral-element simulation of the earthquake-tsunami coupling and bathymetry effects on oceanic wavefields",
GJI, 246, 1-22. [doi.org/10.1093/gji/ggag228](https://doi.org/10.1093/gji/ggag228)

The output seismogram `AA.S0010.BXZ.semd` is the station at the ocean surface with an epicentral distance of 50 km to the source.
This seismogram can be compared against Fig. 5(a) of Hou et al. (2026).
The corresponding station with the same epicentral distance located at the sea floor is `AA.S0005.BXZ.semd`.

Reference outputs are provided in folder `REF_SEIS/`.

