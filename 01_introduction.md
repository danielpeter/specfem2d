**Table of Contents**

-   [Introduction](#introduction)
    -   [Citation](#citation)
    -   [Support](#support)

Introduction
============

SPECFEM2D allows users to perform 2D and 2.5D (i.e., axisymmetric) simulations of acoustic, elastic, viscoelastic, and poroelastic seismic wave propagation as well as full waveform imaging (FWI) or adjoint tomography.

In fluids, SPECFEM2D uses the classical linearized Euler equation; thus if you have sharp local variations of density in the fluid (highly heterogeneous fluids in terms of density) or if density becomes extremely small in some regions of your model (e.g. for upper-atmosphere studies), before using the code please make sure the linearized Euler equation is a valid approximation in the case you want to study. For more details on that see e.g. Jensen et al. (2011).

The 2D spectral-element solver accommodates regular and unstructured meshes, generated for example by Cubit , Gmsh or GiD . Even mesh creation packages that generate triangles, for instance Delaunay-Voronoi triangulation codes, can be used because each triangle can then easily be decomposed into three quadrangles by linking the barycenter to the center of each edge; while this approach does not generate quadrangles of optimal quality, it can ease mesh creation in some situations and it has been shown that the spectral-element method can very accurately handle distorted mesh elements.

With version 7.0, the 2D spectral-element solver accommodates Convolution PML absorbing layers and well as higher-order time schemes (4th order Runge-Kutta and LDDRK4-6). Convolution or Auxiliary Differential Equation Perfectly Matched absorbing Layers (C-PML or ADE-PML) are described in Martin, Komatitsch, and Ezziani (2008; Martin, Komatitsch, and Gedney 2008; Martin and Komatitsch 2009; Martin et al. 2010; Komatitsch and Martin 2007).

The solver has adjoint capabilities and can calculate finite-frequency sensitivity kernels (Tromp, Komatitsch, and Liu 2008; Peter et al. 2011) for acoustic, (an)elastic, and poroelastic media. The package also considers 2D SH and P-SV wave propagation. Finally, the solver can run both in serial and in parallel. See SPECFEM2D for the source code.

The SEM is a continuous Galerkin technique (Tromp, Komatitsch, and Liu 2008; Peter et al. 2011), which can easily be made discontinuous (Bernardi, Maday, and Patera 1994; Chaljub 2000; Kopriva, Woodruff, and Hussaini 2002; Chaljub, Capdeville, and Vilotte 2003; Legay, Wang, and Belytschko 2005; Kopriva 2006; Wilcox et al. 2010; Acosta Minolia and Kopriva 2011); it is then close to a particular case of the discontinuous Galerkin technique (Reed and Hill 1973; Lesaint and Raviart 1974; Arnold 1982; Johnson and Pitkäranta 1986; Bourdel, Mazet, and Helluy 1991; Falk and Richter 1999; Hu, Hussaini, and Rasetarinera 1999; Cockburn, Karniadakis, and Shu 2000; Giraldo, Hesthaven, and Warburton 2002; Rivière and Wheeler 2003; Monk and Richter 2005; Grote, Schneebeli, and Schötzau 2006; Ainsworth, Monk, and Muniz 2006; Bernacki, Lanteri, and Piperno 2006; Dumbser and Käser 2006; De Basabe, Sen, and Wheeler 2008; Puente, Ampuero, and Käser 2009; Wilcox et al. 2010; <span>De Basabe</span> and Sen 2010; Étienne et al. 2010), with optimized efficiency because of its tensorized basis functions (Wilcox et al. 2010; Acosta Minolia and Kopriva 2011). In particular, it can accurately handle very distorted mesh elements (Oliveira and Seriani 2011).

It has very good accuracy and convergence properties (Maday and Patera 1989; Seriani and Priolo 1994; Deville, Fischer, and Mund 2002; Cohen 2002; <span>De Basabe</span> and Sen 2007; Seriani and Oliveira 2008; Ainsworth and Wajid 2009; Ainsworth and Wajid 2010; Melvin, Staniforth, and Thuburn 2012). The spectral element approach admits spectral rates of convergence and allows exploiting \(hp\)-convergence schemes. It is also very well suited to parallel implementation on very large supercomputers (Komatitsch et al. 2003; Tsuboi et al. 2003; Komatitsch, Labarta, and Michéa 2008; Carrington et al. 2008; Komatitsch, Vinnik, and Chevrot 2010) as well as on clusters of GPU accelerating graphics cards (Komatitsch 2011; Michéa and Komatitsch 2010; Komatitsch, Michéa, and Erlebacher 2009; Komatitsch et al. 2010). Tensor products inside each element can be optimized to reach very high efficiency (Deville, Fischer, and Mund 2002), and mesh point and element numbering can be optimized to reduce processor cache misses and improve cache reuse (Komatitsch, Labarta, and Michéa 2008). The SEM can also handle triangular (in 2D) or tetrahedral (in 3D) elements (Wingate and Boyd 1996; Taylor and Wingate 2000; Komatitsch et al. 2001; Cohen 2002; Mercerat, Vilotte, and Sánchez-Sesma 2006) as well as mixed meshes, although with increased cost and reduced accuracy in these elements, as in the discontinuous Galerkin method.

Note that in many geological models in the context of seismic wave propagation studies (except for instance for fault dynamic rupture studies, in which very high frequencies or supershear rupture need to be modeled near the fault, see e.g. Benjemaa et al. (2007; Benjemaa et al. 2009; Puente, Ampuero, and Käser 2009; Tago et al. 2010)) a continuous formulation is sufficient because material property contrasts are not drastic and thus conforming mesh doubling bricks can efficiently handle mesh size variations (Komatitsch and Tromp 2002; Komatitsch et al. 2004; Lee et al. 2008; Lee, Chan, et al. 2009; Lee, Komatitsch, et al. 2009).

For a detailed introduction to the SEM as applied to regional seismic wave propagation, please consult Peter et al. (2011; Tromp, Komatitsch, and Liu 2008; Komatitsch and Vilotte 1998; Komatitsch and Tromp 1999; Chaljub et al. 2007) and in particular Lee, Komatitsch, et al. (2009; Lee, Chan, et al. 2009; Lee et al. 2008; Godinho et al. 2009; Wijk et al. 2004; Komatitsch et al. 2004). A detailed theoretical analysis of the dispersion and stability properties of the SEM is available in Cohen (2002), <span>De Basabe</span> and Sen (2007), Seriani and Oliveira (2007), Seriani and Oliveira (2008) and Melvin, Staniforth, and Thuburn (2012).

The SEM was originally developed in computational fluid dynamics (Patera 1984; Maday and Patera 1989) and has been successfully adapted to address problems in seismic wave propagation. Early seismic wave propagation applications of the SEM, utilizing Legendre basis functions and a perfectly diagonal mass matrix, include Cohen, Joly, and Tordjman (1993), Komatitsch (1997), Faccioli et al. (1997), Casadei and Gabellini (1997), Komatitsch and Vilotte (1998) and Komatitsch and Tromp (1999), whereas applications involving Chebyshev basis functions and a non-diagonal mass matrix include Seriani and Priolo (1994), Priolo, Carcione, and Seriani (1994) and Seriani, Priolo, and Pregarz (1995). In the Legendre version that we use in SPECFEM the mass matrix is purposely slightly inexact but diagonal (but can be made exact if needed, see Teukolsky (2015)), while in the Chebyshev version it is exact but non diagonal.

Beware that, in a spectral-element method, some spurious modes (that have some similarities with classical so-called “Hourglass modes” in finite-element techniques, although in the SEM they are not zero-energy modes) can appear in some (but not all) cases in the spectral element in which the source is located. Fortunately, they do not propagate away from the source element. However, this means that if you put a receiver in the same spectral element as a source, the recorded signals may in some cases be wrong, typically exhibiting some spurious oscillations, which are often even non causal. If that is the case, an easy option is to slightly change the mesh in the source region in order to get rid of these Hourglass-like spurious modes, as explained in Duczek et al. (2014), in which this phenomenon is described in details, and in which practical solutions to avoid it are suggested.

All SPECFEM2D software is written in Fortran2003 with full portability in mind, and conforms strictly to the Fortran2003 standard. It uses no obsolete or obsolescent features of Fortran. The package uses parallel programming based upon the Message Passing Interface (MPI) (Gropp, Lusk, and Skjellum 1994; Pacheco 1997).

This new release of the code includes support for GPU graphics card acceleration (Komatitsch 2011; Michéa and Komatitsch 2010; Komatitsch, Michéa, and Erlebacher 2009; Komatitsch et al. 2010).

The code uses the plane strain convention when the standard P-SV equation case is used, i.e., the off-plane strain \(\epsilon_{zz}\) is zero by definition of the plane strain convention but the off-plane stress \(\sigma_{zz}\) is not equal to zero, one has \(\sigma_{zz} = \lambda (\epsilon_{xx} + \epsilon_{yy})\). This implies, as in any plain strain software package, that the P-SV source is a line source along the direction perpendicular to the plane (see file discussion\_of\_2D\_sources\_and\_approximations\_from\_Pilant\_1979.pdf for more details).

Citation
--------

You can find all the references below in format in file `doc/USER_MANUAL/bibliography.bib`.

If you use this code for your own research, please cite at least one article written by the developers of the package, for instance:

-   Tromp, Komatitsch, and Liu (2008),

-   Peter et al. (2011),

-   Vai et al. (1999),

-   Lee, Chan, et al. (2009),

-   Lee et al. (2008),

-   Lee, Komatitsch, et al. (2009),

-   Komatitsch et al. (2010),

-   Komatitsch, Michéa, and Erlebacher (2009),

-   Liu et al. (2004),

-   Chaljub et al. (2007),

-   Komatitsch and Vilotte (1998),

-   Komatitsch and Tromp (1999),

-   Komatitsch et al. (2004),

-   Morency and Tromp (2008),

-   Blanc et al. (2016),

-   and/or other articles from <https://specfem.github.io/komatitsch.free.fr/publications.html>.

If you use the C-PML absorbing layer capabilities of the code, please cite at least one article written by the developers of the package, for instance:

-   Xie et al. (2014),

-   Xie et al. (2016).

If you use the UNDO\_ATTENUATION option of the code in order to produce full anelastic/viscoelastic sensitivity kernels, please cite at least one article written by the developers of the package, for instance (and in particular):

-   Komatitsch et al. (2016).

More generally, if you use the attenuation (anelastic/viscoelastic) capabilities of the code, please cite at least one article written by the developers of the package, for instance:

-   Komatitsch et al. (2016),

-   Blanc et al. (2016).

If you use the kernel capabilities of the code, please cite at least one article written by the developers of the package, for instance:

-   Tromp, Komatitsch, and Liu (2008),

-   Peter et al. (2011),

-   Liu and Tromp (2006),

-   Morency, Luo, and Tromp (2009).

If you use the SCOTCH / CUBIT non-structured capabilities, please cite:

-   Martin et al. (2008).

If you use axisymmetric geometries please also cite:

-   Bottero et al. (2016)

The corresponding BibTeX entries may be found in file `doc/USER_MANUAL/bibliography.bib`.

Support
-------

This material is based upon work supported by the USA National Science Foundation under Grants No. EAR-0406751 and EAR-0711177, by the French CNRS, French Inria Sud-Ouest MAGIQUE-3D, French ANR NUMASIS under Grant No. ANR-05-CIGC-002, and European FP6 Marie Curie International Reintegration Grant No. MIRG-CT-2005-017461. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the USA National Science Foundation, CNRS, Inria, ANR or the European Marie Curie program.

References
----------

Acosta Minolia, Cesar A., and David A. Kopriva. 2011. “Discontinuous Galerkin Spectral Element Approximations on Moving Meshes.” *J. Comput. Phys.* 230 (5): 1876–1902. doi:[10.1016/j.jcp.2010.11.038](http://dx.doi.org/10.1016/j.jcp.2010.11.038).

Ainsworth, M., and H. Wajid. 2009. “Dispersive and Dissipative Behavior of the Spectral Element Method.” *SIAM Journal on Numerical Analysis* 47 (5): 3910–37. doi:[10.1137/080724976](http://dx.doi.org/10.1137/080724976).

———. 2010. “Optimally Blended Spectral-Finite Element Scheme for Wave Propagation and NonStandard Reduced Integration.” *SIAM Journal on Numerical Analysis* 48 (1): 346–71. doi:[10.1137/090754017](http://dx.doi.org/10.1137/090754017).

Ainsworth, M., P. Monk, and W. Muniz. 2006. “Dispersive and Dissipative Properties of Discontinuous Galerkin Finite Element Methods for the Second-Order Wave Equation.” *Journal of Scientific Computing* 27 (1): 5–40. doi:[10.1007/s10915-005-9044-x](http://dx.doi.org/10.1007/s10915-005-9044-x).

Arnold, Douglas N. 1982. “An Interior Penalty Finite Element Method with Discontinuous Elements.” *SIAM Journal on Numerical Analysis* 19 (4): 742–60. doi:[10.1137/0719052](http://dx.doi.org/10.1137/0719052).

Benjemaa, M., N. Glinsky-Olivier, V. M. Cruz-Atienza, and J. Virieux. 2009. “3D Dynamic Rupture Simulation by a Finite Volume Method.” *Geophys. J. Int.* 178 (1): 541–60. doi:[10.1111/j.1365-246X.2009.04088.x](http://dx.doi.org/10.1111/j.1365-246X.2009.04088.x).

Benjemaa, M., N. Glinsky-Olivier, V. M. Cruz-Atienza, J. Virieux, and S. Piperno. 2007. “Dynamic Non-Planar Crack Rupture by a Finite Volume Method.” *Geophys. J. Int.* 171 (1): 271–85. doi:[10.1111/j.1365-246X.2006.03500.x](http://dx.doi.org/10.1111/j.1365-246X.2006.03500.x).

Bernacki, Marc, Stéphane Lanteri, and Serge Piperno. 2006. “Time-Domain Parallel Simulation of Heterogeneous Wave Propagation on Unstructured Grids Using Explicit, Nondiffusive, Discontinuous Galerkin Methods.” *J. Comput. Acoust.* 14 (1): 57–81.

Bernardi, C., Y. Maday, and A. T. Patera. 1994. “A New Nonconforming Approach to Domain Decomposition: The Mortar Element Method.” In *Nonlinear Partial Differential Equations and Their Applications*, edited by H. Brezis and J. L. Lions, 13–51. Séminaires Du Collège de France. Paris: Pitman.

Blanc, Émilie, Dimitri Komatitsch, Emmanuel Chaljub, Bruno Lombard, and Zhinan Xie. 2016. “Highly Accurate Stability-Preserving Optimization of the Zener Viscoelastic Model, with Application to Wave Propagation in the Presence of Strong Attenuation.” *Geophys. J. Int.* 205 (1): 427–39. doi:[10.1093/gji/ggw024](http://dx.doi.org/10.1093/gji/ggw024).

Bottero, Alexis, Paul Cristini, Dimitri Komatitsch, and Mark Asch. 2016. “An Axisymmetric Time-Domain Spectral-Element Method for Full-Wave Simulations: Application to Ocean Acoustics.” *J. Acoust. Soc. Am.* 140 (5): 3520–30. doi:[10.1121/1.4965964](http://dx.doi.org/10.1121/1.4965964).

Bourdel, Françoise, Pierre-Alain Mazet, and Philippe Helluy. 1991. “Resolution of the Non-Stationary or Harmonic Maxwell Equations by a Discontinuous Finite Element Method: Application to an E.M.I. (Electromagnetic Impulse) Case.” In *Proceedings of the 10th International Conference on Computing Methods in Applied Sciences and Engineering*, 405–22. Commack, NY, USA: Nova Science Publishers.

Carrington, Laura, Dimitri Komatitsch, Michael Laurenzano, Mustafa Tikir, David Michéa, Nicolas <span>Le Goff</span>, Allan Snavely, and Jeroen Tromp. 2008. “High-Frequency Simulations of Global Seismic Wave Propagation Using SPECFEM3D\_GLOBE on 62 Thousand Processor Cores.” In *Proceedings of the SC’08 ACM/IEEE Conference on Supercomputing*, 60:1–60:11. Austin, Texas, USA: IEEE Press. doi:[10.1145/1413370.1413432](http://dx.doi.org/10.1145/1413370.1413432).

Casadei, F., and E. Gabellini. 1997. *Implementation of a 3D Coupled Spectral Element Solver for Wave Propagation and Soil-Structure Interaction Simulations*. Ispra, Italy: European Commission Joint Research Center Report EUR17730EN.

Chaljub, E. 2000. “Modélisation Numérique de La Propagation d’ondes Sismiques En Géométrie Sphérique : Application À La Sismologie Globale (Numerical Modeling of the Propagation of Seismic Waves in Spherical Geometry: application to Global Seismology).” PhD thesis, Paris, France: Université Paris VII Denis Diderot.

Chaljub, E., Y. Capdeville, and J. P. Vilotte. 2003. “Solving Elastodynamics in a Fluid-Solid Heterogeneous Sphere: A Parallel Spectral-Element Approximation on Non-Conforming Grids.” *J. Comput. Phys.* 187 (2): 457–91.

Chaljub, Emmanuel, Dimitri Komatitsch, Jean Pierre Vilotte, Yann Capdeville, Bernard Valette, and Gaetano Festa. 2007. “Spectral Element Analysis in Seismology.” In *Advances in Wave Propagation in Heterogeneous Media*, edited by Ru-Shan Wu and Valérie Maupin, 48:365–419. Advances in Geophysics. Elsevier - Academic Press, London, UK.

Cockburn, Bernardo, George E. Karniadakis, and Chi-Wang Shu. 2000. *Discontinuous Galerkin Methods: Theory, Computation and Applications*. Heidelberg, Germany: Springer-Verlag.

Cohen, G., P. Joly, and N. Tordjman. 1993. “Construction and Analysis of Higher-Order Finite Elements with Mass Lumping for the Wave Equation.” In *Proceedings of the Second International Conference on Mathematical and Numerical Aspects of Wave Propagation*, edited by R. Kleinman, 152–60. SIAM, Philadelphia, Pennsylvania, USA.

Cohen, Gary. 2002. *Higher-Order Numerical Methods for Transient Wave Equations*. Berlin, Germany: Springer-Verlag.

De Basabe, J. D., M. K. Sen, and M. F. Wheeler. 2008. “The Interior Penalty Discontinuous Galerkin Method for Elastic Wave Propagation: Grid Dispersion.” *Geophys. J. Int.* 175 (1): 83–93. doi:[10.1111/j.1365-246X.2008.03915.x](http://dx.doi.org/10.1111/j.1365-246X.2008.03915.x).

<span>De Basabe</span>, Jon<span>á</span>s D., and Mrinal K. Sen. 2007. “Grid Dispersion and Stability Criteria of Some Common Finite-Element Methods for Acoustic and Elastic Wave Equations.” *Geophysics* 72 (6): T81–95. doi:[10.1190/1.2785046](http://dx.doi.org/10.1190/1.2785046).

———. 2010. “Stability of the High-Order Finite Elements for Acoustic or Elastic Wave Propagation with High-Order Time Stepping.” *Geophys. J. Int.* 181 (1): 577–90. doi:[10.1111/j.1365-246X.2010.04536.x](http://dx.doi.org/10.1111/j.1365-246X.2010.04536.x).

Deville, M. O., P. F. Fischer, and E. H. Mund. 2002. *High-Order Methods for Incompressible Fluid Flow*. Cambridge, United Kingdom: Cambridge University Press.

Duczek, Sascha, Steffen Liefold, David Schmicker, and Ulrich Gabbert. 2014. “Wave Propagation Analysis Using High-Order Finite Element Methods: Spurious Oscillations Excited by Internal Element Eigenfrequencies.” *Technische Mechanik* 34: 51–71.

Dumbser, M., and M. Käser. 2006. “An Arbitrary High-Order Discontinuous Galerkin Method for Elastic Waves on Unstructured Meshes-II. The Three-Dimensional Isotropic Case.” *Geophys. J. Int.* 167 (1): 319–36. doi:[10.1111/j.1365-246X.2006.03120.x](http://dx.doi.org/10.1111/j.1365-246X.2006.03120.x).

Étienne, V., E. Chaljub, J. Virieux, and N. Glinsky. 2010. “An \(hp\)-Adaptive Discontinuous Galerkin Finite-Element Method for 3-D Elastic Wave Modelling.” *Geophys. J. Int.* 183 (2): 941–62. doi:[10.1111/j.1365-246X.2010.04764.x](http://dx.doi.org/10.1111/j.1365-246X.2010.04764.x).

Faccioli, E., F. Maggio, R. Paolucci, and A. Quarteroni. 1997. “2D and 3D Elastic Wave Propagation by a Pseudo-Spectral Domain Decomposition Method.” *J. Seismol.* 1: 237–51.

Falk, Richard S., and Gerard R. Richter. 1999. “Explicit Finite Element Methods for Symmetric Hyperbolic Equations.” *SIAM Journal on Numerical Analysis* 36 (3): 935–52. doi:[10.1137/S0036142997329463](http://dx.doi.org/10.1137/S0036142997329463).

Giraldo, F. X., J. S. Hesthaven, and T. Warburton. 2002. “Nodal High-Order Discontinuous Galerkin Methods for the Spherical Shallow Water Equations.” *J. Comput. Phys.* 181 (2): 499–525. doi:[10.1006/jcph.2002.7139](http://dx.doi.org/10.1006/jcph.2002.7139).

Godinho, L., P. Amado Mendes, A. Tadeu, A. Cadena-Isaza, C. Smerzini, F. J. Sánchez-Sesma, R. Madec, and Dimitri Komatitsch. 2009. “Numerical Simulation of Ground Rotations Along 2D Topographical Profiles Under the Incidence of Elastic Plane Waves.” *Bull. Seism. Soc. Am.* 99 (2B): 1147–61. doi:[10.1785/0120080096](http://dx.doi.org/10.1785/0120080096).

Gropp, W., E. Lusk, and A. Skjellum. 1994. *Using MPI, Portable Parallel Programming with the Message-Passing Interface*. Cambridge, USA: MIT Press.

Grote, Marcus J., Anna Schneebeli, and Dominik Schötzau. 2006. “Discontinuous Galerkin Finite Element Method for the Wave Equation.” *SIAM Journal on Numerical Analysis* 44 (6): 2408–31. doi:[10.1137/05063194X](http://dx.doi.org/10.1137/05063194X).

Hu, Fang Q., M. Y. Hussaini, and Patrick Rasetarinera. 1999. “An Analysis of the Discontinuous Galerkin Method for Wave Propagation Problems.” *J. Comput. Phys.* 151 (2): 921–46. doi:[10.1006/jcph.1999.6227](http://dx.doi.org/10.1006/jcph.1999.6227).

Jensen, F. B., W. A. Kuperman, M. Porter, and H. Schmidt. 2011. *Computational Ocean Acoustics*. 2nd ed. Berlin, Germany: Springer-Verlag.

Johnson, C., and J. Pitkäranta. 1986. “An Analysis of the Discontinuous Galerkin Method for a Scalar Hyperbolic Equation.” *Math. Comp.* 46: 1–26. doi:[10.1090/S0025-5718-1986-0815828-4](http://dx.doi.org/10.1090/S0025-5718-1986-0815828-4).

Komatitsch, D., and J. Tromp. 1999. “Introduction to the Spectral-Element Method for 3-D Seismic Wave Propagation.” *Geophys. J. Int.* 139 (3): 806–22. doi:[10.1046/j.1365-246x.1999.00967.x](http://dx.doi.org/10.1046/j.1365-246x.1999.00967.x).

———. 2002. “Spectral-Element Simulations of Global Seismic Wave Propagation-I. Validation.” *Geophys. J. Int.* 149 (2): 390–412. doi:[10.1046/j.1365-246X.2002.01653.x](http://dx.doi.org/10.1046/j.1365-246X.2002.01653.x).

Komatitsch, D., and J. P. Vilotte. 1998. “The Spectral-Element Method: An Efficient Tool to Simulate the Seismic Response of 2D and 3D Geological Structures.” *Bull. Seism. Soc. Am.* 88 (2): 368–92.

Komatitsch, D., R. Martin, J. Tromp, M. A. Taylor, and B. A. Wingate. 2001. “Wave Propagation in 2-D Elastic Media Using a Spectral Element Method with Triangles and Quadrangles.” *J. Comput. Acoust.* 9 (2): 703–18. doi:[10.1142/S0218396X01000796](http://dx.doi.org/10.1142/S0218396X01000796).

Komatitsch, D., L. P. Vinnik, and S. Chevrot. 2010. “SHdiff/SVdiff Splitting in an Isotropic Earth.” *J. Geophys. Res.* 115 (B7): B07312. doi:[10.1029/2009JB006795](http://dx.doi.org/10.1029/2009JB006795).

Komatitsch, Dimitri. 1997. “Méthodes Spectrales et Éléments Spectraux Pour L’équation de L’élastodynamique 2D et 3D En Milieu Hétérogène (Spectral and Spectral-Element Methods for the 2D and 3D Elastodynamics Equations in Heterogeneous Media).” PhD thesis, Paris, France: Institut de Physique du Globe.

———. 2011. “Fluid-Solid Coupling on a Cluster of GPU Graphics Cards for Seismic Wave Propagation.” *C. R. Acad. Sci., Ser. IIb Mec.* 339: 125–35. doi:[10.1016/j.crme.2010.11.007](http://dx.doi.org/10.1016/j.crme.2010.11.007).

Komatitsch, Dimitri, and Roland Martin. 2007. “An Unsplit Convolutional Perfectly Matched Layer Improved at Grazing Incidence for the Seismic Wave Equation.” *Geophysics* 72 (5): SM155–M167. doi:[10.1190/1.2757586](http://dx.doi.org/10.1190/1.2757586).

Komatitsch, Dimitri, Gordon Erlebacher, Dominik Göddeke, and David Michéa. 2010. “High-Order Finite-Element Seismic Wave Propagation Modeling with MPI on a Large GPU Cluster.” *J. Comput. Phys.* 229 (20): 7692–7714. doi:[10.1016/j.jcp.2010.06.024](http://dx.doi.org/10.1016/j.jcp.2010.06.024).

Komatitsch, Dimitri, Jesús Labarta, and David Michéa. 2008. “A Simulation of Seismic Wave Propagation at High Resolution in the Inner Core of the Earth on 2166 Processors of MareNostrum.” *Lecture Notes in Computer Science* 5336: 364–77.

Komatitsch, Dimitri, Qinya Liu, Jeroen Tromp, Peter Süss, Christiane Stidham, and John H. Shaw. 2004. “Simulations of Ground Motion in the Los Angeles Basin Based Upon the Spectral-Element Method.” *Bull. Seism. Soc. Am.* 94 (1): 187–206. doi:[10.1785/0120030077](http://dx.doi.org/10.1785/0120030077).

Komatitsch, Dimitri, David Michéa, and Gordon Erlebacher. 2009. “Porting a High-Order Finite-Element Earthquake Modeling Application to NVIDIA Graphics Cards Using CUDA.” *Journal of Parallel and Distributed Computing* 69 (5): 451–60. doi:[10.1016/j.jpdc.2009.01.006](http://dx.doi.org/10.1016/j.jpdc.2009.01.006).

Komatitsch, Dimitri, Seiji Tsuboi, Chen Ji, and Jeroen Tromp. 2003. “A 14.6 Billion Degrees of Freedom, 5 Teraflops, 2.5 Terabyte Earthquake Simulation on the Earth Simulator.” In *Proceedings of the SC’03 ACM/IEEE Conference on Supercomputing*, 4–11. Phoenix, Arizona, USA: ACM. doi:[10.1145/1048935.1050155](http://dx.doi.org/10.1145/1048935.1050155).

Komatitsch, Dimitri, Zhinan Xie, Ebru Bozda<span>ğ</span>, Elliott Sales de Andrade, Daniel Peter, Qinya Liu, and Jeroen Tromp. 2016. “Anelastic Sensitivity Kernels with Parsimonious Storage for Adjoint Tomography and Full Waveform Inversion.” *Geophys. J. Int.* 206 (3): 1467–78. doi:[10.1093/gji/ggw224](http://dx.doi.org/10.1093/gji/ggw224).

Kopriva, D. A. 2006. “Metric Identities and the Discontinuous Spectral Element Method on Curvilinear Meshes.” *Journal of Scientific Computing* 26 (3): 301–27. doi:[10.1007/s10915-005-9070-8](http://dx.doi.org/10.1007/s10915-005-9070-8).

Kopriva, D. A., S. L. Woodruff, and M. Y. Hussaini. 2002. “Computation of Electromagnetic Scattering with a Non-Conforming Discontinuous Spectral Element Method.” *Int. J. Numer. Methods Eng.* 53 (1): 105–22. doi:[10.1002/nme.394](http://dx.doi.org/10.1002/nme.394).

Lee, Shiann Jong, Yu Chang Chan, Dimitri Komatitsch, Bor Shouh Huang, and Jeroen Tromp. 2009. “Effects of Realistic Surface Topography on Seismic Ground Motion in the Yangminshan Region of Taiwan Based Upon the Spectral-Element Method and LiDAR DTM.” *Bull. Seism. Soc. Am.* 99 (2A): 681–93. doi:[10.1785/0120080264](http://dx.doi.org/10.1785/0120080264).

Lee, Shiann Jong, How Wei Chen, Qinya Liu, Dimitri Komatitsch, Bor Shouh Huang, and Jeroen Tromp. 2008. “Three-Dimensional Simulations of Seismic Wave Propagation in the Taipei Basin with Realistic Topography Based Upon the Spectral-Element Method.” *Bull. Seism. Soc. Am.* 98 (1): 253–64. doi:[10.1785/0120070033](http://dx.doi.org/10.1785/0120070033).

Lee, Shiann Jong, Dimitri Komatitsch, Bor Shouh Huang, and Jeroen Tromp. 2009. “Effects of Topography on Seismic Wave Propagation: An Example from Northern Taiwan.” *Bull. Seism. Soc. Am.* 99 (1): 314–25. doi:[10.1785/0120080020](http://dx.doi.org/10.1785/0120080020).

Legay, A., H. W. Wang, and T. Belytschko. 2005. “Strong and Weak Arbitrary Discontinuities in Spectral Finite Elements.” *Int. J. Numer. Methods Eng.* 64 (8): 991–1008. doi:[10.1002/nme.1388](http://dx.doi.org/10.1002/nme.1388).

Lesaint, P., and P. A. Raviart. 1974. “On a Finite-Element Method for Solving the Neutron Transport Equation (Proc. Symposium, Mathematical Research Center).” In *Mathematical Aspects of Finite Elements in Partial Differential Equations*, edited by Univ. of Wisconsin-Madison, 33:89–123. New York, USA: Academic Press.

Liu, Qinya, and Jeroen Tromp. 2006. “Finite-Frequency Kernels Based on Adjoint Methods.” *Bull. Seism. Soc. Am.* 96 (6): 2383–97. doi:[10.1785/0120060041](http://dx.doi.org/10.1785/0120060041).

Liu, Qinya, Jascha Polet, Dimitri Komatitsch, and Jeroen Tromp. 2004. “Spectral-Element Moment Tensor Inversions for Earthquakes in Southern California.” *Bull. Seism. Soc. Am.* 94 (5): 1748–61. doi:[10.1785/012004038](http://dx.doi.org/10.1785/012004038).

Maday, Y., and A. T. Patera. 1989. “Spectral-Element Methods for the Incompressible Navier-Stokes Equations.” In *State of the Art Survey in Computational Mechanics*, 71–143.

Martin, Roland, and Dimitri Komatitsch. 2009. “An Unsplit Convolutional Perfectly Matched Layer Technique Improved at Grazing Incidence for the Viscoelastic Wave Equation.” *Geophys. J. Int.* 179 (1): 333–44. doi:[10.1111/j.1365-246X.2009.04278.x](http://dx.doi.org/10.1111/j.1365-246X.2009.04278.x).

Martin, Roland, Dimitri Komatitsch, and Abdelaâziz Ezziani. 2008. “An Unsplit Convolutional Perfectly Matched Layer Improved at Grazing Incidence for Seismic Wave Equation in Poroelastic Media.” *Geophysics* 73 (4): T51–61. doi:[10.1190/1.2939484](http://dx.doi.org/10.1190/1.2939484).

Martin, Roland, Dimitri Komatitsch, and Stephen D. Gedney. 2008. “A Variational Formulation of a Stabilized Unsplit Convolutional Perfectly Matched Layer for the Isotropic or Anisotropic Seismic Wave Equation.” *Comput. Model. Eng. Sci.* 37 (3): 274–304. doi:[10.3970/cmes.2008.037.274](http://dx.doi.org/10.3970/cmes.2008.037.274).

Martin, Roland, Dimitri Komatitsch, Céline Blitz, and Nicolas <span>Le Goff</span>. 2008. “Simulation of Seismic Wave Propagation in an Asteroid Based Upon an Unstructured MPI Spectral-Element Method: Blocking and Non-Blocking Communication Strategies.” *Lecture Notes in Computer Science* 5336: 350–63.

Martin, Roland, Dimitri Komatitsch, Stephen D. Gedney, and Émilien Bruthiaux. 2010. “A High-Order Time and Space Formulation of the Unsplit Perfectly Matched Layer for the Seismic Wave Equation Using Auxiliary Differential Equations (ADE-PML).” *Comput. Model. Eng. Sci.* 56 (1): 17–42. doi:[10.3970/cmes.2010.056.017](http://dx.doi.org/10.3970/cmes.2010.056.017).

Melvin, Thomas, Andrew Staniforth, and John Thuburn. 2012. “Dispersion Analysis of the Spectral-Element Method.” *Quarterly Journal of the Royal Meteorological Society* 138 (668): 1934–47. doi:[10.1002/qj.1906](http://dx.doi.org/10.1002/qj.1906).

Mercerat, E. D., J. P. Vilotte, and F. J. Sánchez-Sesma. 2006. “Triangular Spectral-Element Simulation of Two-Dimensional Elastic Wave Propagation Using Unstructured Triangular Grids.” *Geophys. J. Int.* 166 (2): 679–98.

Michéa, David, and Dimitri Komatitsch. 2010. “Accelerating a 3D Finite-Difference Wave Propagation Code Using GPU Graphics Cards.” *Geophys. J. Int.* 182 (1): 389–402. doi:[10.1111/j.1365-246X.2010.04616.x](http://dx.doi.org/10.1111/j.1365-246X.2010.04616.x).

Monk, Peter, and Gerard R. Richter. 2005. “A Discontinuous Galerkin Method for Linear Symmetric Hyperbolic Systems in Inhomogeneous Media.” *Journal of Scientific Computing* 22-23 (1-3): 443–77. doi:[10.1007/s10915-004-4132-5](http://dx.doi.org/10.1007/s10915-004-4132-5).

Morency, C., and J. Tromp. 2008. “Spectral-Element Simulations of Wave Propagation in Poroelastic Media.” *Geophys. J. Int.* 175: 301–45.

Morency, C., Y. Luo, and J. Tromp. 2009. “Finite-Frequency Kernels for Wave Propagation in Porous Media Based Upon Adjoint Methods.” *Geophys. J. Int.* 179: 1148–68. doi:[10.1111/j.1365-246X.2009.04332](http://dx.doi.org/10.1111/j.1365-246X.2009.04332).

Oliveira, S. P., and G. Seriani. 2011. “Effect of Element Distortion on the Numerical Dispersion of Spectral-Element Methods.” *Communications in Computational Physics* 9 (4): 937–58.

Pacheco, P. S. 1997. *Parallel Programming with MPI*. San Francisco, USA: Morgan Kaufmann Press.

Patera, Anthony T. 1984. “A Spectral Element Method for Fluid Dynamics: Laminar Flow in a Channel Expansion.” *J. Comput. Phys.* 54 (3): 468–88. doi:[10.1016/0021-9991(84)90128-1](http://dx.doi.org/10.1016/0021-9991(84)90128-1).

Peter, Daniel, Dimitri Komatitsch, Yang Luo, Roland Martin, Nicolas <span>Le Goff</span>, Emanuele Casarotti, Pieyre <span>Le Loher</span>, et al. 2011. “Forward and Adjoint Simulations of Seismic Wave Propagation on Fully Unstructured Hexahedral Meshes.” *Geophys. J. Int.* 186 (2): 721–39. doi:[10.1111/j.1365-246X.2011.05044.x](http://dx.doi.org/10.1111/j.1365-246X.2011.05044.x).

Priolo, E., J. M. Carcione, and G. Seriani. 1994. “Numerical Simulation of Interface Waves by High-Order Spectral Modeling Techniques.” *J. Acoust. Soc. Am.* 95 (2): 681–93.

Puente, J. de la, J. P. Ampuero, and M. Käser. 2009. “Dynamic Rupture Modeling on Unstructured Meshes Using a Discontinuous Galerkin Method.” *J. Geophys. Res.* 114: B10302. doi:[10.1029/2008JB006271](http://dx.doi.org/10.1029/2008JB006271).

Reed, W. H., and T. R. Hill. 1973. *Triangular Mesh Methods for the Neutron Transport Equation*. LA-UR-73-479. Los Alamos, USA: Los Alamos Scientific Laboratory.

Rivière, B., and M. F. Wheeler. 2003. “Discontinuous Finite Element Methods for Acoustic and Elastic Wave Problems.” *Contemporary Mathematics* 329: 271–82.

Seriani, G., and S. P. Oliveira. 2008. “Dispersion Analysis of Spectral-Element Methods for Elastic Wave Propagation.” *Wave Motion* 45: 729–44. doi:[10.1016/j.wavemoti.2007.11.007](http://dx.doi.org/10.1016/j.wavemoti.2007.11.007).

Seriani, G., and E. Priolo. 1994. “A Spectral Element Method for Acoustic Wave Simulation in Heterogeneous Media.” *Finite Elements in Analysis and Design* 16: 337–48.

Seriani, G., E. Priolo, and A. Pregarz. 1995. “Modelling Waves in Anisotropic Media by a Spectral Element Method.” In *Proceedings of the Third International Conference on Mathematical and Numerical Aspects of Wave Propagation*, edited by G. Cohen, 289–98. SIAM, Philadephia, PA.

Seriani, G<span>é</span>za, and Saulo P. Oliveira. 2007. “Optimal Blended Spectral-Element Operators for Acoustic Wave Modeling.” *Geophysics* 72 (5): SM95–M106. doi:[10.1190/1.2750715](http://dx.doi.org/10.1190/1.2750715).

Tago, J., V. M. Cruz-Atienza, V. Étienne, J. Virieux, M. Benjemaa, and F. J. Sánchez-Sesma. 2010. “3D Dynamic Rupture with Anelastic Wave Propagation Using an Hp-Adaptive Discontinuous Galerkin Method.” In *Abstract S51A-1915 Presented at 2010 AGU Fall Meeting*. San Francisco, California, USA.

Taylor, M. A., and B. A. Wingate. 2000. “A Generalized Diagonal Mass Matrix Spectral Element Method for Non-Quadrilateral Elements.” *Appl. Num. Math.* 33: 259–65.

Teukolsky, Saul A. 2015. “Short Note on the Mass Matrix for Gauss-Lobatto Grid Points.” *J. Comput. Phys.* 283: 408–13. doi:[10.1016/j.jcp.2014.12.012](http://dx.doi.org/10.1016/j.jcp.2014.12.012).

Tromp, Jeroen, Dimitri Komatitsch, and Qinya Liu. 2008. “Spectral-Element and Adjoint Methods in Seismology.” *Communications in Computational Physics* 3 (1): 1–32.

Tsuboi, Seiji, Dimitri Komatitsch, Chen Ji, and Jeroen Tromp. 2003. “Broadband Modeling of the 2002 Denali Fault Earthquake on the Earth Simulator.” *Phys. Earth Planet. Inter.* 139 (3-4): 305–13. doi:[10.1016/j.pepi.2003.09.012](http://dx.doi.org/10.1016/j.pepi.2003.09.012).

Vai, R., J. M. Castillo-Covarrubias, F. J. Sánchez-Sesma, D. Komatitsch, and J. P. Vilotte. 1999. “Elastic Wave Propagation in an Irregularly Layered Medium.” *Soil Dynamics and Earthquake Engineering* 18 (1): 11–18. doi:[10.1016/S0267-7261(98)00027-X](http://dx.doi.org/10.1016/S0267-7261(98)00027-X).

Wijk, Kasper van, Dimitri Komatitsch, John A. Scales, and Jeroen Tromp. 2004. “Analysis of Strong Scattering at the Micro-Scale.” *J. Acoust. Soc. Am.* 115 (3): 1006–11. doi:[10.1121/1.1647480](http://dx.doi.org/10.1121/1.1647480).

Wilcox, Lucas C., Georg Stadler, Carsten Burstedde, and Omar Ghattas. 2010. “A High-Order Discontinuous Galerkin Method for Wave Propagation Through Coupled Elastic-Acoustic Media.” *J. Comput. Phys.* 229 (24): 9373–96. doi:[10.1016/j.jcp.2010.09.008](http://dx.doi.org/10.1016/j.jcp.2010.09.008).

Wingate, B. A., and J. P. Boyd. 1996. “Spectral Element Methods on Triangles for Geophysical Fluid Dynamics Problems.” In *Proceedings of the Third International Conference on Spectral and High-Order Methods*, edited by A. V. Ilin and L. R. Scott, 305–14. Houston, Texas: Houston J. Mathematics.

Xie, Zhinan, Dimitri Komatitsch, Roland Martin, and René Matzen. 2014. “Improved Forward Wave Propagation and Adjoint-Based Sensitivity Kernel Calculations Using a Numerically Stable Finite-Element PML.” *Geophys. J. Int.* 198 (3): 1714–47. doi:[10.1093/gji/ggu219](http://dx.doi.org/10.1093/gji/ggu219).

Xie, Zhinan, René Matzen, Paul Cristini, Dimitri Komatitsch, and Roland Martin. 2016. “A Perfectly Matched Layer for Fluid-Solid Problems: Application to Ocean-Acoustics Simulations with Solid Ocean Bottoms.” *J. Acoust. Soc. Am.* 140 (1): 165–75. doi:[10.1121/1.4954736](http://dx.doi.org/10.1121/1.4954736).

-----
> This documentation has been automatically generated by [pandoc](http://www.pandoc.org)
> based on the User manual (LaTeX version) in folder doc/USER_MANUAL/
> (Jun 28, 2021)

