**Table of Contents**

-   [Getting Started](#getting-started)
    -   [Visualizing the subroutine calling tree of the source code](#visualizing-the-subroutine-calling-tree-of-the-source-code)
    -   [Becoming a developer of the code, or making small modifications in the source code](#becoming-a-developer-of-the-code-or-making-small-modifications-in-the-source-code)

Getting Started
===============

To download the SPECFEM2D software package, type this:

    git clone --recursive --branch devel https://github.com/geodynamics/specfem2d.git

Note: for people who would like to run the package on Windows rather than on Unix machines, you can install Docker or VirtualBox (installing a Linux in VirtualBox in that latter case) and run it easily from inside that.

We recommend that you add `ulimit -S -s unlimited` to your `.bash_profile` file and/or `limit stacksize unlimited` to your `.cshrc` file to suppress any potential limit to the size of the Unix stack.

Then, to configure the software for your system, run the `configure` shell script. This script will attempt to guess the appropriate configuration values for your system. However, at a minimum, it is recommended that you explicitly specify the appropriate command names for your Fortran compiler (another option is to define FC, CC and MPIF90 in your .bash\_profile or your .cshrc file):

        ./configure FC=gfortran CC=gcc

If you want to run in parallel, i.e., using more than one processor core, then you would type

        ./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi

You can replace the GNU compilers above (gfortran and gcc) with other compilers if you want to; for instance for Intel ifort and icc use FC=ifort CC=icc instead.

Before running the `configure` script, you should probably edit file `flags.guess` to make sure that it contains the best compiler options for your system. Known issues or things to check are:

Intel ifort compiler  
See if you need to add `-assume byterecl` for your machine.
In the case of that compiler, we have noticed that initial release versions sometimes have bugs or issues that can lead to wrong results when running the code, thus we *strongly* recommend using a version for which at least one service pack or update has been installed. *In particular, for version 17 of that compiler, users have reported problems (making the code crash at run time) with the `-assume buffered_io` option; if you notice problems, remove that option from file `flags.guess` or change it to `-assume nobuffered_io` and try again.*

IBM compiler  
See if you need to add `-qsave` or `-qnosave` for your machine.

Mac OS  
You will probably need to install `XCODE`.

IBM Blue Gene machines  
Please refer to the manual of SPECFEM3D\_Cartesian, which contains detailed instructions on how to run on Blue Gene.

The SPECFEM2D software package relies on the SCOTCH library to partition meshes. The SCOTCH library (Pellegrini and Roman 1996) provides efficient static mapping, graph and mesh partitioning routines. SCOTCH is a free software package developed by François Pellegrini et al. from LaBRI and Inria in Bordeaux, France, downloadable from the web page <https://gforge.inria.fr/projects/scotch/>. In case no SCOTCH libraries can be found on the system, the configuration will bundle the version provided with the source code for compilation. The path to an existing SCOTCH installation can to be set explicitly with the option `--with-scotch-dir`. Just as an example:

        ./configure FC=ifort MPIFC=mpif90 --with-mpi --with-scotch-dir=/opt/scotch

If you use the Intel ifort compiler to compile the code, we recommend that you use the Intel icc C compiler to compile Scotch, i.e., use:

        ./configure CC=icc FC=ifort MPIFC=mpif90

For further details about the installation of SCOTCH, go to subdirectory `scotch_5.1.11/` and read `INSTALL.txt`. You may want to download more recent versions of SCOTCH in the future from . Support for the METIS graph partitioner has been discontinued because SCOTCH is more recent and performs better.

When compiling the SCOTCH source code, if you get a message such as: “ld: cannot find -lz”, the Zlib compression development library is probably missing on your machine and you will need to install it or ask your system administrator to do so. On Linux machines the package is often called “zlib1g-dev” or similar. (thus “sudo apt-get install zlib1g-dev” would install it)

You may edit the `Makefile` for more specific modifications. Especially, there are several options available:

-   `-DUSE_MPI` compiles with use of an MPI library.

-   `-DUSE_SCOTCH` enables use of graph partitioner SCOTCH.

After these steps, go back to the main directory of SPECFEM2D/ and type

        make

to create all executables which will be placed into the folder `./bin/`.

By default, the solver runs in single precision. This is fine for most application, but if for some reason you want to run the solver in double precision, run the `configure` script with option “`--enable-double-precision`”. Keep in mind that this will of course double total memory size and will also make the solver around 20 to 30% slower on many processors.

If your compiler has problems with the `use mpi` statements that are used in the code, use the script called `replace_use_mpi_with_include_mpif_dot_h.pl` in the root directory to replace all of them with `include mpif.h` automatically.

If you have problems configuring the code on a Cray machine, i.e. for instance if you get an error message from the `configure` script, try exporting these two variables: `MPI_INC=$CRAY_MPICH2_DIR/include and FCLIBS= `, and for more details if needed you can refer to the `utils/Cray_compiler_information` directory. You can also have a look at the configure script called
`utils/Cray_compiler_information/configure_SPECFEM_for_Piz_Daint.bash`.

Visualizing the subroutine calling tree of the source code
----------------------------------------------------------

Packages such as `doxywizard` can be used to visualize the subroutine calling tree of the source code. `Doxywizard` is a GUI front-end for configuring and running `doxygen`.

Becoming a developer of the code, or making small modifications in the source code
----------------------------------------------------------------------------------

If you want to develop new features in the code, and/or if you want to make small changes, improvements, or bug fixes, you are very welcome to contribute.

To do so, i.e. to access the development branch of the source code with read/write access (in a safe way, no need to worry too much about breaking the package, there is a robot called BuildBot that is in charge of checking and validating all new contributions and changes), please visit this Web page:
<https://github.com/geodynamics/specfem2d/wiki>

To visualize the call tree (calling tree) of the source code, you can see the Doxygen tool available in directory `doc/call_trees_of_the_source_code`.

References
----------

Pellegrini, F., and J. Roman. 1996. “SCOTCH: A Software Package for Static Mapping by Dual Recursive Bipartitioning of Process and Architecture Graphs.” *Lecture Notes in Computer Science* 1067: 493–98.

-----
> This documentation has been automatically generated by [pandoc](http://www.pandoc.org)
> based on the User manual (LaTeX version) in folder doc/USER_MANUAL/
> (Oct  1, 2020)

