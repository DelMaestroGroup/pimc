# Documentation  

## Introduction 

This webpage contains the details of a worm algorithm path integral quantum Monte Carlo (WA-PIMC) code actively developed in c++ since 2008 in the [Del Maestro group](http://delmaestro.org/adrian) based on:

- T > 0: [M. Boninsegni, N. V. Prokofiev, and B. Svistunov, Phys. Rev. E <b>74</b>, 036701 (2006)](http://link.aps.org/doi/10.1103/PhysRevE.74.036701)
- T = 0: [A. Sarsa, K.E. Schmidt and W. R. Magro, J. Chem. Phys. <b>113</b>, 1366 (2000)](http://aip.scitation.org/doi/abs/10.1063/1.481926)

It can be used to simulate indistinguishable bosons with various types of realistic interactions in one, two and three spatial dimensions. As written, it takes a large number of command line options and allows for the measurement of essentially any physical observable of interest.

The design philosophy included the goal of abstracting the actual implementation of the WA-PIMC method to a kernel that will never need to be touched by the end user.  The code can be easily extended to study a wide variety of situations by including new types of containers, potentials estimators and communicators.

If you have questions, bug reports or plan to use this code for scientific research, please contact me at Adrian.DelMaestro@utk.edu.

The development and maintenance of this code base has been supported in part by the National Science Foundation under Award Nos. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991) and [DMR-1808440](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1808440).

[<img width="100" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)


## Installation

This program has been successfully compiled and run on both Intel and AMD systems using clang, g++, pathscale and icpc. Before installing, one needs to ensure that all dependencies are met.  We recommend that the required libraries (boost and blitz) are installed in a `local` folder inside your home directory: `$HOME/local`.

## Dependencies 

The code is written in c++ and makes use of both the <a href="https://github.com/blitzpp/blitz">blitz++</a> and <a href="http://www.boost.org/">boost</a> libraries.  You should be able to grab `blitz` from github and compile from source via the instructions below.

We use many of the boost header-only libraries, but two libraries will need to be compiled: boost_program_options and boost_filesystem libraries.  Let us assume that you will be installing both blitz and boost in the folder `$HOME/local` using the GNU C++ compiler.  For icpc or clang, the changes should be obvious, and in particular for the Intel compiler you will need to use `intel-linux` as the toolset while for clang you will use `darwin`.

If you don't have a `$HOME/local` you should create this directory now via

```bash
mkdir $HOME/local
```

### Blitz ###

Unless you need to use the blitz++'s internal debug functionality initiated through \#`define BZ_DEBUG` which is set by including `debug=1` when compiling, blitz can be used as a 'header only' library and does not need to be compiled.  This is the most common use case.  However, as it doesn't take very long to compile one can proceed as follows:

1. Move into your source directory (create if necessary).
```bash
cd $HOME/local/src
```
2. Clone the latest version of blitz++ from  [github](https://github.com/blitzpp/blitz) into `$HOME/local/src`
3. Move into the blitz source directory
4. Read the instructions in the `INSTALL` file to determine if there is anything special you need to do on your system.
5. Execute

```bash
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=PREFIX ..
make lib
make install
```

where `PREFIX` is the location you want to install the libraries, we suggest `$HOME/local` where `$HOME` is your expanded home directory.

<!-- *Note:* If attempting to compile the old version of blitz-0.9 with gcc version 4.3 or later you may encounter errors
when attempting to build blitz++.  To fix this, before issuing `make lib` and/or
`make install` one needs to add headers to a couple of files.  Move to
`$HOME/local/src/blitz-0.9/blitz` (or similarly, `PREFIX/src/blitz-0.9/blitz`) and
add the line

    #include <cstdlib>

to the top of the files `funcs.h` and `mathfunc.h` and save.  -->

### Boost ###

For detailed instructions on installing boost with compiled libraries please see <a href="https://www.boost.org/more/getting_started/index.html">Section 5.2</a> of the official Boost documentation.

1. Download and decompress boost into `$HOME/local/src/`
2. Change to the boost source directory
3. Execute `bootstrap.sh`
If you want to compile for a specific toolset you could add `--with-toolset=gcc`.  Now you are ready to install.  Execute

    ```bash
    ./b2 install --prefix=PREFIX --with-program_options --with-filesystem --with-system --with-serialization cxxflags="-std=c++14" linkflags="-std=c++14"
    ```
    or if you are using the `clang` compiler on mac os

    ```bash
    ./b2 install --prefix=PREFIX --toolset=darwin --with-program_options --with-filesystem --with-system --with-serialization cxxflags="-std=c++14 -stdlib=libc++" linkflags="-std=c++14 -stdlib=libc++" 
    ```

4. If you want to have multiple versions of the library compiled with different compilers you can use the `--layout=versioned` flag above, or you could add `option.set layout : versioned ;` to your `project-config.jam`.  Note: you may have to rename the `$HOME/include/blitz_VER` directory to remove the version number.
5. You should now have a `PREFIX/include` directory containing the header files for `blitz`, `boost` and `random` and your `PREFIX/lib` directory will contain the following files (the `.dylib` files will only appear on Mac OS X)
    ```bash
    libblitz.a   libboost_filesystem.a      libboost_program_options.a     libboost_system.a 
    libblitz.la  libboost_filesystem.dylib  libboost_program_options.dylib libboost_system.dylib
    ```

6. Update the `LD_LIBRARY_PATH` (or `DYLD_LIBRARY_PATH` on mac os) variable inside your `.bahsrc` or `.bash_profile` to include `PREFIX/lib` e.g.

    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:PREFIX/lib
    ```

7. Source your `.bashrc` or `.bash_profile`.

    ```bash
    source ~/.bashrc
    ```

## Path Integral Monte Carlo

After successfully installing blitz and boost you are now ready to compile the
main pimc program on your system.
PIMC uses CMake for build, test and installation automation. For details on using CMake consult https://cmake.org/documentation/. In short, the following steps should work on UNIX-like systems:

  ```
  mkdir build
  cd build
  cmake ../
  make
  sudo make install
  ```

On Windows try:

  ```
  md build
  cd build
  cmake ../
  cmake --build . --config Release
  cmake --build . --target install
  ```

As above, and with further details below, but you should consider using the following CMake options with the appropriate value instead of xxx :

- `-D NDIM=1|2|3` the number of spatial dimensions
- `-D GPU_BLOCK_SIZE=xxx` equal to the maximum threadblock size and enables GPU acceleration (using [AMD's HIP language](https://github.com/ROCm-Developer-Tools/HIP))
- `-D MAX_GPU_STREAMS=xxx` equal to maximum number of concurrent streams on GPU device
- `-D CMAKE_C_COMPILER=xxx` equal to the name of the C99 Compiler you wish to use (or the environment variable `CC`)
- `-D CMAKE_CXX_COMPILER=xxx` equal to the name of the C++17 compiler you wish to use (or the environment variable `CXX`)
- `-D CMAKE_PREFIX_PATH=xxx` to add a non-standard location for CMake to search for libraries, headers or programs
- `-D CMAKE_INSTALL_PREFIX=xxx` to install pimc to a non-standard location
- `-D BOOST_ROOT=xxx` to add non-standard location for Boost install
- `-D BLITZ_ROOT=xxx` to add non-standard location for Blitz++ install
- `-D STATIC=1` to enable a static build
- `-D CMAKE_BUILD_TYPE=Debug` to build pimc in debug mode
- `-D CMAKE_BUILD_TYPE=PIGS` to build pigs
- `-D CMAKE_BUILD_TYPE=PIGSDebug` to build pigs in debug mode
- `-D DEFAULT_CXX_FLAGS="xxx"` to overwrite default compiler flags for Release and PIGS builds (default `-Wall -fno-math-errno -O3`)
- `-D CMAKE_CXX_FLAGS="xxx" to set additional compiler flags for Release and PIGS builds
- `-D CMAKE_CXX_FLAGS_DEBUG="xxx" to set additional compiler flags for PIMCDebug and PIGSDebug builds
- `-E env CXXFLAGS="xxx"` add additional compiler flags
- `-E env LDFLAGS="xxx"` add additional linker flags

Executables will be installed to `CMAKE_INSTALL_PREFIX` location or if the install is skiped will be located in `build/pimc`.
Executables produced are `pimc.e`, `pimcd.e`, `pigs.e`, and `pigsd.e` for `CMAKE_BUILD_TYPE=Release|Debug|PIGS|PIGSDebug` respectively.

If you run into problems, failures with linking etc., common errors may include
not properly setting your `LD_LIBRARY_PATH` or not starting from a clean build
directory (issue `make clean` or `rm -rf ./*` inside the build directory).

## Usage

In order to get a quick idea of the options which the code accepts type:
```bash
pimc.e --help
```

The code requires various combinations of these options to run, and the help message should give
you an idea about which ones are mandatory.

### Quick Start 

If you want to perform a quick test-run for bulk helium you could try something like:
```bash
./pimc.e -T 5 -N 16 -n 0.02198 -t 0.01 -M 8 -C 1.0 -I aziz -X free -E 10000 -S 20 -l 7 -u 0.02 --relax
```

In order for this to work, you will need a folder named `OUTPUT` in the directory where you type the command as it will produce output files in `OUTPUT` that contain all the results of the code.  Each run of the code is associated with a unique identifying integer: the `PIMCID`.  The options used in this demo include a subset of all the possible options:

| Code Option | Description |
| :-----------: | ----------- |
|`T`     |  temperature in kelvin |
|`N`     |  number of particles |
|`n`     |  density in &Aring;<sup>-ndim</sup> (ndim=spatial dimension) |
|`t`     |  the imaginary time step tau |
|`M`     |  number of time slices involved in a bisection move |
|`C`     |  worm prefactor constant |
|`I`     |  interaction potential |
|`X`     |  external potential |
|`E`     |  number of equilibration steps |
|`S`     |  number of production bins to output |
|`l`     |  potential cutoff length in &Aring; |
|`u`     |  chemical potential in kelvin |
|`relax` |  adjust the worm constant to ensure we are in the diagonal ensemble ~75% of the simulation |
|`o`     |  the number of configurations to be stored to disk|
|`p`     |  process or cpu number|
|`R`     |  restart the simulation with a PIMCID|
|`W`     |  the wall clock run limit in hours|
|`s`     |  supply a gce-state-* file to start the simulation from|
|`P`     |  number of imaginary time slices|
|`D`     |  size of the center of mass move in &Aring;|
|`d`     |  size of the single slice displace move in &Aring;|
|`m`     |  mass of the particles in AMU |
|`b`     |  the type of simulation cell|
|`L`     |  linear system size in &Aring;|
|`a`     |  scattering length  in &Aring;|
|`c`     |  strength of the integrated delta function interaction|
|`Lx`     |  linear system size in the x-direction &Aring;|
|`Ly`     |  linear system size in the y-direction in &Aring;|
|`Lz`     |  linear system size in the z-direction in &Aring;|
|`action`     |  the type of effective action used in the simulation |
|`canonical`     |  restrict to the canonical ensemble |
|`window`     |  the particle number window for restricting number fluctuations in the canonical ensemble|
|`imaginary_time_length`  |  the imaginary time extent in K<sup>-1</sup>|
|`wavefunction`  |  the type of trial wavefunction|
|`dimension`  |  output the spatial dimension that the code was compiled with |
|`pigs`     |  perform a simulation at T = 0 K|
|`max_wind`     |  The maximum winding sector to be sampled.  Default=1|
|`staging`     |  Use staging instead of bisection for diagonal updates.|

All options, including lists of possible values and default values can be seen
by using the `--help flag`.

The output of the above command should yield something like:
```bash
  _____    _____   __  __    _____
 |  __ \  |_   _| |  \/  |  / ____|
 | |__) |   | |   | \  / | | |
 |  ___/    | |   | |\/| | | |
 | |       _| |_  | |  | | | |____
 |_|      |_____| |_|  |_|  \_____|

[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Equilibration Stage.
0.72	 1.00000	 0.95000	   19	0.026101
0.72	 0.95000	 0.90250	   18	0.024728
0.59	 0.90250	 0.81225	   18	0.024728
0.60	 0.81225	 0.73102	   19	0.026101
0.75	 0.73102	 0.69447	   16	0.021980
0.63	 0.69447	 0.65975	   17	0.023354
0.74	 0.65975	 0.62676	   20	0.027475
0.81	 0.62676	 0.62676	   18	0.024728
0.72	 0.62676	 0.59542	   17	0.023354
0.77	 0.59542	 0.59542	   18	0.024728
0.74	 0.59542	 0.56565	   17	0.023354
0.79	 0.56565	 0.56565	   19	0.026101
0.72	 0.56565	 0.53737	   13	0.017859
0.73	 0.53737	 0.51050	   17	0.023354
0.75	 0.51050	 0.51050	   15	0.020606
0.74	 0.51050	 0.48498	   17	0.023354
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Measurement Stage.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   1 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   2 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   3 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   4 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   5 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   6 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   7 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   8 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #   9 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  10 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  11 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  12 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  13 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  14 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  15 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  16 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  17 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  18 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  19 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Bin #  20 stored to disk.
[PIMCID: c5555b0b-a259-49bd-a4b1-7a12b8214fd4] - Measurement complete.
```

during the relaxation process where `PIMCID` is a uuid, and 20 measurements will be output to disk.  To analyze the results the code, you will need to obtain a number of python programs located in a `SCRIPTS` directory which can be obtained via:

```bash
svn checkout --username=SVNID http://svn.delmaestro.org/projects/SCRIPTS/ $HOME/local/pimcscripts
```

Which will place them in a folder `pimcscripts` in your `$HOME/local/`
directory.  Many of these depend on some general utility modules that should be
added to this directory on your local machine.

1. Move in to the `pimcscripts` directory
2. Download the relevant scripts (replacing `svnID` with your svn username)

    ```bash
    svn export --username=svnID http://svn.delmaestro.org/pyutils/pyutils.py
    svn export --username=svnID http://svn.delmaestro.org/pyutils/loadgmt.py
    svn export --username=svnID http://svn.delmaestro.org/pyutils/kevent.py
    ```

It may be advantageous to add a new environment variable for the location of
this folder to your `.bashrc` as you will use these scripts extensively.  In
order to take advantage of many of the plotting options you will need to have
various python libraries installed such as
[Matplotlib](http://matplotlib.sourceforge.net/).  For the extra color options
you will need to download and install the gradient files from
[CPT-City](http://soliton.vm.bytemark.co.uk/pub/cpt-city/pkg/)

After this has been completed, you can analyze the results of your run via

```bash
python $HOME/local/pimcsripts/pimcave.py OUTPUT/gce-estimator-05.000-008.996-+000.020-0.01000-c5555b0b-a259-49bd-a4b1-7a12b8214fd4.dat
```

where `c5555b0b-a259-49bd-a4b1-7a12b8214fd4` needs to be replaced with the unique identifier generated on your machine.  The results should yield something like:

```bash
# PIMCID c5555b0b-a259-49bd-a4b1-7a12b8214fd4
# Number Samples     20
K                  342.70210	    16.30687	 4.76
V                 -480.38334	    17.01402	 3.54
E                 -137.68124	    11.58631	 8.42
E_mu              -138.03494	    11.58918	 8.40
K/N                 19.18999	     0.72609	 3.78
V/N                -26.93371	     0.53264	 1.98
E/N                 -7.74372	     0.58986	 7.62
N                   17.68500	     0.32215	 1.82
N^2                315.36300	    10.87876	 3.45
density              0.02429	     0.00044	 1.82
us                1178.62311	    23.20618	 1.97
mcsteps            127.30000	     2.28738	 1.80
diagonal             0.79007	     0.01326	 1.68
```

The basic idea of running the program is that one needs to setup the simulation
cell, by defining either its specific geometry via the size (`L`) flag, or by a
combination of density (`n`) and number of particles (`N`).  At present, two
types of simulation cells are possible, a hypercube in 1,2 or 3 dimensions with
periodic boundary conditions and a cylinder in 3 dimensions, that is obtained
by defining a radius (`r`). One then needs to setup the details of the
simulation, including the temperature (`T`), chemical potential (`u`),
interaction (`I`) and external (`X`) potential.  The simulation details are then
set via the imaginary time step (`t`), worm parameter (`C`) and number of
equilibration (`E`) steps and production output bins (`S`). A more detailed
grasp of all possible program options can be obtained by reading the main
driver file `pdrive.cpp`.

### Output

The results of running the code are a number of data, state and log files that reside in the `OUTPUT` directory.  If the code is run for the cylinder geometry, there will be an additional copy of the files in `OUTPUT/CYLINDER` which contain measurements that have been restricted to some cutoff radius indicated by including the `w` flag when running.  The generic output files are:

| Output File | Description |
| ----------- | ----------- |
|`gce-estimator-T-L-u-t-PIMCID.dat` |  The main estimator file.  Includes binned averages of various non-vector estimators like the energy and density of particles.|
|`gce-log-T-L-u-t-PIMCID.dat` |  The log file, which includes all the details of the simulation (including the command needed to restart it) and details on acceptance and output. |
|`gce-number-T-L-u-t-PIMCID.dat` |  The number probability distribution |
|`gce-obdm-T-L-u-t-PIMCID.dat` |  The one body density matrix |
|`gce-pair-T-L-u-t-PIMCID.dat` | The pair correlation function |
|`gce-pcycle-T-L-u-t-PIMCID.dat` | The permutation cycle distribution |
|`gce-radial-T-L-u-t-PIMCID.dat` | The radial density |
|`gce-state-T-L-u-t-PIMCID.dat` | The state file (used to restart the simulation) |
|`gce-super-T-L-u-t-PIMCID.dat` |  Contains all superfluid estimators |

Each line in either the scalar or vector estimator files contains a bin which is the average of some measurement over a certain number of Monte Carlo steps.  By averaging bins, one can get the final result along with its uncertainty via the variance.

## General Description

A full understanding of this path integral Monte Carlo code requires an
understanding of the WA-PIMC algorithm alluded to in the introduction.  In this
section, we describe the large building blocks of the code.  The actual
specific details of the implementation can be understood by reading the doxygen
documentation included here as well as reading the actual source code.

Any Monte Carlo simulation whether quantum or classical shares a number of
features in common.  Some type of simulation cell is created with a set of
parameters that describe its physical environment.  The initial state of the
system is guessed, and a series of Moves are performed on the constituents
of the system in such a way that detailed balance is maintained.  After some
suitable equilibration period, measurements are made and their results are
stored to disk.

As discussed above, the driver file for this PIMC program is called pdrive.cpp.
It takes a series of command line options, which are used by the Setup class to
initialize ConstantParameters, Container, LookupTable and Communicator objects.
Next, a Potential object is created which describes the potential environment
(any walls etc.) and the interactions between bosons. A Path object is then
instantiated which holds all the details of the actual world lines of the
quantum particles. An Action object is created based on the Potential which
holds an approximation of the action to be discretized in the path integral
decomposition of the partition function. Finally, the main operating object of
the program, of type PathIntegralMonteCarlo is created, which requires both the
Path and the [Action](@ref ActionBase).  This object performs the actual
simulation via a series of [Moves](@ref MoveBase), all of which generate trial
world line configurations that exactly sample the kinetic part of the density
matrix.  All measurements are made via specific [Estimators](@ref EstimatorBase)
with the results being output to disk.

The main kernel of this program should remain relatively untouched, as it has
been extensively tested and optimized.  Generality can come from modifying just
a few things.  For example, in order to implement a new type of measurement,
one would need to write a derived [Estimator](@ref EstimatorBase) class along
with modifying the Communicator class to define an output path.  New types of
particles and external environments can be added by adding new
[Potential](@ref PotentialBase) then updating Setup to allow for their
specification at the command line.  Finally, radically different systems can be
studied by modifying the [Container](@ref Container) class.
