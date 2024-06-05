[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7271913.svg)](https://doi.org/10.5281/zenodo.7271913)

# Documentation  

Additional information about the code can be found at https://code.delmaestro.org

## Introduction 

This webpage contains the details of a worm algorithm path integral quantum Monte Carlo (WA-PIMC) code actively developed in c++ since 2008 in the [Del Maestro group](http://delmaestro.org/adrian) based on:

- T > 0: [M. Boninsegni, N. V. Prokofiev, and B. Svistunov, Phys. Rev. E <b>74</b>, 036701 (2006)](http://link.aps.org/doi/10.1103/PhysRevE.74.036701)
- T = 0: [A. Sarsa, K.E. Schmidt and W. R. Magro, J. Chem. Phys. <b>113</b>, 1366 (2000)](http://aip.scitation.org/doi/abs/10.1063/1.481926)

It can be used to simulate indistinguishable bosons with various types of realistic interactions in one, two and three spatial dimensions. As written, it takes a large number of command line options and allows for the measurement of essentially any physical observable of interest.

The design philosophy included the goal of abstracting the actual implementation of the WA-PIMC method to a kernel that will never need to be touched by the end user.  The code can be easily extended to study a wide variety of situations by including new types of containers, potentials estimators and communicators.

If you have questions, bug reports or plan to use this code for scientific research, please contact me at Adrian.DelMaestro@utk.edu.

## Installation

This program has been successfully compiled and run on both Intel and AMD systems using clang, g++, pathscale and icpc. Before installing, one needs to ensure that all dependencies are met.  We recommend that the required libraries (boost and blitz) are installed in a `local` folder inside your home directory: `$HOME/local`.

## Dependencies 

The code is written in c++ and makes use of both the <a href="https://github.com/blitzpp/blitz">blitz++</a> and <a href="http://www.boost.org/">boost</a> libraries.  You should be able to grab `blitz` from github and compile from source via the instructions below.

We use many of the boost header-only libraries, but two libraries will need to be compiled: boost_program_options and boost_serialization libraries.  Let us assume that you will be installing both blitz and boost in the folder `$HOME/local` using the GNU C++ compiler.  For icpc or clang, the changes should be obvious, and in particular for the Intel compiler you will need to use `intel-linux` as the toolset while for clang you will use `darwin`.

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
If you want to compile for a specific toolset you could add `--with-toolset=clang`.  Now you are ready to install.  Execute

    ```bash
    ./b2 install --prefix=PREFIX --with-program_options --with-serialization cxxflags="-std=c++17" linkflags="-std=c++17"
    ```
    or if you are using the `clang` compiler on mac os

    ```bash
    ./b2 install --prefix=PREFIX --toolset=clang --with-program_options --with-serialization cxxflags="-std=c++17 -stdlib=libc++" linkflags="-std=c++17 -stdlib=libc++" 
    ```

4. If you want to have multiple versions of the library compiled with different compilers you can use the `--layout=versioned` flag above, or you could add `option.set layout : versioned ;` to your `project-config.jam`.  Note: you may have to rename the `$HOME/include/blitz_VER` directory to remove the version number.
5. You should now have a `PREFIX/include` directory containing the header files for `blitz`, `boost` and `random` and your `PREFIX/lib` directory will contain the following files (the `.dylib` files will only appear on Mac OS X)
    ```bash
    libblitz.a   libboost_program_options.a  libblitz.la  libboost_program_options.dylib 
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
- `-D BOLTZMANNONS=1` constrain the simulation to distinguishable quantum particles (i.e. boltzmannons)
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
- `-D CMAKE_CXX_FLAGS="xxx"` to set additional compiler flags for Release and PIGS builds
- `-D CMAKE_CXX_FLAGS_DEBUG="xxx"` to set additional compiler flags for PIMCDebug and PIGSDebug builds
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

      _____    _____   __  __    _____
     |  __ \  |_   _| |  \/  |  / ____|
     | |__) |   | |   | \  / | | |
     |  ___/    | |   | |\/| | | |
     | |       _| |_  | |  | | | |____
     |_|      |_____| |_|  |_|  \_____|

    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Pre-Equilibration Stage.
    [▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇] 	. Diagonal Pre-Equilibration.
    [▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇] 	. Off-Diagonal Pre-Equilibration.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Relax Worm Constant.


    Converged on C0 =  1.00000

    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Equilibration Stage.
    [▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇▇]
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Measurement Stage.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    1 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    2 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    3 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    4 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    5 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    6 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    7 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    8 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #    9 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   10 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   11 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   12 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   13 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   14 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   15 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   16 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   17 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   18 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   19 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Bin #   20 stored to disk.
    [PIMCID: 395a20fc-f1aa-4169-81e0-97a77e251bc4] - Measurement complete.

during the relaxation process where the string following `PIMCID:` is a uuid that uniquely tags your run, and 20 measurements will be output to disk.  To analyze the results the code, you will need to obtain and install the [pimcscripts](https://github.com/DelMaestroGroup/pimcscripts) package via: 


	pip install --upgrade git+git://github.com/DelMaestroGroup/pimcscripts.git#egg=pimcscripts

Which will add the `pimcscripts` library and a number of useful python analysis programs to your path. 

After this has been completed, you can analyze the results of your run via

```bash
pimcave.py OUTPUT/gce-estimator-05.000-008.996-+000.020-0.01000-395a20fc-f1aa-4169-81e0-97a77e251bc4.dat
```

where `395a20fc-f1aa-4169-81e0-97a77e251bc4` needs to be replaced with the unique identifier generated on your machine.  The results should yield something like:

```bash
# PIMCID 395a20fc-f1aa-4169-81e0-97a77e251bc4
# Number Samples     20
K                  315.54530	    12.70789	 4.03
V                 -509.32392	    11.22932	 2.20
V_ext                0.00000	     0.00000	 0.00
V_int             -507.58296	    10.92684	 2.15
E                 -193.77862	    10.05201	 5.19
E_mu              -194.11968	    10.05299	 5.18
K/N                 18.40800	     0.62058	 3.37
V/N                -29.71109	     0.35392	 1.19
E/N                -11.30309	     0.56667	 5.01
N                   17.05300	     0.17947	 1.05
N^2                292.21000	     6.07779	 2.08
density              0.02343	     0.00025	 1.05
us                1290.97997	    22.17522	 1.72
mcsteps            157.20000	     3.86251	 2.46
diagonal             0.64350	     0.01594	 2.48
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

## Support

The development and maintenance of this code has been supported in part by the National Science Foundation under Award Nos. [DMR-1553991](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991) and [DMR-1808440](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1808440) and in part by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Number DE-SC0024333

<img height="200px" src="https://new.nsf.gov/themes/custom/nsf_theme/components/images/logo/logo-desktop.svg"> <img height="200px" src="https://science.osti.gov/assets/img/doe-logos/logo.png">

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

As discussed above, the driver file for this PIMC program is called pdrive.cpp.  It takes a series of command line options, which are used by the Setup class to initialize ConstantParameters, Container, LookupTable and Communicator objects.  Next, a Potential object is created which describes the potential environment (any walls etc.) and the interactions between bosons. A Path object is then instantiated which holds all the details of the actual world lines of the quantum particles. An Action object is created based on the Potential which holds an approximation of the action to be discretized in the path integral decomposition of the partition function. Finally, the main operating object of the program, of type PathIntegralMonteCarlo is created, which requires both the Path and the [Action](https://code.delmaestro.org/classActionBase.html).  This object performs the actual simulation via a series of [Moves](https://code.delmaestro.org/classMoveBase.html), all of which generate trial world line configurations that exactly sample the kinetic part of the density matrix.  All measurements are made via specific [Estimators](https://code.delmaestro.org/classEstimatorBase.html) with the results being output to disk.

The main kernel of this program should remain relatively untouched, as it has been extensively tested and optimized.  Generality can come from modifying just a few things.  For example, in order to implement a new type of measurement, one would need to write a derived [Estimator](https://code.delmaestro.org/classEstimatorBase.html) class along with modifying the Communicator class to define an output path.  New types of particles and external environments can be added by adding new [Potential](https://code.delmaestro.org/classPotentialBase.html) then updating Setup to allow for their specification at the command line.  Finally, radically different systems can be studied by modifying the [Container](https://code.delmaestro.org/classContainer.html) class.
