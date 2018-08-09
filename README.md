Documentation  {#mainpage}
=============

Introduction {#introduction}
============

This webpage contains the details of a worm algorithm path integral quantum
Monte Carlo (WA-PIMC) code actively developed in c++ since 2008 in the [Del Maestro
group](http://delmaestro.org/adrian) based on:

- T>0: [M. Boninsegni, N. V. Prokofiev, and B. Svistunov, Phys. Rev. E <b>74</b>, 036701 (2006)](http://link.aps.org/doi/10.1103/PhysRevE.74.036701)
- T=0: [A. Sarsa, K.E. Schmidt and W. R. Magro, J. Chem. Phys. <b>113</b>, 1366 (2000)] (http://aip.scitation.org/doi/abs/10.1063/1.481926)

It can be used to simulate indistinguishable bosons with various types of
realistic interactions in one, two and three spatial dimensions. As written, it
takes a large number of command line options and allows for the measurement of
essentially any physical observable of interest. 

The design philosophy included the goal of abstracting the actual
implementation of the WA-PIMC method to a kernel that will never need to be
touched by the end user.  The code can be easily extended to study a wide
variety of situations by including new types of containers, potentials
estimators and communicators.

If you have questions, bug reports or plan to use this code for scientific
research, please contact me at Adrian.DelMaestro@uvm.edu.

Since July 2016, the development and maintenance of this code base has been
supported in part by the National Science Foundation under Award No.
DMR-1553991.

![NSF](https://www.nsf.gov/images/logos/nsf1.gif)

Installation {#installation}
============

This program has been successfully compiled and run on both Intel and AMD systems using
g++, pathscale and icpc. Before installing, one needs to ensure that all
dependencies are met.  We recommend that the required libraries (boost and
blitz) are installed in a `local` folder inside your home directory: `~/local`.

Dependencies {#dependencies}
------------

The code is written in c++ and makes use of both the <a
href="http://sourceforge.net/projects/blitz/">blitz++</a> and <a
href="http://www.boost.org/">boost</a> libraries.  You should be able to
download `blitz-0.10` and compile from source via the instructions below.

We use many of the boost header-only libraries, but two libraries will need to
be compiled: boost_program_options and boost_filesystem libraries.  Let us
assume that you will be installing both blitz and boost in the folder
`$HOME/local` using the GNU C++ compiler.  For icpc or clang, the changes
should be obvious, and in particular for the Intel compiler you will need to use
`intel-linux` as the toolset.

If you don't have a `$HOME/local` you should create this directory now via

    mkdir $HOME/local

### Blitz ###

Unless you need to use the blitz++'s internal debug functionality initiated
through \#`define BZ_DEBUG` which is set by including `debug=1` when compiling
the path integral code, blitz can be used as a 'header only' library and does
not need to be compiled.  This is the most common use case.  However, as it
doesn't take very long to compile one can proceed as follows:

1. Move into your source directory (create if necessary).
~~~
cd $HOME/local/src 
~~~
2. Get the latest version of blitz++ from the sourceforge <a
href="http://sourceforge.net/projects/blitz/">website</a>
3. Move into the blitz source directory
4. Read the instructions in the `INSTALL` file to determine if there is
anything special you need to do on your system.
5. Execute
~~~
./configure cxx=g++ --prefix=PREFIX
make lib
make install
~~~
where `PREFIX` is the location you want to install the libraries, we suggest
`$HOME/local` where `$HOME` is your expanded home directory.

*Note:* If attempting to compile the old version of blitz-0.9 with gcc version 4.3 or later you may encounter errors
when attempting to build blitz++.  To fix this, before issuing `make lib` and/or
`make install` one needs to add headers to a couple of files.  Move to
`$HOME/local/src/blitz-0.9/blitz` (or similarly, `PREFIX/src/blitz-0.9/blitz`) and
add the line 

    #include <cstdlib>

to the top of the files `funcs.h` and `mathfunc.h` and save. 

### Boost ###

For detailed instructions on installing boost with compiled libraries please see <a href="https://www.boost.org/more/getting_started/index.html">Section 5.2</a> of the official Boost documentation.

1. Download and decompress boost into `$HOME/local/src/`
2. Change to the boost source directory 
3. Execute
~~~
./bootstrap.sh --with-toolset=gcc
~~~
4. Open `project-config.jam` and make the section:
~~~
if ! gcc in [ feature.values <toolset> ]
{
    using gcc : 5.5.0 : /opt/local/bin/g++-mp-5 ;
}

project : default-build <toolset>gcc ;

~~~
includes all the correct details and location of your gcc compiler, see [here](https://solarianprogrammer.com/2016/03/06/compiling-boost-gcc-5-clang-mac-os-x/) for more details
5. Execute
~~~
./b2 install --prefix=PREFIX --toolset=gcc cxxflags=-std=c++14 linkflags=-std=c++14 --with-program_options --with-filesystem 
~~~
or if you are using the clang compiler
~~~
./b2 install --prefix=PREFIX --toolset=clang --with-program_options --with-filesystem cxxflags="-std=c++14 -stdlib=libc++" linkflags="-std=c++14 -stdlib=libc++" --layout=versioned
~~~
The `b2` executable may also be in `tools/build/bin/` depending on your
machine's configuration.  If you would like to compile boost with different
compilers on your machine and would like to enforce a detailed labelling scheme
for the libraries include `--layout=versioned` when calling `b2` above. See 
[here](http://stackoverflow.com/questions/8940249/boost-how-bjam-constructs-a-library-name "Versioned Library Layout")
for more detail.   You may also need to add `-stdlib=libc++`  to the linkflags when compiling for clang.
6.  You should now have a `PREFIX/include` directory containing the header files
for `blitz`, `boost` and `random` and your `PREFIX/lib` directory will
contain the following files (the `.dylib` files will only appear on Mac OS X)
~~~
libblitz.a   libboost_filesystem.a      libboost_program_options.a libboost_system.a
libblitz.la  libboost_filesystem.dylib  libboost_program_options.dylib libboost_system.dylib
~~~
7. Update the `LD_LIBRARY_PATH` (or `DYLD_LIBRARY_PATH` on Mac OS X) variable
inside your `.bahsrc` or `.bash_profile` to include `PREFIX/lib` eg.  
~~~
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:PREFIX/lib
~~~
8. Source your `.bashrc` or `.bash_profile`.
~~~
source ~/.bashrc
~~~

Path Integral Monte Carlo {#pimc}
-------------------------

After successfully installing blitz and boost you are now ready to compile the
main pimc program on your system.  There are currently four makefiles
available: `Makefile`, `makefile.g++`, `makefile.icc` and `Makefile.path`.  We only include
details of compiling with `g++` here.  For other machines, clusters or
architectures, please read the details in the makefiles.

In order to compile with g++:
1. Open up `Makefile.g++`, find the comment:
~~~
Edit below to include details on your specific host
~~~
and copy and paste the following:
~~~
ifdef target
OPT = -Wall -O3 -fno-math-errno
BOOSTVER =
CODEDIR = $$HOME/local
CFLAGS  = $(OPT) $(DIM) $(DEBUG) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib $(LDEBUG) -lboost_program_options$(BOOSTVER) -lboost_filesystem$(BOOSTVER) -lboost_system$(BOOSTVER)
endif
~~~
where `target` is replaced with a unique identifier for your machine.
If you want to run blitz in debug mode you will need to
explicitly link to the blitz library with `-lblitz` added to `LDFLAGS` above.
2. Edit the `CODEDIR` variable to point to the location where you have
   installed blitz and boost above.  We suggest `$HOME/local`
3. Edit the `OPT` variable to reflect yoru local compile options.
4. If you installed boost with the `--layout=versioned` command above and you
   have multiple versions installed on your machine, you may need to append
   the particular version you want to link to in the names of the boost
   libraries.  This is most easily done by updating the `BOOSTVAR` variable in your
   Makefile: `BOOSTVER = -gcc42-mt-1_49` where here we have compiled boost
   v1.49 with gcc v4.2.  This will need to be updated for your particular
   configuration.  
5. The make process will then take three options:
 - `debug=1` turn on debugging options
 - `ndim=1,2,3` the number of spatial dimensions
 - `target=1` compile for host `target`
5. To compile for bulk (3D) helium in production mode on host target:
~~~
make -f Makefile.g++ ndim=3 target=1
~~~
which will produce the executable `pimc.e`

If you run into problems, failures with linking etc., common errors may include
not properly setting your `LD_LIBRARY_PATH` or not starting from a clean build
directory (issue `make -f Makefile.g++ clean`).

Usage       {#usage}
=====

In order to get a quick idea of the options which the code accepts type:

    pimc.e --help 

The code requires various combinations of these options to run, and the help message should give
you an idea about which ones are mandatory.

Quick Start     {#quickstart}
-----------

If you want to perform a quick test-run for bulk helium you could try something like:

    ./pimc.e -T 5 -N 16 -n 0.02198 -t 0.01 -M 8 -C 1.0 -I aziz -X free -E 10000 -S 20 -l 7 -u 0.02 --relax

In order for this to work, you will need a folder named `OUTPUT` in the directory where you type
the command as it will produce output files in `OUTPUT` that contain all the results of the
code.  Each run of the code is associated with a unique identifying integer:
the `PIMCID`.  The options used in this demo include a subset of all the possible options:

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

The output of the above command should yield:

    [PIMCID: XXXXXXXXX] - Equilibration Stage.
    0.66     1.00000         0.95000           16   0.021980
    0.53     0.95000         0.85500           12   0.016485
    0.56     0.85500         0.76950           14   0.019233
    0.69     0.76950         0.73102           15   0.020606
    0.63     0.73102         0.69447           15   0.020606
    0.69     0.69447         0.65975           17   0.023354
    0.82     0.65975         0.65975           18   0.024728
    0.85     0.65975         0.69274           18   0.024728
    0.75     0.69274         0.65810           16   0.021980
    0.65     0.65810         0.62520           16   0.021980
    0.62     0.62520         0.59394           15   0.020606
    0.63     0.59394         0.56424           14   0.019233
    0.70     0.56424         0.53603           17   0.023354
    0.76     0.53603         0.53603           18   0.024728
    0.75     0.53603         0.50923           16   0.021980
    0.78     0.50923         0.50923           18   0.024728
    [PIMCID: XXXXXXXXX] - Measurement Stage.
    [PIMCID: XXXXXXXXX] - Bin #   1 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   2 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   3 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   4 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   5 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   6 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   7 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   8 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #   9 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  10 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  11 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  12 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  13 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  14 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  15 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  16 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  17 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  18 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  19 stored to disk.
    [PIMCID: XXXXXXXXX] - Bin #  20 stored to disk.
    [PIMCID: XXXXXXXXX] - Measurement complete.

during the relaxation process where `PIMCID` will be replaced with an integer
and 20 measurements will be output to disk.  To analyze the results the code,
you will need to obtain a number of python programs located in a `SCRIPTS`
directory which can be obtained via:

    svn checkout --username=SVNID http://svn.delmaestro.org/projects/SCRIPTS/ $HOME/local/pimcscripts

Which will place them in a folder `pimcscripts` in your `$HOME/local/`
directory.  Many of these depend on some general utility modules that should be
added to this directory on your local machine.   

1. Move in to the `pimcsripts` directory
2. Download the relevant scripts (replacing `svnID` with your svn username)
~~~
svn export --username=svnID http://svn.delmaestro.org/pyutils/pyutils.py
svn export --username=svnID http://svn.delmaestro.org/pyutils/loadgmt.py
svn export --username=svnID http://svn.delmaestro.org/pyutils/kevent.py
~~~
It may be advantageous to add a new environment variable for the location of
this folder to your `.bashrc` as you will use these scripts extensively.  In
order to take advantage of many of the plotting options you will need to have
various python libraries installed such as
[Matplotlib](http://matplotlib.sourceforge.net/).  For the extra color options
you will need to download and install the gradient files from 
[CPT-City](http://soliton.vm.bytemark.co.uk/pub/cpt-city/pkg/)

After this has been completed, you can analyze the results of your run via

    python $HOME/local/pimcsripts/pimcave.py OUTPUT/gce-estimator-05.000-008.996-+000.020-0.01000-XXXXXXXXX.dat

where `XXXXXXXXX` needs to be replaced with the unique identifier generated on
your machine.  The results should yield something like:

    # PIMCID: XXXXXXXXX
    # Number Samples     20
    K                  332.92686	    15.77642
    V                 -438.56116	    19.73930
    E                 -105.63430	    14.58256
    E_mu              -105.97187	    14.58693
    K/N                 19.65034	     0.70454
    V/N                -25.71091	     0.67103
    E/N                 -6.06056	     0.80983
    N                   16.87850	     0.34935
    N^2                287.46950	    11.62243
    density              0.02319	     0.00048
    diagonal             0.77904	     0.01467

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

Output      {#output}
======

The results of running the code are a number of data, state and log files that reside in the
`OUTPUT` directory.  If the code is run for the cylinder geometry, there will be an additional
copy of the files in `OUTPUT/CYLINDER` which contain measurements that have been restricted to
some cutoff radius indicated by including the `w` flag when running.  The generic output files
are:

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
|`gce-worm-T-L-u-t-PIMCID.dat` | Contains details on the worm |

Each line in either the scalar or vector estimator files contains a bin which is the average of
some measurement over a certain number of Monte Carlo steps.  By averaging bins, one can get the
final result along with its uncertainty via the variance.

General Description     {#description}
===================

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

Python Script User Guide {#scripts}
=========================

[PIMC Scripts User Guide](sphinx/index.html)
