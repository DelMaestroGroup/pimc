Documentation  {#mainpage}
=============

Introduction {#introduction}
============

This program implements the Worm Algorithm Path Integral Quantum Monte Carlo (WA-PIMC) technique
introduced in [M. Boninsegni, N. V. Prokofiev, and B. Svistunov, Phys. Rev. E <b>74</b>, 036701
(2006)](http://link.aps.org/doi/10.1103/PhysRevE.74.036701). It can be used to
simulate indistinguishable bosons with various types of realistic interactions
in one, two and three spatial dimensions. As written, it takes a large number
of command line options and allows for the measurement of essentially any
physical observable of interest. 

The design philosophy included the goal of abstracting the actual implementation of the WA-PIMC
method to a kernel that will never need to be touched by the end user.  The code can be easily
extended to study a wide variety of situations by including new types of containers, potentials
estimators and communicators.

If you have questions, bug reports or plan to use this code for any scientific research, please
contact me at Adrian.DelMaestro@uvm.edu.

Installation {#installation}
============

This program has been successfully compiled and run on both Intel and AMD systems using
g++, pathscale and icpc. Before installing, one needs to ensure that all
dependencies are met.  We recommend that the required libraries (boost and
blitz) are installed in a `local` folder inside your home directory: `~/local`.

Dependencies {#dependencies}
------------

The code is written in C++ and makes use of both the <a
href="http://www.oonumerics.org/blitz/">blitz++</a> and <a
href="http://www.boost.org/">boost</a> libraries.  Unfortunately the blitz++
library hosted on the official website (version 0.9) is somewhat out of date and
a more recent version with bug-fixes is developed using a mercurial repository.
If you don't have <a href="http://mercurial.selenic.com/"
title="mercurial">mercurial</a> installed on your machine you should do that
now.

We use many of the boost header-only libraries, but two libraries will need to
be compiled: boost_program_options and boost_filesystem libraries.  Let us
assume that you will be installing both blitz and boost in the folder
`$HOME/local` using the GNU C++ compiler.  For icpc or pathscale, the changes
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
2. Get the latest version of blitz++ from the mercurial repository
~~~
hg clone http://blitz.hg.sourceforge.net:8000/hgroot/blitz/blitz blitz
~~~
3. Move into the blitz source directory
4. Generate a configure script by running
~~~
autoreconf -fiv
~~~
5. Execute
~~~
./configure cxx=g++ --prefix=PREFIX
make lib
make install
~~~
where `PREFIX` is the location you want to install the libraries, we suggest
`$HOME/local` where `$HOME` is your expanded home directory.

*Note:* If attempting to compile blitz-0.9 with gcc version 4.3 or later you may encounter errors
when attempting to build blitz++.  To fix this, before issuing `make lib` and/or
`make install` one needs to add headers to a couple of files.  Move to
`$HOME/local/src/blitz-0.9/blitz` (or similarly, `PREFIX/src/blitz-0.9/blitz`) and
add the line 

    #include <cstdlib>

to the top of the files `funcs.h` and `mathfunc.h` and save. 

### Boost ###

For detailed instructions on installing boost with compiled libraries please see <a
href="http://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html#or-build-custom-binaries">Section 5.2</a> 
of the official Boost documentation.

1. Download and decompress boost into `~/local/src/`
2. Change to the directory `tools/build/v2/` inside the boost source directory
3. Execute
~~~
./bootstrap.sh --with-toolset=gcc
~~~
4. Move up to the top level of the boost source directory
5. Execute
~~~
tools/build/v2/b2 install --prefix=PREFIX --toolset=gcc --with-program_options --with-filesystem
~~~
The `b2` executable may also be in `tools/build/v2/bin/` depending on your
machine's configuration.  If you would like to compile boost with different
compilers on your machine and would like to enforce a detailed labelling scheme
for the libraries include `--layout=versioned` when calling `b2` above. See 
[here](http://stackoverflow.com/questions/8940249/boost-how-bjam-constructs-a-library-name "Versioned Library Layout")
for more detail.  
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
main pimc program on your system.  there are currently three makefiles
available: `makefile.g++`, `makefile.icc` and `Makefile.path`.  We only include
details of compiling with g++ here.  For other machines, clusters or
architectures, please read the details in the makefiles.

In order to compile with g++:
1. Open up `Makefile.g++`, find the comment:
~~~
# Edit below to include details on your specific host
~~~
and copy and paste the following:
~~~
ifdef target_hostname
OPT = -Wall -O3 -fno-math-errno
CODEDIR = $$HOME/local
CFLAGS  = $(OPT) $(DIM) $(DEBUG) -I$(CODEDIR)/include
LDFLAGS = -L$(CODEDIR)/lib -lboost_program_options -lboost_filesystem -lboost_system
endif
~~~
where `target_hostname` is replaced with a unique identifier for your machine,
let's call it foo. If you want to run blitz in debug mode you will need to
explicitely link to the blitz library with `-lblitz` added to `LDFLAGS` above.
2. Edit the `CODEDIR` variable to point to the location where you have
   installed blitz and boost above.  We suggest `$HOME/local`
3. Edit the `OPT` variable to reflect yoru local compile options.
4. The make process will then take three options:
 - `debug=1` turn on debugging options
 - `ndim=1,2,3` the number of spatial dimensions
 - `foo=1` the hostname 
5. To compile for bulk (3D) helium in production mode on host foo:
~~~
make -f Makefile.g++ ndim=3 foo=1
~~~
which will produce the executable `pimc.e`
 

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
the `PIMCID`.  The options used in this demo include:

| Code Option | Description |
| :-----------: | ----------- |
|`T`     |  temperature in kelvin |
|`N`     |  number of bosons |
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

All options, including lists of possible values and default values can be seen
by using the `--help flag`.

The output of the above command should yield:

    [PIMCID: pimcid] - Equilibration Stage.
    0.57     1.00000         0.90000           16   0.021980
    0.60     0.90000         0.81000           16   0.021980
    0.66     0.81000         0.76950           17   0.023354
    0.79     0.76950         0.76950           16   0.021980
    0.82     0.76950         0.76950           19   0.026101
    0.75     0.76950         0.73103           16   0.021980
    0.70     0.73103         0.69447           17   0.023354
    0.80     0.69447         0.69447           16   0.021980
    0.72     0.69447         0.65975           16   0.021980
    0.76     0.65975         0.65975           18   0.024728
    0.77     0.65975         0.65975           16   0.021980
    0.74     0.65975         0.62676           15   0.020606
    0.69     0.62676         0.59542           14   0.019233
    0.73     0.59542         0.56565           16   0.021980
    0.85     0.56565         0.59394           18   0.024728
    0.81     0.59394         0.59394           15   0.020606


during the relaxation process where `pimcid` will be replaced with an integer
and 20 measurements will be output to disk.  To analyze the results the code
includes a number of python programs located in the `SCRIPTS` directory.  Many
of these depend on some general utility modules that should be added to this
directory on your local machine.   This can be accomplished by checking them
out of the subversion repository.

1. Move in to the `SCRIPTS` directory
2. Download the relevent scripts (replacing `svnID` with your svn username)
~~~
svn export --username=svnID http://svn.delmaestro.org/pyutils/pyutils.py
svn export --username=svnID http://svn.delmaestro.org/pyutils/loadgmt.py
svn export --username=svnID http://svn.delmaestro.org/pyutils/kevent.py
~~~
In order to take advantage of many of the plotting options you will need to
have various python libraries installed such as
[Matplotlib](http://matplotlib.sourceforge.net/).  For the extra color options
you will need to download and install the gradient files from 
[CPT-City](http://soliton.vm.bytemark.co.uk/pub/cpt-city/pkg/)

After this has been completed, you can analyze the results of your run via

    python SCRIPTS/pimcave.py OUTPUT/gce-estimator-05.000-008.996-+000.020-0.01000-pimcid.dat

where `pimcid` needs to be replaced with the unique identifier generated on
your machine.  The results should yield:

    # Number Samples     20
    K                  347.06803        12.39602
    V                 -469.38790        12.51445
    E                 -122.31988         7.88706
    E_mu              -122.67036         7.88832
    K/N                 19.69534         0.50710
    V/N                -26.62204         0.38998
    E/N                 -6.92669         0.41966
    N                   17.52425         0.22106
    N^2                308.90575         7.74563
    density              0.02407         0.00030
    diagonal             0.77606         0.01076


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
[Potential](@ref PotentialBase) then updating Setup to allow for their specification at the
command line.  Finally, radically different systems can be studied by modifying
the [Container](@ref Container) class.

