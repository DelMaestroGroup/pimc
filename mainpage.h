/** 
 * @mainpage Documentation
 *
 * @section sec_intro Introduction
 *
 * This program implements the Worm Algorithm Path Integral Quantum Monte Carlo (WA-PIMC) technique
 * introduced in M. Boninsegni, N. V. Prokofiev, and B. Svistunov, Phys. Rev. E <b>74</b>, 036701
 * (2006). It can be used to simulate indistinguishable bosons with various types of realistic
 * interactions in one, two and three spatial dimensions. As written, it takes a large number of
 * command line options and allows for the measurement of essentially any physical observable of
 * interest. 
 *
 * The design philosophy included the goal of abstracting the actual implementation of the WA-PIMC
 * method to a kernel that will never need to be touched by the end user.  The code can be easily
 * extended to study a wide variety of situations by including new types of containers, potentials
 * estimators and communicators.
 *
 * If you have questions, bug reports or plan to use this code for any scientific research, please
 * contact me at adrian __at__ delmaestro __dot__ org.
 *
 * @section sec_install Installation
 *
 * This program has been successfully compiled and run on both Intel and AMD systems using
 * g++, pathscale and icpc. Before installing, one needs to ensure that all dependencies are met.
 *
 * @subsection subsec_depend Dependencies 
 *
 * The code is written in C++ and makes use of both the <a
 * href="http://www.oonumerics.org/blitz/">blitz++</a> and and <a
 * href="http://www.boost.org/">boost</a> libraries which can be downloaded from their respective
 * websites.  Specifically, you will need to compile the boost_program_options library.  Let us
 * assume that you will be installing both blitz and boost in the folder @p $HOME/local using the GNU
 * C++ compiler.  for icpc or pathscale, the changes should be obvious.
 *
 * @subsubsection blitz
 *
 * -# download and decompress blitz++
 * -# move into the source directory
 * -# @code ./configure cxx=g++ --prefix=$home/local
 *  make lib
 *  make install
 *  @endcode
 *
 * @subsubsection boost
 * for detailed instructions on installing boost with compiled libraries please see <a
 * href="http://www.boost.org/doc/libs/1_44_0/more/getting_started/unix-variants.html#prepare-to-use-a-boost-library-binary">here</a>.
 *
 * -# download and decompress boost as well as boost-jam
 * -# move into the boost-jam source directory
 * -# @code ./build.sh gcc @endcode
 * -# move the newly created bjam file in the bin.arch directory to the top level boost directory
 * -# move into the boost directory
 * -# @code ./bjam install --prefix=$home/local --toolset=gcc --with-program_options @endcode
 * -# update your @p ld_library_path variable to include @p $home/local/lib
 *
 * @subsection subsec_pimc Path Integral Monte Carlo
 *
 * After successfully installing blitz and boost you are now ready to compile the main pimc program
 * on your system.  there are currently four makefiles available: @p makefile.g++ ,  @p makefile.CC
 * , @p makefile.icc and @p Makefile.path.  We only include details of compiling with g++ here.  For
 * other machines, clusters or architectures, please read the details in the makefiles.
 *
 * In order to compile with g++:
 * -# open up @p Makefile.g++ and edit the target hostname with the name of your machines, let's
 *  call it foo.
 * -# edit the CODEDIR variable to point to the location where you have installed blitz and boost
 *  above.  We suggest @p $HOME/local
 * -# the make process will then take four options:
 *  - @p debug=1 turn on debugging options
 *  - @p ndim=1,2,3 the number of spatial dimensions
 *  - @p foo=1 the hostname 
 * -# for example, to compile for bulk (3D) helium in production mode on host foo:
 * @code make -f Makefile.g++ ndim=3 foo=1 @endcode
 * which will produce the executable @p pimc.e
 *  
 * @section sec_usage Usage
 *
 * In order to get a quick idea of the options which the code accepts type:
 *
 * @code pimc.e --help @endcode
 *
 * The code requires various combinations of these options to run, and the help message should give
 * you an idea about which ones are mandatory.
 * 
 * @subsection subsec_quickstart Quick Start
 *
 * If you want to perform a quick test-run for bulk helium you could try something like:
 *
 * @code ./pimc.e -T 5 -N 16 -n 0.02198 -t 0.01 -M 8 -C 1.0 -I aziz -X free -E 10000 -S 20 -l 7 -u 0.02 --relax @endcode
 * 
 * In order for this to work, you will need a folder named @p OUTPUT in the directory where you type
 * the command as it will produce output files in @p OUTPUT that contain all the results of the
 * code.  The options used in this demo include:
 * - @p T  temperature in kelvin
 * - @p N  number of bosons
 * - @p n  density in &Aring;<sup>-NDIM</sup> (NDIM=spatial dimension)
 * - @p t  the imaginary time step tau
 * - @p M  number of time slices involved in a bisection move
 * - @p C  worm prefactor constant
 * - @p I  interaction potential
 * - @p X  external potential
 * - @p E  number of equilibration steps
 * - @p S  number of production bins to output
 * - @p l  potential cutoff length in &Aring;
 * - @p u  chemical potential in kelvin
 * - @p relax  perform relaxation in the worm constant C to ensure we are in the diagonal ensemble
 *   ~75% of the simulation
 *
 * All options, including lists of possible values and default values can be seen by using the @p
 * --help flag.
 *
 * The output of the above command should yield:
 *
 * @code
 * [PIMCID: 072366333] - Equilibration Stage.
 * 0.57     1.00000         0.90000           16   0.021980
 * 0.60     0.90000         0.81000           16   0.021980
 * 0.66     0.81000         0.76950           17   0.023354
 * 0.79     0.76950         0.76950           16   0.021980
 * 0.82     0.76950         0.76950           19   0.026101
 * 0.75     0.76950         0.73103           16   0.021980
 * 0.70     0.73103         0.69447           17   0.023354
 * 0.80     0.69447         0.69447           16   0.021980
 * 0.72     0.69447         0.65975           16   0.021980
 * 0.76     0.65975         0.65975           18   0.024728
 * 0.77     0.65975         0.65975           16   0.021980
 * 0.74     0.65975         0.62676           15   0.020606
 * 0.69     0.62676         0.59542           14   0.019233
 * 0.73     0.59542         0.56565           16   0.021980
 * 0.85     0.56565         0.59394           18   0.024728
 * 0.81     0.59394         0.59394           15   0.020606
 * @endcode
 *
 * during the relaxation process and analyzing the results with @code pimcave @endcode will give:
 *
 * @code
 * # Number Samples     20
 * K                  347.06803        12.39602
 * V                 -469.38790        12.51445
 * E                 -122.31988         7.88706
 * E_mu              -122.67036         7.88832
 * K/N                 19.69534         0.50710
 * V/N                -26.62204         0.38998
 * E/N                 -6.92669         0.41966
 * N                   17.52425         0.22106
 * N^2                308.90575         7.74563
 * density              0.02407         0.00030
 * diagonal             0.77606         0.01076
 * @endcode
 *
 * The basic idea of running the program is that one need to setup the simulation cell, by defining
 * either its specific geometry via the size (@p L) flag, or by a combination of density (@p n) and
 * number of particles (@p N).  At present, two types of simulation cells are possible, a hypercube
 * in 1,2 or 3 dimensions with periodic boundary conditions and a cylinder in 3 dimensions, that is
 * obtained by defining a radius (@p r). One then needs to setup the details of the simulation,
 * including the temperature (@p T), chemical potential (@p u), interaction (@p I) and external (@p
 * X) potential.  The simulation details are then set via the imaginary time step (@p t), worm
 * parameter (@p C) and number of equilibration (@p E) steps and production output bins (@p S). A
 * more detailed grasp of all possible program options can be obtained by reading the main driver
 * file @p pdrive.cpp.
 *
 * @section sec_output Output
 *
 * The results of running the code are a number of data, state and log files that reside in the
 * @p OUTPUT directory.  If the code is run for the cylinder geometry, there will be an additional
 * copy of the files in @p OUTPUT/CYLINDER which contain measurements that have been restricted to
 * some cutoff radius indicated by including the @p w flag when running.  The generic output files
 * are:
 *
 * - @p gce-estimator-T-L-u-t-PIMCID.dat:  The main estimator file.  Includes binned averages of
 *   various non-vector estimators like the energy and density of particles.
 * - @p gce-log-T-L-u-t-PIMCID.dat:  The log file, which includes all the details of the simulation
 *   (including the command needed to restart it) and details on acceptance and output.
 * - @p gce-number-T-L-u-t-PIMCID.dat:  The number probability distribution
 * - @p gce-obdm-T-L-u-t-PIMCID.dat:  The one body density matrix
 * - @p gce-pair-T-L-u-t-PIMCID.dat: The pair correlation function
 * - @p gce-pcycle-T-L-u-t-PIMCID.dat: The permutation cycle distribution
 * - @p gce-radial-T-L-u-t-PIMCID.dat: The radial density 
 * - @p gce-state-T-L-u-t-PIMCID.dat: The state file (used to restart the simulation)
 * - @p gce-super-T-L-u-t-PIMCID.dat:  Contains all superfluid estimators
 * - @p gce-worm-T-L-u-t-PIMCID.dat: Contains details on the worm
 *
 * Each line in either the scalar or vector estimator files contains a bin which is the average of
 * some measurement over a certain number of Monte Carlo steps.  By averaging bins, one can get the
 * final result along with its uncertainty via the variance.
 *
 * @section sec_general General Description
 *
 * A full understanding of this path integral Monte Carlo code requires an understanding of the
 * WA-PIMC algorithm alluded to in the introduction.  In this section, we describe the large
 * building blocks of the code.  The actual specific details of the implementation can be understood
 * by reading the doxygen documentation included here as well as reading the actual source code.
 *
 * Any Monte Carlo simulation whether quantum or classical shares a number of features in common.
 * Some type of simulation cell is created with a set of parameters that describe its physical
 * environment.  The initial state of the system is guessed, and a series of @e moves are performed
 * on the constituents of the system in such a way that detailed balance is maintained.  After some
 * suitable equilibration period, measurements are made and their results are stored to disk.
 *
 * As discussed above, the driver file for this PIMC program is called pdrive.cpp.  It takes a series
 * of command line options, then uses them to setup ConstantParameters, Container, LookupTable and
 * Communicator objects.  Next, a Potential object is created which describes the potential
 * environment (any walls etc.) and the interactions between bosons. A Path object is then
 * instantiated which holds all the details of the actual world lines of the quantum particles. An
 * Action object is created based on the Potential which holds an approximation of the action to be
 * discretized in the path integral decomposition of the partition function. Finally, the main
 * operating object of the program, of type PathIntegralMonteCarlo is created, which requires both the
 * Path and the @link ActionBase Action @endlink.  This object performs the actual simulation via a
 * series of @link MoveBase Moves @endlink, all of which generate trial world line configurations that
 * exactly sample the kinetic part of the density matrix.  All measurements are made via specific
 * @link EstimatorBase Estimators @endlink with the results being output to disk.
 *
 * The main kernel of this program should remain relatively untouched, as it has been extensively
 * tested and optimized.  Generality can come from modifying just a few things.  For example, in
 * order to implement a new type of measurement, one would need to write a derived @link
 * EstimatorBase Estimator @endlink class along with modifying the Communicator class to define an
 * output path.  New types of particles and external environments can be added by adding new @link
 * PotentialBase Potentials @endlink then updating @p pdrive.cpp to allow for their specification at
 * the command line.  Finally, radically different systems can be studied by modifying the Container
 * class.
 *
 */
