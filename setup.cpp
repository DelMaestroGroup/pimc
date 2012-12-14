/** 
 * @file setup.cpp
 * @author Adrian Del Maestro
 * @date 04.06.2011
 * 
 * @brief Implements the setup stage of the simulation.
 */

#include "setup.h"
#include "container.h"
#include "constants.h"
#include "communicator.h"
#include "potential.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SETUP CLASS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
/**************************************************************************//**
 * Setup the program_options variables.
 *
 * We initialize all variables and define the names of all allowed interaction 
 * and external potentials.
******************************************************************************/
Setup::Setup() :
	params(),
    generalOptions("General Options"),
    cellOptions("Simulation Cell Options"),
    potentialOptions("Potential Options"),
    physicalOptions("Physical Options"),
    algorithmicOptions("Algorithmic Options"),
    cmdLineOptions("Combined Command Line Options")
{

	/* Initialize the volume and sides of the simulation cell. */
	side = 0.0;
	volume = 0.0;

	/* An iterator to move through names */
	vector<string>::iterator it;

	/* Define the allowed interaction potential names */
	interactionPotentialName.push_back("aziz");
	interactionPotentialName.push_back("delta");
	interactionPotentialName.push_back("lorentzian");
	interactionPotentialName.push_back("free");
	
	/* Create the interaction potential name string */
	interactionNames = "";
	for (it = interactionPotentialName.begin(); it != interactionPotentialName.end(); ++it) 
		interactionNames += *it + ",";
	interactionNames.erase(interactionNames.end()-1);

	/* Define the allowed interaction potential names */
	externalPotentialName.push_back("free");
	externalPotentialName.push_back("harmonic");
	externalPotentialName.push_back("osc_tube");
	externalPotentialName.push_back("lj_tube");
	externalPotentialName.push_back("hard_tube");
	externalPotentialName.push_back("fixed_aziz");

	/* Create the external potential name string */
	externalNames = "";
	for (it = externalPotentialName.begin(); it != externalPotentialName.end(); ++it)
		externalNames += *it + ",";
	externalNames.erase(externalNames.end()-1);
}

/**************************************************************************//**
 * Define all command line options and get them from the command line.
 *
 * We use boost::program options to get simulation parameters from the 
 * command line.
 * @param argc number of command line arguments
 * @param argv command line string
******************************************************************************/
void Setup::getOptions(int argc, char *argv[])
{
	/* Setup all possible command line options */
	generalOptions.add_options()
		("help,h", "produce help message")
		("output_config,o", po::value<int>()->default_value(0), "number of output configurations")
		("process,p", po::value<uint32>()->default_value(0), "process or cpu number")
		("restart,R", po::value<uint32>()->default_value(0), 
		 "restart running simulation with PIMCID")
		("start_with_state,s", po::value<string>()->default_value(""), 
		 "start simulation with a supplied state file.")
		;

	cellOptions.add_options()
		("dimension", "output currently compiled dimension")
		("cell_type,b", po::value<string>()->default_value("prism"),
		 "simulation cell type [prism,cylinder]")
		("size,L", po::value<double>(), "linear system size [angstroms]")
		("Lx", po::value<double>(), "system size in x-direction [angstroms]")
		("Ly", po::value<double>(), "system size in y-direction [angstroms]")
		("Lz", po::value<double>(), "system size in z-direction [angstroms]")
		("radius,r", po::value<double>(), "tube radius [angstroms]")
		("estimator_radius,w", po::value<double>()->default_value(2.0),
		 "maximum radius for cylinder estimators") 
		;

	potentialOptions.add_options()
		("interaction_potential,I", po::value<string>()->default_value("aziz"), 
		 str(format("interaction potential type {%s}") % interactionNames).c_str())
		("external_potential,X", po::value<string>()->default_value("free"),
		 str(format("external potential type {%s}") % externalNames).c_str())
		("delta_width,a", po::value<double>()->default_value(1.0E-3),
		 "delta function potential width")
		("delta_strength,c", po::value<double>()->default_value(10.0),
		 "delta function potential integrated strength") 
		("fixed,f", po::value<string>()->default_value(""), 
		 "input file name for fixed atomic positions.")
		("potential_cutoff,l", po::value<double>()->default_value(6.0), 
		 "interaction potential cutoff length [angstroms]")
		;

	physicalOptions.add_options()
		("canonical", "perform a canonical simulation")
		("mass,m", po::value<double>()->default_value(4.0030),"particle mass [amu]")
		("density,n", po::value<double>(), str(format("initial density [angstroms^(-%d)]") % NDIM).c_str())
		("number_particles,N", po::value<int>(), "number of particles")
		("temperature,T", po::value<double>(), "temperature [kelvin]")
		("chemical_potential,u", po::value<double>()->default_value(0.0), "chemical potential [kelvin]")
		;

	algorithmicOptions.add_options()
		("relax", "perform a worm constant relaxation")
		("number_time_slices,P", po::value<int>(), "number of time slices")
		("imaginary_time_step,t", po::value<double>(), "imaginary time step [kelvin^(-1)]")
		("worm_constant,C", po::value<double>()->default_value(1.0), "worm acceptance constant")
		("Delta,D", po::value<double>(),"center of mass shift")
		("Mbar,M", po::value<int>(), "worm update length, Mbar")
		("number_eq_steps,E", po::value<uint32>()->default_value(1), 
		 "number of equilibration steps")
		("number_bins_stored,S", po::value<int>()->default_value(1), 
		 "number of estimator bins stored")
		;

	cmdLineOptions.add(generalOptions).add(algorithmicOptions).add(physicalOptions);
	cmdLineOptions.add(cellOptions).add(potentialOptions);

	po::store(po::parse_command_line(argc, argv, cmdLineOptions), params);
	po::notify(params);    

}

/**************************************************************************//**
 * Parse the command line options for obvious errors and return values.
 *
 * Here we go through the commmand line options and test for any problems.  
 * This probably needs more work to test all possible outcomes.
 * @return true if we exit, false if we continue
******************************************************************************/
bool Setup::parseOptions() {

	/* Do we need help? */
	if (params.count("help")) {
		cout << cmdLineOptions << endl;
		return true;
	}

	/* Output the dimension the code was compiled with then exit */
	if (params.count("dimension")) {
		cout << endl << format("Code was compiled for a %d-dimensional system.") % NDIM 
			<< endl << endl;
		return true;
	}

    /* Have we defined a temperature? */
    if (!params.count("temperature")) {
		cerr << endl << "PIMC ERROR: No temperature defined!" << endl << endl;
		cerr << "Action: specify temperature (T)" << endl;
    }

	/* Have we physically defined a simulation cell? */
	definedCell = false;
	if (params.count("size")) {
		definedCell = true;
		side = params["size"].as<double>();
	}
	else if ((params.count("Lx") + params.count("Ly") + params.count("Lz")) == NDIM) {
		definedCell = true;
		if (params.count("Lx"))
			side[0] = params["Lx"].as<double>();
		if (params.count("Ly"))
			side[1] = params["Ly"].as<double>();
		if (params.count("Lz"))
			side[2] = params["Lz"].as<double>();
	}

	/* Make sure we have defined enough options to create the simulation cell */
	if (!( (params.count("density") && params.count("number_particles")) ||
	 	 (definedCell && params.count("number_particles")) ||
		 (definedCell && params.count("density")) ) ) {
		cerr << endl << "PIMC ERROR: Cannot create the simulation cell!" << endl << endl;
		cerr << "Action: define a valid simulation cell." << endl;
		cerr << "Need: [number_particles (N) AND density (n)] OR " << endl;
		cerr << "      [number_particles (N) AND size (L) or Lx,Ly,Lz] OR" << endl;
		cerr << "      [size (L) or Lx,Ly,Lz AND density (n)]" << endl << endl;
		cerr << cellOptions << endl;
		return true;
	}

	/* Make sure we have selected a valid cell type */
	if (!( (params["cell_type"].as<string>() == "cylinder")  || 
	       (params["cell_type"].as<string>() == "prism") ))
	{
		cerr << endl << "PIMC ERROR: Invalid simulation cell type." << endl << endl;
		cerr << "Action: change cell_type (b) to one of:" << endl
			<< "\t[prism,cylinder]" << endl;
		return true;
	}

	/* Can we create the worldlines? */
	if (!(params.count("number_time_slices") || (params.count("imaginary_time_step")))) {
		cerr << endl << "PIMC ERROR: Cannot create imaginary time paths!" << endl << endl;
		cerr << "Action: define number_time_slices (P) OR imaginary_time_step (t)" << endl << endl;
		cerr << algorithmicOptions << endl;
		return true;
	}

	/* Make sure we have selected a valid interaction potential */
	bool validPotential = false;
	vector<string>::iterator it;
	for (it = interactionPotentialName.begin(); it != interactionPotentialName.end(); ++it) {
		if (params["interaction_potential"].as<string>() == *it) {
			validPotential = true;
			break;
		}
	}
	if (!validPotential) {
		cerr << endl << "PIMC ERROR: Invalid interaction potential!" << endl << endl;
		cerr << "Action: set interaction_potential (I) to one of:" << endl
			 << "\t[" << interactionNames << "]" <<  endl;
		return true;
	}

	/* Make sure we have selected a valid external potential */
	validPotential = false;
	for (it = externalPotentialName.begin(); it != externalPotentialName.end(); ++it) {
		if (params["external_potential"].as<string>() == *it) {
			validPotential = true;
			break;
		}
	}
	if (!validPotential) {
		cerr << endl << "PIMC ERROR: Invalid external potential!" << endl << endl;
		cerr << "Action: set external_potential (X) must to one of:" << endl
			 << "\t[" << externalNames << "]" << endl;
		return true;
	}

	/* We can only use the cylinder potentials for a 3D system */
	if ((params["external_potential"].as<string>().find("tube") != string::npos) && (NDIM != 3)) {
		cerr << endl << "PIMC ERROR: Can only use tube potentials for a 3D system!" << endl << endl;
		cerr << "Action: change the potential or recompile with ndim=3." << endl;
		return 1;
	}

	/* Need to specify a radius for the tube potentials */
	if ( (params["external_potential"].as<string>().find("tube") != string::npos) && (!params.count("radius")) ) {
		cerr << endl << "PIMC ERROR: Incomplete specification for external potential!" << endl << endl;
		cerr << "Action: specfity a radius (r) for the tube potentials." << endl;
		return 1;
	}

	return false;
}

/**************************************************************************//**
 * Return the random seed.
 *
 * We add the process number to a fixed initial random seed.
 * @param startSeed The fixed initial seed
 * @return A seed shifted by the process number
******************************************************************************/
uint32 Setup::seed (const uint32 startSeed) {
	return startSeed + params["process"].as<uint32>();
}

/**************************************************************************//**
 * Setup the simulation cell.
 *
 * We setup the simulation cell, and return a pointer to a container opject
 * with a type that depends on the specified simulation cell.
******************************************************************************/
Container * Setup::cell() {

	/* Initialize the local box pointer */
	Container *boxPtr = NULL;

	/* Setup a cylindrical simulation cell */
	if (params["cell_type"].as<string>() == "cylinder") {
		if (definedCell && params.count("number_particles"))
			boxPtr = new Cylinder(params["radius"].as<double>(),side[NDIM-1]);
		else if (definedCell && params.count("density")) {
			boxPtr = new Cylinder(params["radius"].as<double>(),side[NDIM-1]);
			insertOption("number_particles", int(boxPtr->volume*params["density"].as<double>()));
		}
		else
			boxPtr = new Cylinder(params["density"].as<double>(),
					params["radius"].as<double>(),params["number_particles"].as<int>());
	}
	/* Setup a hyperprism */
	else if (params["cell_type"].as<string>() == "prism") {
		if (definedCell && params.count("number_particles")) 
			boxPtr = new Prism(side);
		else if (definedCell && params.count("density")) {
			boxPtr = new Prism(side);
			insertOption("number_particles", int(boxPtr->volume*params["density"].as<double>()));
		}
		else
			boxPtr = new Prism(params["density"].as<double>(),params["number_particles"].as<int>());
	}

	/* Set local copies of box and side */
	volume = boxPtr->volume;
	side = boxPtr->side;

	return boxPtr;
}

/**************************************************************************//**
 * Setup the worldlines.
 *
 * Depending on whether we have defined the size of the imaginary time step
 * tau or the number of time slices we setup the imaginary time extent of the
 * worldlines.
******************************************************************************/
bool Setup::worldlines() {

	/* We determine if we have fixed P or tau.  We require that the number of time 
	 * slices is always even. */
	int numTimeSlices;
	double tau;
	if (!params.count("number_time_slices")) {
		tau = params["imaginary_time_step"].as<double>();
		numTimeSlices = static_cast<int>(1.0/(params["temperature"].as<double>() * tau) + EPS);
		if ((numTimeSlices % 2) != 0)
			numTimeSlices--;
		insertOption("number_time_slices",numTimeSlices);
	}
	else {
		numTimeSlices = params["numTimeSlices"].as<int>();
		if ((numTimeSlices % 2) != 0)
			numTimeSlices--;
		tau = 1.0/(params["temperature"].as<double>() * numTimeSlices);
		setOption("number_time_slices",numTimeSlices);
		insertOption("imaginary_time_step",tau);
	}

	/* If we haven't fixed Mbar, do so now. We use the value defined in
	 * PRE 74, 036701 (2006). We make sure it is always >= 8 */
	int Mbar = 0;
	if (!params.count("Mbar")) {
		/* Compute the average inter-particle separation */
		double ell = pow((1.0*params["number_particles"].as<int>()/volume),-1.0/(1.0*NDIM));
		double lambda = 24.24 / params["mass"].as<double>();
		Mbar = int(ell*ell/(16*lambda*params["imaginary_time_step"].as<double>()));
		if (Mbar < 8)
			Mbar = 8;
		Mbar = 8;
		insertOption("Mbar",Mbar);
	}
	else 
		Mbar = params["Mbar"].as<int>();

	/* Now we make sure it is even and not too large*/
	if (Mbar%2)
		Mbar++;
	if (Mbar > numTimeSlices) {
		cerr << endl << "PIMC ERROR: Update length > number time slices!" << endl << endl;
        cerr << "Action: Increase number_time_slices (P) OR" <<  endl;
        cerr << "        Increase Mbar (M) OR"  << endl;
        cerr << "        Decrease imaginary_time_step (t)" << endl; 
        return true;
	}
	setOption("Mbar",Mbar);

	return false;
}

/**************************************************************************//**
 * Setup the simulation constants.
 *
 * Fix all simulation constants.
******************************************************************************/
void Setup::setConstants() {

	/* Set the required constants */
	constants()->initConstants(
			params.count("canonical"),
			params["temperature"].as<double>(),
			params["chemical_potential"].as<double>(),
			params["mass"].as<double>(),
			params["potential_cutoff"].as<double>(),
			params["worm_constant"].as<double>(),
			volume,
			side[NDIM-1],
			params["number_particles"].as<int>(),
			params["Mbar"].as<int>(),
			params["number_time_slices"].as<int>(),
			params["restart"].as<uint32>(),
			params["process"].as<uint32>(),
			params["number_eq_steps"].as<uint32>(),
			params["interaction_potential"].as<string>(),
			params["external_potential"].as<string>());

	/* If we have specified the center of mass shift delta on the command line, we update
	 * its value. */
	if (params.count("Delta"))
		constants()->setDelta(params["Delta"].as<double>());
	else {
		insertOption("Delta",constants()->Delta());
	}
}

/**************************************************************************//**
 * Setup the communicator.
 *
 * Initialize the communicator, we need to know if we are outputing any config
 * files to disk.  The files are labelled differently depending on whether we
 * are in the canonical or grand-canonical ensemble.  We also need to initialize
 * a possible initial state file and a fixed position file.  Also, since the
 * value of tau we might specifiy at the command line is not the actual one
 * used in the simulation (since the number of time slices must be an integer)
 * we pass it to the communicator for propper labelling of output files.
******************************************************************************/
void Setup::communicator() {
		
    communicate()->init(params["imaginary_time_step"].as<double>(),
            (params["output_config"].as<int>() > 0),params["start_with_state"].as<string>(),
            params["fixed"].as<string>());

    ccommunicate()->init(params["imaginary_time_step"].as<double>(),
            (params["output_config"].as<int>() > 0),params["start_with_state"].as<string>(),
            params["fixed"].as<string>());
}

/*************************************************************************//**
 * Setup the interaction potential.
 *
 * Based on the user's choice we create a new interaction potential pointer
 * which is returned to the main program.  
******************************************************************************/
PotentialBase * Setup::interactionPotential() {

	PotentialBase *interactionPotentialPtr = NULL;
	if (constants()->intPotentialType() == "free")
		interactionPotentialPtr = new FreePotential();
	else if (constants()->intPotentialType() == "delta")
		interactionPotentialPtr = new DeltaPotential(params["delta_width"].as<double>(),
				params["delta_strength"].as<double>());
	else if (constants()->intPotentialType() == "lorentzian")
		interactionPotentialPtr = new LorentzianPotential(params["delta_width"].as<double>(),
				params["delta_strength"].as<double>());
	else if (constants()->intPotentialType() == "aziz")
		interactionPotentialPtr = new AzizPotential(side);

	return interactionPotentialPtr;
}

/*************************************************************************//**
 * Setup the external potential.
 *
 * Based on the user's choice we create a new external potential pointer
 * which is returned to the main program.  
******************************************************************************/
PotentialBase * Setup::externalPotential(const Container* boxPtr) {

	PotentialBase *externalPotentialPtr = NULL;

	if (constants()->extPotentialType() == "harmonic")
		externalPotentialPtr = new HarmonicPotential();
	else if (constants()->extPotentialType() == "free")
		externalPotentialPtr = new FreePotential();
	else if (constants()->extPotentialType() == "osc_tube")
		externalPotentialPtr = new HarmonicCylinderPotential(params["radius"].as<double>());
	else if (constants()->extPotentialType() == "lj_tube")
		externalPotentialPtr = new LJCylinderPotential(params["radius"].as<double>());
	else if (constants()->extPotentialType() == "hard_tube") 
		externalPotentialPtr = new HardCylinderPotential(params["radius"].as<double>());
	else if (constants()->extPotentialType() == "single_well")
		externalPotentialPtr = new SingleWellPotential();
	else if (constants()->extPotentialType() == "fixed_aziz") 
		externalPotentialPtr = new FixedAzizPotential(boxPtr);

	return externalPotentialPtr;
}

/*************************************************************************//**
 * Output the simulation parameters to a log file.
 *
 * After we have finished equilibrating, we output all the simulation
 * parameters to disk in addition to a command that can be used to 
 * restart the simulation.
 * @param argc The number of command line arguments
 * @param argv The commmand line string
 * @param _seed The random seed
 * @param boxPtr A pointer to the container
 * @param nnGrid The lookup table nearest neighbor grid
******************************************************************************/
void Setup::outputOptions(int argc, char *argv[], const uint32 _seed, 
		const Container *boxPtr, const iVec &nnGrid) {

	communicate()->logFile() << endl << "# ";

	/* Construct the command that would be required to restart the simulation */
	bool outputC0 = false;
	bool outputD = false;
	for (int n = 0; n < argc; n++) {

		if ((argv[n][0] == '-') && (argv[n][1] == 's'))
			n++;
		else if ((argv[n][0] == '-') && (argv[n][1] == 'C')) {
			communicate()->logFile() << format("-C %10.4e ") % constants()->C0();
			n++;
			outputC0 = true;
		}
		else if ((argv[n][0] == '-') && (argv[n][1] == 'p')) {
			communicate()->logFile() << format("-p %03d ") % params["process"].as<uint32>();
			n++;
		}
		else if ((argv[n][0] == '-') && (argv[n][1] == 'D')) {
            communicate()->logFile() << format("-D %10.4e ") % constants()->Delta();
			outputD = true;
			n++;
		}
        else if ((argv[n][0] == '-') && (argv[n][1] == '-') && (argv[n][2] == 'r') && 
                (argv[n][3] == 'e') && (argv[n][4] == 'l') && (argv[n][5] == 'a') && 
                (argv[n][6] == 'x')) {
            // Do nothing
        }
		else 
			communicate()->logFile() << argv[n] << " ";
	}

    /* Output the restart flag */
    communicate()->logFile() << format("-R %09d ") % constants()->id();

	/* If we haven't specified the worm constant, output it now */
	if (!outputC0)
		communicate()->logFile() << format("-C %10.4e ") % constants()->C0();

	/* If we haven't specified Delta, output it now */
	if (!outputD)
        communicate()->logFile() << format("-D %10.4e ") % constants()->Delta();

	communicate()->logFile() << endl << endl;
	communicate()->logFile() << "---------- Begin Simulation Parameters ----------" << endl;
	communicate()->logFile() << endl;
	if (constants()->canonical())
		communicate()->logFile() << format("%-24s\t:\t%s\n") % "Ensenble" % "canonical";
	else
		communicate()->logFile() << format("%-24s\t:\t%s\n") % "Ensenble" % "grand canonical";
	communicate()->logFile() << format("%-24s\t:\t%s\n") % "Interaction Potential" % 
		params["interaction_potential"].as<string>();
	if ( (params["interaction_potential"].as<string>().find("delta") != string::npos) ||
			(params["interaction_potential"].as<string>().find("lorentzian") != string::npos) ) {
		communicate()->logFile() << format("%-24s\t:\t%8.3e\n") % "Delta Width" 
			% params["delta_width"].as<double>();
		communicate()->logFile() << format("%-24s\t:\t%-7.2f\n") % "Delta Strength" 
			% params["delta_strength"].as<double>();
	}
	communicate()->logFile() << format("%-24s\t:\t%s\n") 
		% "External Potential" % params["external_potential"].as<string>();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Temperature" % params["temperature"].as<double>();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Chemical Potential" % constants()->mu();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Particle Mass" % params["mass"].as<double>();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Number Time Slices" % constants()->numTimeSlices();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Imaginary Time Step" % constants()->tau();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Initial Number Particles" % params["number_particles"].as<int>();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Initial Density" 
		% (1.0*params["number_particles"].as<int>()/boxPtr->volume);
	communicate()->logFile() << 
		format("%-24s\t:\t%s\n") % "Container Type" % boxPtr->name;
	communicate()->logFile() << format("%-24s\t:\t") % "Container Dimensions";
	for (int i = 0; i < NDIM; i++) {
		communicate()->logFile() << format("%7.5f") % boxPtr->side[i];
		if (i < (NDIM-1))
			communicate()->logFile() << " x ";
		else
			communicate()->logFile() << endl;
	}
	communicate()->logFile() << format("%-24s\t:\t%7.5f\n") % "Container Volume" 
		% boxPtr->volume;
	communicate()->logFile() << format("%-24s\t:\t") % "Lookup Table";
	for (int i = 0; i < NDIM; i++) {
		communicate()->logFile() << format("%d") % nnGrid[i];
		if (i < (NDIM-1))
			communicate()->logFile() << " x ";
		else
			communicate()->logFile() << endl;
	}
	communicate()->logFile() << format("%-24s\t:\t%7.5f\n") % "Initial Worm Constant" % 
		params["worm_constant"].as<double>();
	communicate()->logFile() << format("%-24s\t:\t%7.5f\n") % "Worm Constant" % constants()->C0();
	communicate()->logFile() << format("%-24s\t:\t%7.5f\n") % "Inital CoM Shift" % params["Delta"].as<double>();
	communicate()->logFile() << format("%-24s\t:\t%7.5g\n") % "CoM Shift" % constants()->Delta();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Bisection Parameter" % constants()->b();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Update Slices (Mbar)" % constants()->Mbar();
	communicate()->logFile() << 
		format("%-24s\t:\t%7.5f\n") % "Potential Cutoff Length" 
		% params["potential_cutoff"].as<double>();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Number EQ Steps" % params["number_eq_steps"].as<uint32>();
	communicate()->logFile() << 
		format("%-24s\t:\t%d\n") % "Number Bins Stored" % params["number_bins_stored"].as<int>();
	communicate()->logFile() << format("%-24s\t:\t%d\n") % "Random Number Seed" % _seed;
	communicate()->logFile() << endl;
	communicate()->logFile() << "---------- End Simulation Parameters ------------" << endl;

}
