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
#include "wavefunction.h"
#include "action.h"
#include "move.h"
#include "estimator.h"

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

	/* Define the allowed interaction potential names */
	interactionPotentialName.push_back("aziz");
	interactionPotentialName.push_back("delta");
	interactionPotentialName.push_back("lorentzian");
	interactionPotentialName.push_back("sutherland");
	interactionPotentialName.push_back("hard_sphere");
	interactionPotentialName.push_back("hard_rod");
	interactionPotentialName.push_back("free");
    interactionPotentialName.push_back("delta1D");
    interactionPotentialName.push_back("harmonic");

    interactionNames = getOptionList(interactionPotentialName);

	/* Define the allowed external  potential names */
	externalPotentialName.push_back("free");
	externalPotentialName.push_back("harmonic");
	externalPotentialName.push_back("osc_tube");
	externalPotentialName.push_back("lj_tube");
	externalPotentialName.push_back("hard_tube");
	externalPotentialName.push_back("hg_tube");
	externalPotentialName.push_back("fixed_aziz");
    externalPotentialName.push_back("gasp_prim");

    externalNames = getOptionList(externalPotentialName);

	/* Define the allowed action names */
	actionName.push_back("primitive");
	actionName.push_back("li_broughton");
	actionName.push_back("gsf");
	actionName.push_back("pair_product");

    actionNames = getOptionList(actionName);

    /* Define the allowed trial wave function names */
	waveFunctionName.push_back("constant");
	waveFunctionName.push_back("sech");
	waveFunctionName.push_back("jastrow");

    waveFunctionNames = getOptionList(waveFunctionName);
}

/**************************************************************************//**
 * Make a list of possible option strings
 *
 * @param option the stl vector of options
 * @return a comma separated list of options
******************************************************************************/
string Setup::getOptionList(const vector<string> option) {

    ostringstream optionList;
    std::copy(option.begin(), option.end() - 1, std::ostream_iterator<string>(optionList, ", "));
    optionList << option.back();
    return optionList.str();
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
		("version", "output svn version")
		("output_config,o", po::value<int>()->default_value(0), "number of output configurations")
		("process,p", po::value<uint32>()->default_value(0), "process or cpu number")
		("restart,R", po::value<uint32>()->default_value(0), 
		 "restart running simulation with PIMCID")
        ("wall_clock,W", po::value<double>()->default_value(0),
         "set wall clock limit in hours")
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
		("delta_radius", po::value<double>(), "differential radius for hourglass potential [angstroms]")
		("estimator_radius,w", po::value<double>()->default_value(2.0),
		 "maximum radius for cylinder estimators") 
        ("empty_width_y,y", po::value<double>(), "how much space (in y-) around barrier [Gasparini]")
        ("empty_width_z,z", po::value<double>(), "how much space (in z-) around barrier [Gasparini]")
		;

	potentialOptions.add_options()
		("interaction_potential,I", po::value<string>()->default_value("aziz"), 
		 str(format("interaction potential type {%s}") % interactionNames).c_str())
		("external_potential,X", po::value<string>()->default_value("free"),
		 str(format("external potential type {%s}") % externalNames).c_str())
		("scattering_length,a", po::value<double>()->default_value(1.0),
		 "scattering length [angstroms]")
		("delta_width", po::value<double>()->default_value(1.0E-3),
		 "delta function potential width")
		("delta_strength,c", po::value<double>()->default_value(10.0),
		 "delta function potential integrated strength")
		("interaction_strength,g", po::value<double>()->default_value(1.0),
		 "interaction parameter")
        ("omega", po::value<double>()->default_value(1.0),
         "harmonic interaction potential frequency")
		("fixed,f", po::value<string>()->default_value(""),
		 "input file name for fixed atomic positions.")
		("potential_cutoff,l", po::value<double>(), "interaction potential cutoff length [angstroms]")
		;

	physicalOptions.add_options()
		("canonical", "perform a canonical simulation")
		("pigs", "perform a path integral ground state (PIGS) simulation at T = 0")
        ("window", po::value<int>()->default_value(-1),"set particle number window")
        ("gaussian_ensemble_SD", po::value<double>()->default_value(-1.0),"set gaussian ensemble weight")
		("mass,m", po::value<double>()->default_value(4.0030),"particle mass [amu]")
		("density,n", po::value<double>(), str(format("initial density [angstroms^(-%d)]") % NDIM).c_str())
		("number_particles,N", po::value<int>(), "number of particles")
		("temperature,T", po::value<double>(), "temperature [kelvin]")
		("imaginary_time_length", po::value<double>(), "total path length in imaginary time [kelvin^(-1)]")
        ("wavefunction", po::value<string>()->default_value("constant"),
         str(format("trial wave function type {%s}") % waveFunctionNames).c_str())
        ("end_factor", po::value<double>()->default_value(1.0),"end bead potential action multiplicatave factor")
		("chemical_potential,u", po::value<double>()->default_value(0.0), "chemical potential [kelvin]")
        ("number_broken", po::value<int>()->default_value(0), "number of broken world-lines")
        ("spatial_subregion", po::value<double>()->default_value(BIG), "define a spatial subregion")
        ("number_paths", po::value<int>()->default_value(1), "number of paths")
		;

	algorithmicOptions.add_options()
		("relax", "perform a worm constant relaxation")
		("relaxmu", "perform a chemical potential relaxation to target a fixed density")
		("number_time_slices,P", po::value<int>(), "number of time slices")
		("imaginary_time_step,t", po::value<double>(), "imaginary time step [kelvin^(-1)]")
		("worm_constant,C", po::value<double>()->default_value(1.0), "worm acceptance constant")
		("com_delta,D", po::value<double>(),"center of mass update radius[angstroms]")
		("displace_delta,d", po::value<double>(),"displace update radius [angstroms]")
		("Mbar,M", po::value<int>(), "worm update length, Mbar")
		("number_eq_steps,E", po::value<uint32>()->default_value(1), 
		 "number of equilibration steps")
		("number_bins_stored,S", po::value<int>()->default_value(1), 
		 "number of estimator bins stored")
        ("bin_size", po::value<int>()->default_value(100),
         "number of updates per bin")
		("action", po::value<string>()->default_value("gsf"),
		 str(format("action type {%s}") % actionNames).c_str())
		("full_updates", "perform full staging updates")
		("staging", "perform staging instead of bisection for the diagonal update")
		("max_wind", po::value<int>()->default_value(1), "the maximum winding number to sample")
        ("virial_window,V", po::value<int>()->default_value(5), "centroid virial energy estimator window")
        ("time_resolved", "compute time resolved estimators")
		("no_save_state", "Only save a state file at the end of a simulation")
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

	/* Output the svn version that the code was compiled with */
	if (params.count("version")) {
		cout << endl << format("Code was compiled with repo version %s.") % SVN_VERSION
			<< endl << endl;
		return true;
	}

    /* Have we defined a temperature for a PIMC simulation?*/
    if (!params.count("temperature") && !params.count("pigs")) {
		cerr << endl << "ERROR: No temperature defined!" << endl << endl;
		cerr << "Action: specify temperature (T)" << endl;
    }

    /* Have we mistakingly defined one for a PIGS simulation? */
    if (params.count("temperature") && params.count("pigs")) {
		cerr << endl << "ERROR: Did you mean to define a non-zero temperature for PIGS?" << endl << endl;
		cerr << "Action: remove temperature (T)" << endl;
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
		cerr << endl << "ERROR: Cannot create the simulation cell!" << endl << endl;
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
		cerr << endl << "ERROR: Invalid simulation cell type." << endl << endl;
		cerr << "Action: change cell_type (b) to one of:" << endl
			<< "\t[prism,cylinder]" << endl;
		return true;
	}

	/* Can we create the worldlines? */
	if (!((params.count("imaginary_time_length") && params.count("number_time_slices")) || 
          (params.count("imaginary_time_length") && params.count("imaginary_time_step")) || 
	      (params.count("number_time_slices") && params.count("imaginary_time_step")) ||
          (params.count("temperature") && (params.count("number_time_slices") || (params.count("imaginary_time_step")))) ) ) {
		cerr << endl << "ERROR: Cannot create imaginary time paths!" << endl << endl;
		cerr << "Action: define imaginary_time_length AND imaginary_time_step (t) OR" << endl;
		cerr << "        define imaginary_time_length AND number_time_steps (P) OR" << endl; 
		cerr << "        define number_time_slices (P) AND imaginary_time_step (t)" << endl << endl;
		cerr << "        define temperature (T) AND (number_time_slices (P) OR imaginary_time_step (t))" << endl << endl;
		cerr << algorithmicOptions << endl;
		return true;
	}

	/* Make sure we have selected a valid interaction potential */
    if (std::find(interactionPotentialName.begin(), interactionPotentialName.end(), 
                params["interaction_potential"].as<string>()) == interactionPotentialName.end()) {
		cerr << endl << "ERROR: Invalid interaction potential!" << endl << endl;
		cerr << "Action: set interaction_potential (I) to one of:" << endl
			 << "\t[" << interactionNames << "]" <<  endl;
		return true;
	}

	/* Make sure we have selected a valid external potential */
    if (std::find(externalPotentialName.begin(), externalPotentialName.end(), 
                params["external_potential"].as<string>()) == externalPotentialName.end()) {
		cerr << endl << "ERROR: Invalid external potential!" << endl << endl;
		cerr << "Action: set external_potential (X) must to one of:" << endl
			 << "\t[" << externalNames << "]" << endl;
		return true;
	}

	/* Make sure we have selected a valid action type */
    if (std::find(actionName.begin(), actionName.end(), 
                params["action"].as<string>()) == actionName.end()) {
		cerr << endl << "ERROR: Invalid action!" << endl << endl;
		cerr << "Action: set action to one of:" << endl
			 << "\t[" << actionNames << "]" <<  endl;
		return true;
	}
    
    /* Make sure we have selected a valid trial wave fucntion */
    if (std::find(waveFunctionName.begin(), waveFunctionName.end(), 
                params["wavefunction"].as<string>()) == waveFunctionName.end()) {
		cerr << endl << "ERROR: Invalid trial wave function!" << endl << endl;
		cerr << "Action: set wave function to one of :" << endl
        << "\t[" << waveFunctionNames << "]" << endl;
		return true;
	}

    /* Make sure we use the pair product approximation for discontinuous
     * potentials */
    if (((params["interaction_potential"].as<string>() == "hard_sphere") ||
        (params["interaction_potential"].as<string>() == "hard_rod") ||
        (params["interaction_potential"].as<string>() == "delta1D") ) &&
	    (params["action"].as<string>() != "pair_product") ) {
        cout << endl;
		cerr << "ERROR: Need to use the pair product approximation with discontinuous potentials!"; 
        cout << endl << endl;
		cerr << "Action: change the action to pair_product." << endl;
		return 1;
    }

	/* We can only use the cylinder potentials for a 3D system */
	if ((params["external_potential"].as<string>().find("tube") != string::npos) && (NDIM != 3)) {
		cerr << endl << "ERROR: Can only use tube potentials for a 3D system!" << endl << endl;
		cerr << "Action: change the potential or recompile with ndim=3." << endl;
		return 1;
	}

	/* Need to specify a radius for the tube potentials */
	if ( (params["external_potential"].as<string>().find("tube") != string::npos) && 
            (!params.count("radius")) ) {
		cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
		cerr << "Action: specfiy a radius (r) for the tube potentials." << endl;
		return 1;
	}

    /* Need to specify a y- barrier width scale factor for Gasparini potential */
	if ( (params["external_potential"].as<string>().find("gasp_prim") != string::npos) &&
             (!params.count("empty_width_y")) ) {
		cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
		cerr << "Action: specify a y- scale factor (y) for the Gasparini potential." << endl;
		return 1;
	}

    /* Need to specify a z- barrier width scale factor for Gasparini potential */
	if ( (params["external_potential"].as<string>().find("gasp_prim") != string::npos) && 
            (!params.count("empty_width_z")) ) {
		cerr << endl << "ERROR: Incomplete specification for external potential!" << endl << endl;
		cerr << "Action: specify a z- scale factor (z) for the Gasparini potential." << endl;
		return 1;
	}

	/* We can only use the hard sphere potential in a 3D system */
	if ((params["interaction_potential"].as<string>().find("hard_sphere") != string::npos) && (NDIM != 3)) {
		cerr << endl << "ERROR: Can only use hard sphere potentials for a 3D system!" << endl << endl;
		cerr << "Action: change the potential or recompile with ndim=3." << endl;
		return 1;
	}

	/* We can only use the hard rod potential in a 1D system */
	if ((params["interaction_potential"].as<string>().find("hard_rod") != string::npos) && (NDIM != 1)) {
		cerr << endl << "ERROR: Can only use hard rod potentials for a 1D system!" << endl << endl;
		cerr << "Action: change the potential or recompile with ndim=1." << endl;
		return 1;
	}
    
    /* We can only use the delta1D potential in a 1D system */
	if ((params["interaction_potential"].as<string>().find("delta1D") != string::npos) && (NDIM != 1)) {
		cerr << endl << "ERROR: Can only use delta1D potentials for a 1D system!" << endl << endl;
		cerr << "Action: change the potential or recompile with ndim=1." << endl;
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

	int numTimeSlices,numDeltaTau;
    double tau,imaginaryTimeLength;
    bool pathBreak = ((params["number_broken"].as<int>() != 0)||
                                    (params["spatial_subregion"].as<double>() < BIG/2.0));

    /* !!NB!! We separate the worldline building into two sections, one for PIGS and
     * one for PIMC.  This can probably be cleaned up in the future. */

    if (!params.count("pigs")) {

        /* Set the imaginary time length*/
        if (!params.count("imaginary_time_length"))
            insertOption("imaginary_time_length",1.0/params["temperature"].as<double>());

        /* We determine if we have fixed P or tau.  We require that the number of time 
         * slices is always even. */
        if (!params.count("number_time_slices")) {
            tau = params["imaginary_time_step"].as<double>();
            numTimeSlices = static_cast<int>(1.0/(params["temperature"].as<double>() * tau) + EPS);
            if ((numTimeSlices % 2) != 0)
                numTimeSlices--;
            insertOption("number_time_slices",numTimeSlices);
        }
        else {
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) != 0)
                numTimeSlices--;
            tau = 1.0/(params["temperature"].as<double>() * numTimeSlices);
            setOption("number_time_slices",numTimeSlices);
            insertOption("imaginary_time_step",tau);
        }
    }
    else {
        /* Fixing the total imaginary time and the number of time slices */
        if ( params.count("imaginary_time_length") && params.count("number_time_slices") ) {

            /* We make sure we have an odd number of time slices */
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;
            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;

            tau = params["imaginary_time_length"].as<double>() / (numTimeSlices-1);
            setOption("number_time_slices",numTimeSlices);
            insertOption("imaginary_time_step",tau);
        }
        /* Fixing the total imaginary time and the time step.  We need to make sure
         * we can get an odd number of slices */
        else if ( params.count("imaginary_time_length") && params.count("imaginary_time_step") ) {
            tau = params["imaginary_time_step"].as<double>();
            numTimeSlices = static_cast<int>((params["imaginary_time_length"].as<double>() / tau) + EPS)+1;
            
            /* We make sure we have an odd number of time slices */
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;
            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;
            
            imaginaryTimeLength = (numTimeSlices-1)*tau;
            setOption("imaginary_time_length",imaginaryTimeLength);
            
            insertOption("number_time_slices",numTimeSlices);
        }
        /* Fixing the number of time steps and the size of the time step.  */
        else {
            /* We make sure we have an odd number of time slices */
            numTimeSlices = params["number_time_slices"].as<int>();
            if ((numTimeSlices % 2) == 0)
                numTimeSlices++;
            /* Make sure the center of the path is an odd (even) slice for a
             closed (open) path calucation */
            numDeltaTau = (numTimeSlices-1)/2;
            if ( ((pathBreak)&&(numDeltaTau%2==0)) || ((!pathBreak)&&(numDeltaTau%2==1))    )
                numTimeSlices += 2;

            /* Set the updated number of time slices and imaginary time length */
            setOption("number_time_slices",numTimeSlices);
            imaginaryTimeLength = (numTimeSlices-1)*params["imaginary_time_step"].as<double>();
            insertOption("imaginary_time_length",imaginaryTimeLength);
        }

        /* We set the effective temperature to 1.0/imagTimeLength */
        insertOption("temperature",1.0/params["imaginary_time_length"].as<double>());
    }

	/* If we haven't fixed Mbar, do so now. We use the value defined in
	 * PRE 74, 036701 (2006). We make sure it is always >= 8 */
	int Mbar = 0;
	if (!params.count("Mbar")) {
		/* Compute the average inter-particle separation */
		double ell = pow((1.0*params["number_particles"].as<int>()/volume),-1.0/(1.0*NDIM));
		double lambda = 24.24 / params["mass"].as<double>();
		Mbar = int(ell*ell/(16*lambda*params["imaginary_time_step"].as<double>()));
		if (Mbar < 2)
			Mbar = 2;
        else if (Mbar > numTimeSlices)
            Mbar = numTimeSlices/2;
		insertOption("Mbar",Mbar);
	}
	else 
		Mbar = params["Mbar"].as<int>();

	/* Now we make sure it is even and not too large*/
	if (Mbar%2)
		Mbar++;

	setOption("Mbar",Mbar);

	if (Mbar > numTimeSlices) {
        cerr << Mbar << " " << numTimeSlices << endl;
		cerr << endl << "ERROR: Update length > number time slices!" << endl << endl;
        cerr << "Action: Increase number_time_slices (P) OR" <<  endl;
        cerr << "        Increase Mbar (M) OR"  << endl;
        cerr << "        Decrease imaginary_time_step (t) OR" << endl; 
        cerr << "        Increase imaginary_time_length" << endl; 
        return true;
	}

	return false;
}

/**************************************************************************//**
 * Setup the simulation constants.
 *
 * Fix all simulation constants.
******************************************************************************/
void Setup::setConstants() {

    /* At present, we need to make sure that if a pair_product action has been
     * selected, that we turn off the cuttoff by making it the size of the box */
    if (params["action"].as<string>() == "pair_product"){
        if (params.count("potential_cutoff"))
            setOption("potential_cutoff",side[NDIM-1]);
        else
            insertOption("potential_cutoff",side[NDIM-1]);
    }

    /* If we haven't set a potential_cutoff, set it to the size of the box */
    if (!params.count("potential_cutoff"))
        insertOption("potential_cutoff",side[NDIM-1]);

	/* Set the required constants */
	constants()->initConstants(
            params.count("pigs"),
			params.count("canonical"),
			params["temperature"].as<double>(),
            params["imaginary_time_length"].as<double>(),
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
            params["wall_clock"].as<double>(),
			params["number_eq_steps"].as<uint32>(),
			params["interaction_potential"].as<string>(),
			params["external_potential"].as<string>(),
            params["wavefunction"].as<string>(),
            params["action"].as<string>(),
            params["window"].as<int>(),
            params["gaussian_ensemble_SD"].as<double>(),
            params["max_wind"].as<int>(),
            params["virial_window"].as<int>(),
            params["number_broken"].as<int>(),
            params["spatial_subregion"].as<double>(),
            params["end_factor"].as<double>(),
            params["number_paths"].as<int>(),
            (!params.count("no_save_state"))
            );

    /* If we have specified either the center of mass or displace shift on the
     * command line, we update their values. */
	if (params.count("com_delta"))
		constants()->setCoMDelta(params["com_delta"].as<double>());
	else {
		insertOption("com_delta",constants()->comDelta());
	}
	if (params.count("displace_delta"))
		constants()->setDisplaceDelta(params["displace_delta"].as<double>());
	else {
		insertOption("displace_delta",constants()->displaceDelta());
	}
}

/**************************************************************************//**
 * Setup the communicator.
 *
 * Initialize the communicator, we need to know if we are outputting any config
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
	else if (constants()->intPotentialType() == "sutherland")
		interactionPotentialPtr = new SutherlandPotential(params["interaction_strength"].as<double>());
	else if (constants()->intPotentialType() == "hard_sphere")
        interactionPotentialPtr = new HardSpherePotential(params["scattering_length"].as<double>());
	else if (constants()->intPotentialType() == "hard_rod")
        interactionPotentialPtr = new HardRodPotential(params["scattering_length"].as<double>());
    else if (constants()->intPotentialType() == "delta1D")
        interactionPotentialPtr = new Delta1DPotential(params["delta_strength"].as<double>());
	else if (constants()->intPotentialType() == "lorentzian")
		interactionPotentialPtr = new LorentzianPotential(params["delta_width"].as<double>(),
				params["delta_strength"].as<double>());
	else if (constants()->intPotentialType() == "aziz")
		interactionPotentialPtr = new AzizPotential(side);
    else if (constants()->intPotentialType() == "harmonic")
		interactionPotentialPtr = new HarmonicPotential(params["omega"].as<double>());

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
	else if (constants()->extPotentialType() == "hg_tube")
		externalPotentialPtr = new LJHourGlassPotential(boxPtr->side[NDIM-1],
                params["radius"].as<double>(), params["delta_radius"].as<double>());
	else if (constants()->extPotentialType() == "hard_tube") 
		externalPotentialPtr = new HardCylinderPotential(params["radius"].as<double>());
	else if (constants()->extPotentialType() == "single_well")
		externalPotentialPtr = new SingleWellPotential();
	else if (constants()->extPotentialType() == "fixed_aziz") 
		externalPotentialPtr = new FixedAzizPotential(boxPtr);
    else if (constants()->extPotentialType() == "gasp_prim")
        externalPotentialPtr = new Gasparini_1_Potential(params["empty_width_z"].as<double>(),
                params["empty_width_y"].as<double>(),boxPtr);

	return externalPotentialPtr;
}


/*************************************************************************//**
* Setup the trial wave function.
*
* Based on the user's choice we create a new trial wave function  pointer
* which is returned to the main program.
******************************************************************************/
WaveFunctionBase * Setup::waveFunction(const Path &path) {
    
	WaveFunctionBase *waveFunctionPtr = NULL;
    
	if (constants()->waveFunctionType() == "constant")
		waveFunctionPtr = new WaveFunctionBase(path);
	else if (constants()->waveFunctionType() == "sech")
		waveFunctionPtr = new SechWaveFunction(path);
	else if (constants()->waveFunctionType() == "jastrow")
		waveFunctionPtr = new JastrowWaveFunction(path);
    
	return waveFunctionPtr;
}

/*************************************************************************//**
* Setup the action.
*
* Based on the user's choices we create a new action pointer which is returned
* to the main program.
******************************************************************************/
ActionBase * Setup::action(const Path &path, LookupTable &lookup, 
        PotentialBase *externalPotentialPtr, 
        PotentialBase * interactionPotentialPtr, 
        WaveFunctionBase * waveFunctionPtr) {

    ActionBase *actionPtr = NULL;

    /* Non Local Actions */
    if (constants()->actionType() == "pair_product") {
        actionPtr = new NonLocalAction(path,lookup,externalPotentialPtr,
                interactionPotentialPtr,waveFunctionPtr,false,constants()->actionType());	
    }
    /* Local Actions */
    else {

        /* The factors needed for local actions. */
        TinyVector <double,2> VFactor;	    
        TinyVector <double,2> gradVFactor;	

        VFactor = 1.0;
        gradVFactor = 0.0;

        if (constants()->actionType() == "gsf") {

            /* There are different pre-factors for even/odd slices, we use the 
             * Prokof'ev version here. I have taken alpha = 0 here. */
            VFactor[0] = 2.0/3.0;
            VFactor[1] = 4.0/3.0;

            double alpha = 0.0;
            gradVFactor[0] = 2.0*alpha/9.0;
            gradVFactor[1] = 2.0*(1.0-alpha)/9.0;
        }
        else if (constants()->actionType() == "li_broughton")
            gradVFactor = 1.0 / 12.0;

        /* Do we want to use full staging moves? */
        bool local = true;
        if (params.count("full_updates"))
            local = false;

        actionPtr = new LocalAction(path,lookup,externalPotentialPtr,
                interactionPotentialPtr,waveFunctionPtr,VFactor,gradVFactor,local,
                                    constants()->actionType(),constants()->endFactor());
        
        }

    return actionPtr;
}

/*************************************************************************//**
 * Define the Monte Carlo updates that will be performed
 *
 * @param path A reference to the paths
 * @param actionPtr The action in use
 * @param random The random number generator 
 * @return a list of Monte Carlo updates
******************************************************************************/
auto_ptr< boost::ptr_vector<MoveBase> > Setup::moves(Path &path, 
        ActionBase *actionPtr, MTRand &random) {

    boost::ptr_vector<MoveBase> move;

    /* All simulations include the Center of Mass move */
    move.push_back(new CenterOfMassMove(path,actionPtr,random));

    /* PIGS simulations use staging and displace moves */
    if (params.count("pigs")) {
        move.push_back(new StagingMove(path,actionPtr,random));
        move.push_back(new EndStagingMove(path,actionPtr,random));

        if (constants()->numBroken() > 0) {
            move.push_back(new SwapBreakMove(path,actionPtr,random));
            constants()->setAttemptProb("diagonal",0.5);
	    constants()->setAttemptProb("swap break",0.1);
	    
        } else if (constants()->spatialSubregionOn() ){
            move.push_back(new MidStagingMove(path,actionPtr,random));
            constants()->setAttemptProb("diagonal",0.5);
	    constants()->setAttemptProb("mid-staging",0.1);
        }

        move.push_back(new DisplaceMove(path,actionPtr,random));
    }
    else {
        /* We determine which type of diagonal path update we will use */
        if ( (params.count("full_updates")) || params.count("staging") || 
                (params["action"].as<string>() == "pair_product")  ||
                (params["max_wind"].as<int>() > 1)  ) {
            move.push_back(new StagingMove(path,actionPtr,random));
        }
        else 
            move.push_back(new BisectionMove(path,actionPtr,random));

        /* Include all other moves */
        move.push_back(new OpenMove(path,actionPtr,random));
        move.push_back(new CloseMove(path,actionPtr,random));
        move.push_back(new InsertMove(path,actionPtr,random));
        move.push_back(new RemoveMove(path,actionPtr,random));
        move.push_back(new AdvanceHeadMove(path,actionPtr,random));
        move.push_back(new RecedeHeadMove(path,actionPtr,random));
        move.push_back(new AdvanceTailMove(path,actionPtr,random));
        move.push_back(new RecedeTailMove(path,actionPtr,random));
        move.push_back(new SwapHeadMove(path,actionPtr,random));
        move.push_back(new SwapTailMove(path,actionPtr,random));
    }

    return move.release();
}

/*************************************************************************//**
 * Create a list of estimators to be measured
 *
 * @param path A reference to the paths
 * @param actionPtr The action in use
 * @param random The random number generator 
 * @return a list of estimators
******************************************************************************/
auto_ptr< boost::ptr_vector<EstimatorBase> > Setup::estimators(Path &path, 
        ActionBase *actionPtr, MTRand &random) {

    boost::ptr_vector<EstimatorBase> estimator;
    
    if (params.count("pigs")) {
        estimator.push_back(new PigsEnergyEstimator(path,actionPtr));
        //estimator.push_back(new ParticleResolvedPositionEstimator(path));
        //estimator.push_back(new ParticleCorrelationEstimator(path));
        if (params.count("time_resolved")){
            estimator.push_back(new PotentialEnergyEstimator(path,actionPtr));
            estimator.push_back(new KineticEnergyEstimator(path,actionPtr));
        }
        //estimator.push_back(new SubregionOccupationEstimator(path,actionPtr));
        //estimator.push_back(new TotalEnergyEstimator(path,actionPtr));
        //estimator.push_back(new ThermoPotentialEnergyEstimator(path,actionPtr));
        //estimator.push_back(new VelocityEstimator(path));
        //estimator.push_back(new PIGSOneBodyDensityMatrixEstimator(path,actionPtr,random));
    }
    /* The finite temperature measurments */
    else {

        /* !!NB!! Scalar estimators, the order here is important! */
        estimator.push_back(new VirialEnergyEstimator(path,actionPtr));
        // estimator.push_back(new EnergyEstimator(path,actionPtr));
        estimator.push_back(new NumberParticlesEstimator(path));
        estimator.push_back(new DiagonalFractionEstimator(path));
        estimator.push_back(new SuperfluidFractionEstimator(path));
        estimator.push_back(new WormPropertiesEstimator(path));

        /* Vector estimators */
        //estimator.push_back(new PermutationCycleEstimator(path));
        
        /* We only measure the number distribution function if we are grand
         * canonical */
        // if (!constants()->canonical())
        //     estimator.push_back(new NumberDistributionEstimator(path));

        // estimator.push_back(new PairCorrelationEstimator(path,actionPtr));

        //estimator.push_back(new OneBodyDensityMatrixEstimator(path,actionPtr,random));

        /* Time Consuming local estimators */
    //    /* z-averaged estimators */
    //    estimator.push_back(new ParticlePositionEstimator(path));
    //    estimator.push_back(new PlaneParticlePositionEstimator(path));
    //    estimator.push_back(new LocalSuperfluidDensityEstimator(path));
    //    estimator.push_back(new PlaneWindingSuperfluidDensityEstimator(path));
    //    estimator.push_back(new PlaneAreaSuperfluidDensityEstimator(path));
    //    estimator.push_back(new LocalPermutationEstimator(path));
    
        /* Excluded volume estimators */
        if (constants()->extPotentialType().find("gasp_prim") != string::npos){
            estimator.push_back(new BipartitionDensityEstimator(path,actionPtr));
        }

        /* Cylinder estimators */
        if (constants()->extPotentialType().find("tube") != string::npos) {
            double maxR = params["estimator_radius"].as<double>();
            estimator.push_back(new CylinderEnergyEstimator(path,actionPtr,maxR));
            estimator.push_back(new CylinderNumberParticlesEstimator(path,maxR));
            estimator.push_back(new CylinderSuperfluidFractionEstimator(path,maxR));
            estimator.push_back(new CylinderRadialPotentialEstimator(path,actionPtr,random,maxR));
            estimator.push_back(new CylinderNumberDistributionEstimator(path,maxR));
            estimator.push_back(new CylinderPairCorrelationEstimator(path,actionPtr,maxR));
            estimator.push_back(new CylinderOneBodyDensityMatrixEstimator(path,actionPtr,random,maxR));

            /* Radially averaged estimators */
            estimator.push_back(new RadialWindingSuperfluidDensityEstimator(path));
            estimator.push_back(new RadialAreaSuperfluidDensityEstimator(path));
            estimator.push_back(new RadialDensityEstimator(path));
        }
    } // pimc

    return estimator.release();
}


/*************************************************************************//**
* Create a list of double path estimators to be measured
*
* @param path A reference to the paths
* @param actionPtr The action in use
* @return a list of double path estimators
******************************************************************************/
auto_ptr< boost::ptr_vector<EstimatorBase> > Setup::multiPathEstimators(
        vector<Path *> &pathPtrVec,vector<ActionBase *> &actionPtrVec) {
    
    boost::ptr_vector<EstimatorBase> multiEstimators;
    
    multiEstimators.push_back(new SwapEstimator(*pathPtrVec[0],*pathPtrVec[1],
                                                  actionPtrVec[0],actionPtrVec[1]));
    //doubledEstimators.push_back(new EntPartEstimator(path,path2,actionPtr,actionPtr2));

    return multiEstimators.release();
}


/*************************************************************************//**
 * Compare a char array with an option name
 *
 * @param argv The commmand line string
 * @param target The target option name
******************************************************************************/
bool Setup::checkOption(const string option, const string target) {

    /* check short options first  */
    if (option[0] == '-' && option[1] != '-')
        return (option == target);

    /* Now test for a long option, making sure to ignore any short
     * option targets */
    else if (target[0] != '-')
        return (option.find(target) != string::npos);

    /* otherwise return false */
    return false;
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

	communicate()->file("log")->stream() << endl << "# ";

	/* Construct the command that would be required to restart the simulation */
	bool outputC0 = false;
	bool outputD = false;
	bool outputd = false;

	for (int n = 0; n < argc; n++) {

        /* convert the argument to a string */
        string arg(argv[n]);

		if (checkOption(arg,"-s") || checkOption(arg,"start_with_state") ) {
            n++;
        }
		else if ( checkOption(arg,"-C") || checkOption(arg,"worm_constant") ) {
			communicate()->file("log")->stream() << format("-C %10.4e ") % constants()->C0();
			n++;
			outputC0 = true;
		}
		else if (checkOption(arg,"-p") || checkOption(arg,"process")) {
			communicate()->file("log")->stream() << format("-p %03d ") % params["process"].as<uint32>();
			n++;
		}
		else if (checkOption(arg,"-D") || checkOption(arg,"com_delta")) {
            communicate()->file("log")->stream() << format("-D %10.4e ") % constants()->comDelta();
			outputD = true;
			n++;
		}
		else if (checkOption(arg,"-d") || checkOption(arg,"displace_delta")) {
            communicate()->file("log")->stream() << format("-d %10.4e ") % constants()->displaceDelta();
			outputd = true;
			n++;
		}
		else {
			communicate()->file("log")->stream() << arg << " ";
        }

	}

    /* Output the restart flag */
    communicate()->file("log")->stream() << format("-R %09d ") % constants()->id();

	/* If we haven't specified the worm constant, output it now */
	if (!outputC0)
		communicate()->file("log")->stream() << format("-C %21.15e ") % constants()->C0();

	/* If we haven't specified the center of mass Delta, output it now */
	if (!outputD)
        communicate()->file("log")->stream() << format("-D %21.15e ") % constants()->comDelta();
    
	/* If we haven't specified the displace delta, output it now */
    if (constants()->pigs() && !outputd)
        communicate()->file("log")->stream() << format("-d %21.15e ") % constants()->displaceDelta();

	communicate()->file("log")->stream() << endl << endl;
	communicate()->file("log")->stream() << "---------- Begin Simulation Parameters ----------" << endl;
	communicate()->file("log")->stream() << endl;
	if (constants()->canonical())
		communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Ensenble" % "canonical";
	else
		communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Ensenble" % "grand canonical";
	if (constants()->pigs())
		communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Simulation Type" % "PIGS";
	else
		communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Simulation Type" % "PIMC";
	communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Action Type" % params["action"].as<string>();
    communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Number of paths" % params["number_paths"].as<int>();
	communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") % "Interaction Potential" % 
        params["interaction_potential"].as<string>();

    /* Ouptut a possible delta function width and strength */
	if ( (params["interaction_potential"].as<string>().find("delta") != string::npos) ||
        (params["interaction_potential"].as<string>().find("lorentzian") != string::npos) ) {
		communicate()->file("log")->stream() << format("%-24s\t:\t%8.3e\n") % "Delta Width"
            % params["delta_width"].as<double>();
		communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Delta Strength"
            % params["delta_strength"].as<double>();
	}

    /* Ouptut a possible sutherland model interaction strength*/
	if (params["interaction_potential"].as<string>().find("sutherland") != string::npos) {
		communicate()->file("log")->stream() << format("%-24s\t:\t%8.3e\n") % "Interaction Strength"
            % params["interaction_strength"].as<double>();
	}

    // if ( (params["interaction_potential"].as<string>().find("delta1D") != string::npos) ) {
    //     communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Delta Strength"
    //         % params["delta_strength"].as<double>();
	// }

    /* Output harmonic interaction frequecy */
    if ( (params["interaction_potential"].as<string>().find("harmonic") != string::npos) ) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Harmonic Int. Freq."
            % params["omega"].as<double>();
	}

    /* Output a possible scattering length */
	if ( (params["interaction_potential"].as<string>().find("hard_sphere") != string::npos) ||
         (params["interaction_potential"].as<string>().find("hard_rod") != string::npos) ) {
		communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") % "Scattering Length" 
			% params["scattering_length"].as<double>();
    }
	communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") 
		% "External Potential" % params["external_potential"].as<string>();

    /* output a possible delta radius for the hourglass potential. */
	if (params["external_potential"].as<string>().find("hg_tube") != string::npos) 
        communicate()->file("log")->stream() << format("%-24s\t:\t%-7.2f\n") 
            % "Delta Radius" % params["delta_radius"].as<double>();

	communicate()->file("log")->stream() << format("%-24s\t:\t%s\n") 
		% "Wavefunction Type" % params["wavefunction"].as<string>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Temperature" % params["temperature"].as<double>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Chemical Potential" % constants()->mu();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Particle Mass" % params["mass"].as<double>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Number Time Slices" % constants()->numTimeSlices();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Imaginary Time Step" % constants()->tau();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Imaginary Time Length" % constants()->imagTimeLength();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Initial Number Particles" % params["number_particles"].as<int>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Initial Density" 
		% (1.0*params["number_particles"].as<int>()/boxPtr->volume);
    communicate()->file("log")->stream() <<
        format("%-24s\t:\t%d\n") % "Num. Broken World-lines" % params["number_broken"].as<int>();
    if ( constants()->spatialSubregionOn()){
        communicate()->file("log")->stream() <<
            format("%-24s\t:\t%d\n") % "Spatial Subregion" % params["spatial_subregion"].as<double>();
    }
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%s\n") % "Container Type" % boxPtr->name;
	communicate()->file("log")->stream() << format("%-24s\t:\t") % "Container Dimensions";
	for (int i = 0; i < NDIM; i++) {
		communicate()->file("log")->stream() << format("%7.5f") % boxPtr->side[i];
		if (i < (NDIM-1))
			communicate()->file("log")->stream() << " x ";
		else
			communicate()->file("log")->stream() << endl;
	}
	communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Container Volume" 
		% boxPtr->volume;
	communicate()->file("log")->stream() << format("%-24s\t:\t") % "Lookup Table";
	for (int i = 0; i < NDIM; i++) {
		communicate()->file("log")->stream() << format("%d") % nnGrid[i];
		if (i < (NDIM-1))
			communicate()->file("log")->stream() << " x ";
		else
			communicate()->file("log")->stream() << endl;
	}
	communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Initial Worm Constant" % 
		params["worm_constant"].as<double>();
	communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Worm Constant" % constants()->C0();
	communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Inital CoM Delta" % params["com_delta"].as<double>();
	communicate()->file("log")->stream() << format("%-24s\t:\t%7.5g\n") % "CoM Delta" % constants()->comDelta();
    if (constants()->pigs()) {
        communicate()->file("log")->stream() << format("%-24s\t:\t%7.5f\n") % "Inital Displace Delta" % params["displace_delta"].as<double>();
        communicate()->file("log")->stream() << format("%-24s\t:\t%7.5g\n") % "Displace Delta" % constants()->displaceDelta();
    }
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Bisection Parameter" % constants()->b();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Update Slices (Mbar)" % constants()->Mbar();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%7.5f\n") % "Potential Cutoff Length" % params["potential_cutoff"].as<double>();
    communicate()->file("log")->stream() <<
        format("%-24s\t:\t%d\n") % "Bin Size" % params["bin_size"].as<int>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Number EQ Steps" % params["number_eq_steps"].as<uint32>();
	communicate()->file("log")->stream() << 
		format("%-24s\t:\t%d\n") % "Number Bins Stored" % params["number_bins_stored"].as<int>();
	communicate()->file("log")->stream() << format("%-24s\t:\t%d\n") % "Random Number Seed" % _seed;
    if (params["virial_window"].as<int>() == 0){
        communicate()->file("log")->stream() << 
            format("%-24s\t:\t%d\n") % "Virial Window" % constants()->numTimeSlices();
    }
    else{
        communicate()->file("log")->stream() << 
            format("%-24s\t:\t%d\n") % "Virial Window" % params["virial_window"].as<int>();
    }
	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << "---------- End Simulation Parameters ------------" << endl;
}
