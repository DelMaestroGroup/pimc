/**
 * @file pdrive.cpp
 * @brief A dimensionally independent worm algorithm path integral monte carlo code driver.
 * @author Adrian Del Maestro
 * @date 10.14.2008
 */

#include "common.h"
#include "constants.h"
#include "container.h"
#include "path.h"
#include "potential.h"
#include "action.h"
#include "pimc.h"
#include "move.h"
#include "lookuptable.h"
#include "communicator.h"
#include "setup.h"

/**
 * Main driver.
 * Read in all program options from the user using boost::program_options and setup the simulation
 * cell, initial conditions and both the interaction and external potential. Either equilibrate or
 * restart a simulation, then start measuring. We output all the simulation parameters to disk as a
 * log file so that it can be restart again assigning it a unique PIMCID.
 * @see boost::program_options -- http://www.boost.org/doc/libs/1_43_0/doc/html/program_options.html
 */
int main (int argc, char *argv[]) {

	uint32 seed = 139853;	// The seed for the random number generator

	Setup setup;

	/* Attempt to parse the command line options */
    try {
		setup.getOptions(argc,argv);
    }
    catch(exception& ex) {
        cerr << "error: " << ex.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

	/* Parse the setup options and possibly exit */
	if (setup.parseOptions())
		return 1;

	/* The global random number generator, we add the process number to the seed (for
	 * use in parallel simulations.*/
	seed = setup.seed(seed);
	MTRand random(seed);

	/* Get the simulation box */
	Container *boxPtr = NULL;
	boxPtr = setup.cell(); 

	/* Create the worldliens */
	if (setup.worldlines())
		return 1;

	/* Setup the simulation constants */
	setup.setConstants();

	/* Setup the simulation communicator */
	setup.communicator();

	/* Create and initialize the Nearest Neighbor Lookup Table */
	LookupTable lookup(boxPtr,constants()->numTimeSlices(), constants()->initialNumParticles());

	/* Create and initialize the potential pointers */
	PotentialBase *interactionPotentialPtr = NULL;
	interactionPotentialPtr = setup.interactionPotential();

	PotentialBase *externalPotentialPtr = NULL;
	externalPotentialPtr = setup.externalPotential(boxPtr);

	/* Get the initial conditions associated with the external potential */
	/* Must use the copy constructor as we return a copy */
	Array<dVec,1> initialPos = 
		externalPotentialPtr->initialConfig(boxPtr,random,constants()->initialNumParticles());

	/* Setup the path data variable */
	Path path(boxPtr,lookup,constants()->numTimeSlices(),initialPos);

	/* The potential object, which needs to know about both the external and 
	 * interaction potential */
	Potential potential(path,lookup,externalPotentialPtr,interactionPotentialPtr);

	/* Setup the action */
	GSFAction action(path,&potential);	
	//ActionBase action(path,&potential);

	/* Setup the pimc object */
	PathIntegralMonteCarlo pimc(path,&action,random,setup.params["estimator_radius"].as<double>(),
			!setup.params["start_with_state"].as<string>().empty());

	/* If this is a fresh run, we equilibrate and output simulation parameters to disk */
	if (!constants()->restart()) {

		/* Equilibrate */
		cout << format("[PIMCID: %09d] - Equilibration Stage.") % constants()->id() << endl;
		for (uint32 n = 0; n < constants()->numEqSteps(); n++) 
			pimc.equilStep_test(n,setup.params.count("relax"));

		/* Output simulation details/parameters */
		setup.outputOptions(argc,argv,seed,boxPtr,lookup.getNumNNGrid());
	}

	cout << format("[PIMCID: %09d] - Measurement Stage.") % constants()->id() << endl;

	/* Sample */
	int oldNumStored = 0;
	int outNum = 0;
	int numOutput = setup.params["output_config"].as<int>();
	uint32 n = 0;
	do {
		pimc.step();
		if (pimc.numStoredBins > oldNumStored) {
			oldNumStored = pimc.numStoredBins;
			cout << format("[PIMCID: %09d] - Bin #%4d stored to disk.") % constants()->id() 
				% oldNumStored << endl;
		}
		n++;

		/* Output configurations to disk */
		if ((numOutput > 0) && ((n % numOutput) == 0)) {
			path.outputConfig(outNum);
			outNum++;
		}

	} while (pimc.numStoredBins < setup.params["number_bins_stored"].as<int>());
	cout << format("[PIMCID: %09d] - Measurement complete.") % constants()->id() << endl;

	/* Output Results */
	pimc.finalOutput();

	/* Free up memory */
	delete interactionPotentialPtr;
	delete externalPotentialPtr;
	delete boxPtr;

	initialPos.free();

	return 1;
}
