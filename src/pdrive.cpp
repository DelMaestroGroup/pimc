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
#include "wavefunction.h"
#include "pimc.h"
#include "lookuptable.h"
#include "communicator.h"
#include "setup.h"
#include "cmc.h"
#include "move.h"

/**
 * Main driver.
 * Read in all program options from the user using boost::program_options and setup the simulation
 * cell, initial conditions and both the interaction and external potential. Either equilibrate or
 * restart a simulation, then start measuring. We output all the simulation parameters to disk as a
 * log file so that it can be restart again assigning it a unique PIMCID.
 * @see http://www.boost.org/doc/libs/release/doc/html/program_options.html
 */
int main (int argc, char *argv[]) {

    /* Get initial time */
    time_t start_time = time(NULL);
    time_t current_time; //current time
    bool wallClockReached = false;

    uint32 seed = 139853;   // The seed for the random number generator

    Setup& setup = Setup::instance();

    /* Attempt to parse the command line options */
    try {
        setup.getOptions(argc,argv);
    }
    catch(std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }

    /* Parse the setup options and possibly exit */
    if (setup.parseOptions())
        return 1;

    /* The global random number generator, we add the process number to the seed (for
     * use in parallel simulations.*/
    seed = setup.seed(seed);
    MTRand random(seed);

    /* Get the simulation box */
    setup.set_cell();
    Container *boxPtr = setup.get_cell();

    /* Create the worldlines */
    if (setup.worldlines())
        return 1;

    /* Setup the simulation constants */
    setup.setConstants();

    /* Setup the simulation communicator */
    setup.communicator();

    /* Get number of paths to use */
    int Npaths = constants()->Npaths();
    
    /* Create and initialize the Nearest Neighbor Lookup Table */
    boost::ptr_vector<LookupTable> lookupPtrVec;
    for(int i=0; i<Npaths; i++){
        lookupPtrVec.push_back(
                new LookupTable(boxPtr,constants()->numTimeSlices(),
                    constants()->initialNumParticles()));
    }
    
    /* Create and initialize the potential pointers */
    PotentialBase *interactionPotentialPtr = setup.interactionPotential();
    PotentialBase *externalPotentialPtr = setup.externalPotential();

    /* This functionality allows us to use the code directly to generate
     * graphene potential binary or plain-text lookup tables (and convert
     * between them) */ 
    std::vector<std::string> grapheneLUTOptions= {"graphenelut3dtobinary", "graphenelut3dtotext","graphenelut3dgenerate"};
    if (std::find(grapheneLUTOptions.begin(),grapheneLUTOptions.end(), constants<std::string>("external")) != grapheneLUTOptions.end())
            return 99;


    /* Get the initial conditions associated with the external potential */
    /* Must use the copy constructor as we return a copy */
    DynamicArray<dVec,1> initialPos = 
        externalPotentialPtr->initialConfig(boxPtr,random,constants()->initialNumParticles());

    /* Perform a classical canonical pre-equilibration to obtain a suitable
     * initial state if N > 0*/
    if (!constants()->restart() && (constants()->initialNumParticles() > 0) ) {
       ClassicalMonteCarlo CMC(externalPotentialPtr,interactionPotentialPtr,random,boxPtr,
               initialPos);
       CMC.run(constants()->numEqSteps(),0);
    }

    /* Setup the path data variable */
    boost::ptr_vector<Path> pathPtrVec;
    for(int i=0; i<Npaths; i++){
        pathPtrVec.push_back(
                new Path(boxPtr,lookupPtrVec[i],constants()->numTimeSlices(),
                         initialPos,constants()->numBroken()));
    }
    
    //PotentialBase *externalPotentialPtr = setup.externalPotential();
    /* The Trial Wave Function (constant for pimc) */
    WaveFunctionBase *waveFunctionPtr = setup.waveFunction(pathPtrVec.front(),lookupPtrVec.front());

    /* Setup the action */
    boost::ptr_vector<ActionBase> actionPtrVec;
    for(int i=0; i<Npaths; i++){
        actionPtrVec.push_back(
                setup.action(pathPtrVec[i],lookupPtrVec[i],externalPotentialPtr,
                             interactionPotentialPtr,waveFunctionPtr) );
    }

    /* The list of Monte Carlo updates (moves) that will be performed */
    boost::ptr_vector< boost::ptr_vector<MoveBase> > movesPtrVec;
    for(int i=0; i<Npaths;i++){
        movesPtrVec.push_back(
                setup.moves(pathPtrVec[i],&actionPtrVec[i],random));
    }

    /* The list of estimators that will be performed */
    boost::ptr_vector< boost::ptr_vector<EstimatorBase> > estimatorsPtrVec;
    for(int i=0; i<Npaths;i++){
        estimatorsPtrVec.push_back(
                setup.estimators(pathPtrVec[i],&actionPtrVec[i],random));

        /* Add labels to estimator output files for multiple paths */
        if(i > 0) {
            for(uint32 j = 0; j < estimatorsPtrVec.back().size(); j++)
                estimatorsPtrVec.back().at(j).appendLabel(str(format("%d") % (i+1)));
        }
    }

    /* Setup the multi-path estimators */
    if(Npaths>1){
        estimatorsPtrVec.push_back(setup.estimators(pathPtrVec,actionPtrVec,random));
    }

    /* Setup the pimc object */
    PathIntegralMonteCarlo pimc(pathPtrVec,random,movesPtrVec,estimatorsPtrVec,
                                !setup.params["start_with_state"].as<std::string>().empty());

    /* A silly banner */
    if (PIGS)
        std::cout << std::endl 
             << " _____    _____    _____    _____"   << std::endl 
             << "|  __ \\  |_   _|  / ____|  / ____|" << std::endl
             << "| |__) |   | |   | |  __  | (___"    << std::endl
             << "|  ___/    | |   | | |_ |  \\___ \\" << std::endl
             << "| |       _| |_  | |__| |  ____) |"  << std::endl
             << "|_|      |_____|  \\_____| |_____/"  << std::endl
             << std::endl;  
    else 
        std::cout << std::endl
             << "  _____    _____   __  __    _____"    << std::endl
             << " |  __ \\  |_   _| |  \\/  |  / ____|" << std::endl
             << " | |__) |   | |   | \\  / | | |     "  << std::endl
             << " |  ___/    | |   | |\\/| | | |     "  << std::endl
             << " | |       _| |_  | |  | | | |____ "   << std::endl
             << " |_|      |_____| |_|  |_|  \\_____|"  << std::endl
             << std::endl;

    /* If this is a fresh run, we equilibrate and output simulation parameters to disk */
    if (!constants()->restart()) {

        /* Equilibrate */
        std::cout << format("[PIMCID: %s] - Pre-Equilibration Stage.") % constants()->id() << std::endl;
        for (uint32 n = 0; n < constants()->numEqSteps(); n++) 
            pimc.equilStep(n,setup.params("relax"),setup.params("relaxmu"));

        /* If we have relaxed the chemical potential, need to update file names 
         * in the grand canonical ensemble */
        if (!setup.params("canonical") && setup.params("relaxmu")) {
            communicate()->updateNames();
        }

        /* Output simulation details/parameters */
        setup.outputOptions(argc,argv,seed,boxPtr,lookupPtrVec.front().getNumNNGrid());
    }

    std::cout << format("[PIMCID: %s] - Measurement Stage.") % constants()->id() << std::endl;

    /* Sample */
    int oldNumStored = 0;
    int outNum = 0;
    int numOutput = setup.params["output_config"].as<int>();
    uint32 n = 0;
    do {
        pimc.step();
        if (pimc.numStoredBins > oldNumStored) {
            oldNumStored = pimc.numStoredBins;
            std::cout << format("[PIMCID: %s] - Bin #%5d stored to disk.") % constants()->id() 
                % oldNumStored << std::endl;
        }
        n++;

        /* Output configurations to disk */
        if ((numOutput > 0) && ((n % numOutput) == 0)) {
            pathPtrVec.front().outputConfig(outNum);
            outNum++;
        }
        
        /* Check if we've reached the wall clock limit*/
        if(constants()->wallClockOn()){
            current_time = time(NULL);
            if ( uint32(current_time)  > (uint32(start_time) + constants()->wallClock()) ){
                wallClockReached = true;
                break;
            }
        }
    } while (pimc.numStoredBins < setup.params["number_bins_stored"].as<int>());
    if (wallClockReached)
        std::cout << format("[PIMCID: %s] - Wall clock limit reached.") % constants()->id() << std::endl;
    else
        std::cout << format("[PIMCID: %s] - Measurement complete.") % constants()->id() << std::endl;

    /* Output Results */
    if (!constants()->saveStateFiles())
        pimc.saveState(1);
    pimc.finalOutput();

    /* Free up memory */
    delete interactionPotentialPtr;
    delete externalPotentialPtr;
    delete boxPtr;
    delete waveFunctionPtr;

    return 1;
}
