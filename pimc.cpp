/**
 * @file pimc.cpp
 * @author Adrian Del Maestro
 * 
 * @brief PathIntegralMonteCarlo class implementation.
 */

#include "pimc.h"
#include "path.h"
#include "lookuptable.h"
#include "move.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PATH INTEGRAL MONTE CARLO CLASS -------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
*  Constructor.
*
*  Here we initialize all data structures, moves and estimators that will be
*  required with peforming a path integral quantum monte carlo simulation.
*  The initialization depends on whether or not we are restarting, or starting
*  from a user supplied state.
******************************************************************************/
PathIntegralMonteCarlo::PathIntegralMonteCarlo (boost::ptr_vector<Path> &_pathPtrVec,
        MTRand &_random, boost::ptr_vector<move_vector> &_movePtrVec,
        boost::ptr_vector<estimator_vector> &_estimatorPtrVec,
        const bool _startWithState,uint32 _binSize) :
    random(_random),
    binSize(_binSize),
    Npaths(_pathPtrVec.size()),
    pathPtrVec(_pathPtrVec),
    path(pathPtrVec.front()),
    movePtrVec(_movePtrVec),
    move(movePtrVec.front()),
    estimatorPtrVec(_estimatorPtrVec),
    estimator(estimatorPtrVec.front())
{
	/* Are we starting from a saved state? */
	startWithState = _startWithState;
    
    /* Initialize stateStrings */
    stateStrings.resize(Npaths);
    
	/* We keep simulating until we have stored a certain number of measurements */
	numStoredBins = 0;
    
    /* Initialize the number of sweeps and updates per sweep to zero */
    numSteps = 0;
    numUpdates = 0;
    
	/* Calculate the number of sweeps to make sure we touch every bead */
	numImagTimeSweeps = int(ceil((1.0*constants()->numTimeSlices()/(1.0*constants()->Mbar()))));
    
	/* These counters are used in the equilibration process */
	numConfig = 0;
	numDiagonal = 0;
	numCoMAttempted = 200;
	prevNumCoMAttempted = 200;
	numCoMAccepted = 0;
	numDisplaceAttempted = 200;
	numDisplaceAccepted = 0;
    numMuAttempted = 5000;
    numNAttempted = 0;
    cN = 0;
    
	/* Have we saved a state? */
	savedState  = constants()->restart();
    
	/* Intialize the config number to zero */
	configNumber = 0;
    
	/* Determine the cumulative attempt move probabilities, and the indices of
     * various diagonal moves */
    double cumDiagProb = 0.0;
    double cumOffDiagProb = 0.0;
    string moveName;
	for (move_vector::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {

        /* Get the namne of the move and check if it is the generic diagonal
         * move */
        moveName = movePtr->getName();
        if ((moveName == "bisection") || (moveName == "staging"))
            moveName = "diagonal";

        /* Accumulate the diagonal moves */
        if ( (movePtr->operateOnConfig == DIAGONAL) || movePtr->operateOnConfig == ANY ) {
            cumDiagProb += constants()->attemptProb(moveName);
            attemptDiagProb.push_back(cumDiagProb);
        }
        else 
            attemptDiagProb.push_back(cumDiagProb);

        /* Accumulate the off-diagonal moves */
        if ( (movePtr->operateOnConfig == OFFDIAGONAL) || movePtr->operateOnConfig == ANY) {
            cumOffDiagProb += constants()->attemptProb(moveName);
            attemptOffDiagProb.push_back(cumOffDiagProb);
        }
        else
            attemptOffDiagProb.push_back(cumOffDiagProb);

        /* Find the indices of moves in case we need them */
        moveIndex[moveName] = std::distance(move.begin(), movePtr);
    }
    
    /* Make sure the move cumulative probability arrays add up to 1 */
    attemptDiagProb.back() = 1.0 + EPS;
    attemptOffDiagProb.back() = 1.0 + EPS;
    
	/* If we are restarting, or loading a state from disk, do so now. */
	if (startWithState || constants()->restart())
		loadState();

    /* Setup all the estimators for measurement i/o */
    for (boost::ptr_vector<estimator_vector>::iterator estimatorPtrVecItr = estimatorPtrVec.begin();
         estimatorPtrVecItr != estimatorPtrVec.end(); ++estimatorPtrVecItr) {
        for (estimator_vector::iterator estimatorPtr = estimatorPtrVecItr->begin();
             estimatorPtr != estimatorPtrVecItr->end(); ++estimatorPtr)
            estimatorPtr->prepare();
    }

    /* Make a list of estimator names for the 0th estimator */
    for (estimator_vector::iterator estimatorPtr = estimator.begin();
            estimatorPtr != estimator.end(); ++estimatorPtr) 
        estimatorIndex[estimatorPtr->getName()] = std::distance(estimator.begin(), estimatorPtr);
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
PathIntegralMonteCarlo::~PathIntegralMonteCarlo () {
    // empty destructor
}

/**************************************************************************//**
 *  Performs one of various possible Monte Carlo updats on the path.
******************************************************************************/
string PathIntegralMonteCarlo::update(const double x, const int sweep, const int pathIdx=0) { 

    success = false;
	string moveName = "NONE";
    int index;

    /* If we have no beads, all we can do is insert */
    /* N.B. This is not multi-path safe!!! */
    if (path.worm.getNumBeadsOn() == 0) {
        index = moveIndex["insert"];
        success = move.at(index).attemptMove();
        moveName = move.at(index).getName();
        return moveName;
    }

    /* Determine the index of the move to be performed */
	if (pathPtrVec[pathIdx].worm.isConfigDiagonal)
        index = std::lower_bound(attemptDiagProb.begin(),attemptDiagProb.end(),x)
            - attemptDiagProb.begin();
    else 
        index = std::lower_bound(attemptOffDiagProb.begin(),attemptOffDiagProb.end(),x)
            - attemptOffDiagProb.begin();

    /* Perform the move */
    moveName = movePtrVec[pathIdx].at(index).getName();

    if (moveName == "center of mass"){
        if ( (numParticles > 0) && ((sweep % numParticles)== 0) )
            success = movePtrVec[pathIdx].at(index).attemptMove();
    }
    else
        success = movePtrVec[pathIdx].at(index).attemptMove();

	return moveName;
}

/**************************************************************************//**
 *  Equilibration.
 * 
 *  The equilibration method, where we perform fully diagonal moves 1/3 of the
 *  time, finding an optimal, COM step, then for another 1/3 we perform 
 *  all moves adjusting C0 to yield a diagonal fraction near 0.75 and and for 
 *  the final 1/3 we just equilibrate.
******************************************************************************/
void PathIntegralMonteCarlo::equilStep(const uint32 iStep, const bool relaxC0, const bool relaxmu) {

	double x;

	/* How far are we into the equilibration? */
	double equilFrac = (1.0*iStep) / (1.0*constants()->numEqSteps());
    numParticles = path.getTrueNumParticles();

    /* Always do at least one move */
    if (numParticles == 0)
        ++numParticles;

	/* For the first 1/3 of equilibration steps, we only do diagonal moves */
	if (equilFrac < 1.0/3.0 && !startWithState) {

        for (int n = 0; n < numParticles; n++) {
            
            x = random.rand();
            
            /* Here we do the diagonal pre-equilibration moves, and allow for 
             * optimization of simulation parameters */
            if (x < constants()->attemptProb("center of mass")) {

                /* Get the move index for the center of mass */
                int index = moveIndex["center of mass"];

                for( boost::ptr_vector<move_vector>::iterator movePtrVecItr=movePtrVec.begin();
                    movePtrVecItr!=movePtrVec.end();++movePtrVecItr){
                    movePtrVecItr->at(index).attemptMove();
                }

                /* We check how many CoM moves we have tried.  Every 200 moves, we see if we need
                 * to adjust comDelta, provided we are in the pre-equilibration diagonal state. */
                if ( (move.at(index).getNumAttempted() > 0) 
                        && (move.at(index).getNumAttempted() > prevNumCoMAttempted)
                        && (move.at(index).getNumAttempted() % numCoMAttempted == 0) 
                        && (constants()->comDelta() < 0.5*blitz::min(path.boxPtr->side)) ) {

                    numCoMAccepted  = move.at(index).getNumAccepted() - numCoMAccepted;
                    double CoMRatio = 1.0*numCoMAccepted / numCoMAttempted;
                    if (CoMRatio < 0.2)
                        constants()->shiftCoMDelta(-0.6);
                    else if (CoMRatio < 0.3)
                        constants()->shiftCoMDelta(-0.4);
                    else if (CoMRatio < 0.4)
                        constants()->shiftCoMDelta(-0.2);
                    else if (CoMRatio > 0.6)
                        constants()->shiftCoMDelta(0.2);
                    else if (CoMRatio > 0.7)
                        constants()->shiftCoMDelta(0.4);
                    else if (CoMRatio > 0.8)
                        constants()->shiftCoMDelta(0.6);

                    /* Reset the counters */
                    numCoMAccepted = move.at(index).getNumAccepted();
                    prevNumCoMAttempted = move.at(index).getNumAttempted();
                } // CoM Delta Shift

            } // Center of mass move
            /* Now try a displace move */
            else if (x < constants()->attemptProb("center of mass") + constants()->attemptProb("displace")) {

                int index = moveIndex["displace"];
                for( boost::ptr_vector<move_vector>::iterator movePtrVecItr=movePtrVec.begin();
                    movePtrVecItr!=movePtrVec.end();++movePtrVecItr){
                    movePtrVecItr->at(index).attemptMove();
                }

                /* We check how many displace moves we have tried.  Every numDisplaceAttempted moves, we see if we need
                 * to adjust delta, provided we are in the pre-equilibration diagonal state. */
                if ( (move.at(index).getNumAttempted() > 0) 
                        && (move.at(index).getNumAttempted() % numDisplaceAttempted == 0) ) {

                    numDisplaceAccepted = move.at(index).getNumAccepted() - numDisplaceAccepted;
                    double displaceRatio = 1.0*numDisplaceAccepted / numDisplaceAttempted;
                    if (displaceRatio < 0.2)
                        constants()->shiftDisplaceDelta(-0.6);
                    else if (displaceRatio < 0.3)
                        constants()->shiftDisplaceDelta(-0.4);
                    else if (displaceRatio < 0.4)
                        constants()->shiftDisplaceDelta(-0.2);
                    else if (displaceRatio > 0.6)
                        constants()->shiftDisplaceDelta(0.2);
                    else if (displaceRatio > 0.7)
                        constants()->shiftDisplaceDelta(0.4);
                    else if (displaceRatio > 0.8)
                        constants()->shiftDisplaceDelta(0.6);

                    //cout << "delta: " << constants()->delta() << " " << displaceRatio << endl;
                    /* Reset the counters */
                    numDisplaceAccepted = move.at(index).getNumAccepted();
                }
            } // displace move
            else {
                /* Attemp a diagonal path update*/
                for (int sweep = 0; sweep < numImagTimeSweeps; sweep++){
                    for( boost::ptr_vector<move_vector>::iterator movePtrVecItr=movePtrVec.begin();
                        movePtrVecItr!=movePtrVec.end();++movePtrVecItr){
                        movePtrVecItr->at(moveIndex["diagonal"]).attemptMove();
                    }
                }
            }

        } // numParticles

    } // equilFrac < 1/3 and !startWithState

    /* Start the 2/3 off diagonal portion of the equilibration */
	else {

        /* How many updates should we perform? We always try atleast 1 */
        numUpdates = max(1,int(ceil(path.worm.getNumBeadsOn()/(1.0*constants()->Mbar()))));

        /* Increment the number of off-diagonal steps */
        numSteps++;

        for (int n = 0; n < numUpdates; n++) {

            /* Relax the chemical potential to get a desired average number of
             * particles */
            if (equilFrac < 2.0/3.0) {
                cN += path.getTrueNumParticles();
                numNAttempted += 1;
                if ( (relaxmu) && (numNAttempted == numMuAttempted) ) {

                    double aveN = cN / (1.0*numNAttempted);

                    /* The factor that we shift the chemical potenital by */
                    double factor = (1.0 - aveN/constants()->initialNumParticles());

                    cout << format("mu: %8.5E\t%8.5E\t") % constants()->mu() % factor;
                    constants()->shiftmu(factor);
                    cout << format("%8.5E\t%8.5f\t%5d\n") % constants()->mu() 
                        % aveN % path.getTrueNumParticles();

                    cN = 0;
                    numNAttempted = 0;
                }
            }

            /* Generate random number and run through all moves */
            x = random.rand();
            string mName;
            for(uint32 p=0;p<Npaths;p++){
                mName = update(x,n,p);
            }

            /* We accumulate the number of diagonal configurations */
            numConfig++;
            if (path.worm.isConfigDiagonal) 
                ++numDiagonal;

            /* cout << mName << "\t" << numDiagonal << endl; */

            /* Every 200 steps, we check the diagonal fraction and update the
             * worm constant for the first 1/3 of equilibration */
            if ( (relaxC0) && (numSteps > 200) && (equilFrac < 2.0/3.0) ) {

                double diagFrac = 1.0*numDiagonal / (1.0*numConfig);

                cout << format("%4.2f\t%8.5f\t") % diagFrac % constants()->C0();

                /* Adjust the worm constant to ensure that the diagonal fraction
                 * is between 0.75 and 0.85, but don't let it get too small or too large*/
                if (constants()->C0() > 1.0E-5) {
                    if (diagFrac < 0.2)
                        constants()->shiftC0(-0.5);
                    else if (diagFrac >= 0.2 && diagFrac < 0.3)
                        constants()->shiftC0(-0.4);
                    else if (diagFrac >= 0.3 && diagFrac < 0.4)
                        constants()->shiftC0(-0.3);
                    else if (diagFrac >= 0.4 && diagFrac < 0.5)
                        constants()->shiftC0(-0.2);
                    else if (diagFrac >= 0.5 && diagFrac < 0.6)
                        constants()->shiftC0(-0.1);
                    else if (diagFrac >= 0.6 && diagFrac < 0.75)
                        constants()->shiftC0(-0.05);
                }

                if (constants()->C0() < 1.0E4) {
                    if (diagFrac <= 0.9 && diagFrac > 0.85)
                        constants()->shiftC0(0.05);
                    else if (diagFrac <= 0.95 && diagFrac > 0.9)
                        constants()->shiftC0(0.1);
                    else if (diagFrac > 0.95)
                        constants()->shiftC0(0.2);
                }

                cout << format("%8.5f\t%5d\t%8.6f\n") % constants()->C0() 
                    % path.getTrueNumParticles() 
                    % (1.0*path.getTrueNumParticles()/path.boxPtr->volume);

                /* Reset the counters */
                numDiagonal = 0;
                numConfig = 0;
                numSteps = 0;

            } // relaxC0

        } // numUpdates 

    } // 2/3 off-diagonal equilibration

    /* Save a state to disk every 200 equilibrium steps */
    if ((constants()->saveStateFiles()) && (iStep > 0) && (iStep % 200) == 0)
        saveState();
}

/**************************************************************************//**
 *  PIMC step.
 * 
 *  This method performs the metropolis sampling, which is actually a 
 *  complicated multi-step operation which consists of various types of moves.  
 *  They can in general occur at different frequencies. We also measure all 
 *  estimators.
******************************************************************************/
void PathIntegralMonteCarlo::step() {

	string moveName;

    for (uint32 pIdx=0; pIdx<Npaths; pIdx++) {
    
        numUpdates = max(1,int(ceil(pathPtrVec[pIdx].worm.getNumBeadsOn()/(1.0*constants()->Mbar()))));

        /* We run through all moves, making sure that we could have touched each bead at least once */
        for (int n = 0; n < numUpdates ; n++) 
            moveName = update(random.rand(),n,pIdx);

        numParticles = pathPtrVec[pIdx].getTrueNumParticles();
        if (numParticles > 0) {
        /* Perform all measurements */
        for (estimator_vector::iterator estimatorPtr = estimatorPtrVec[pIdx].begin();
             estimatorPtr != estimatorPtrVec[pIdx].end(); ++estimatorPtr)
            estimatorPtr->sample();
        }
        
        /* Every binSize measurements, we output averages to disk and record the
         * state of the simulation on disk.  */
        if (estimatorPtrVec[pIdx].size() > 0){
            if (estimatorPtrVec[pIdx].front().getNumAccumulated() >= binSize) {
                for (estimator_vector::iterator estimatorPtr = estimatorPtrVec[pIdx].begin();
                     estimatorPtr != estimatorPtrVec[pIdx].end(); ++estimatorPtr) {
                    if (estimatorPtr->getNumAccumulated() >= binSize)
                        estimatorPtr->output();
                }
                if(Npaths==1){
                    if (constants()->saveStateFiles())
                        saveState();
                    else
                        saveStateString();
                }
                if (pIdx == 0){
                    ++numStoredBins;
                }
            }
        }
    }
        
    //////Multi-path Estimators////////////
    if(estimatorPtrVec.size() > Npaths){
        for (estimator_vector::iterator estimatorPtr = estimatorPtrVec.back().begin();
                estimatorPtr != estimatorPtrVec.back().end(); ++estimatorPtr)
            estimatorPtr->sample();

        /* Every binSize measurements, we output averages to disk and record the
         * state of the simulation on disk.  */
        if (estimatorPtrVec.back().front().getNumAccumulated() >= binSize) {
            for (estimator_vector::iterator estimatorPtr = estimatorPtrVec.back().begin();
                    estimatorPtr != estimatorPtrVec.back().end(); ++estimatorPtr) {
                if (estimatorPtr->getNumAccumulated() >= binSize)
                    estimatorPtr->output();
            }

            /* Save to disk or store a state file */
            if (constants()->saveStateFiles())
                saveState();
            else
                saveStateString();

            if (estimator.size() == 0){
                ++numStoredBins;
            }
        }
    }
}

/**************************************************************************//**
 *  Output simulation statistics to disk.
 *
 *  We perform this once at the very end of the simulation.  It saves the 
 *  details of accetance probabilities as well as estimators to a file.
******************************************************************************/
void PathIntegralMonteCarlo::finalOutput() {

	/* Output the acceptance data to the log file */
	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << endl;

	communicate()->file("log")->stream() << "---------- Begin Acceptance Data ---------------" << endl;
	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << format("%-29s\t:\t%7.5f\n") % "Total Rate" 
		% move.front().getTotAcceptanceRatio();
	communicate()->file("log")->stream() << endl;

	/* Ouptut all the move acceptance information to disk */
	string moveName;
	for (move_vector::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		moveName = movePtr->getName();

		/* We only output levels for those moves which have a variable size */
        if (movePtr->variableLength) {
			for (int n = 0; n <= constants()->b(); n++) {
				communicate()->file("log")->stream() << format("%-12s Level %-10d\t:\t%7.5f\t(%d/%d)\n") 
					% moveName % n % movePtr->getAcceptanceRatioLevel(n) 
					% movePtr->numAcceptedLevel(n) % movePtr->numAttemptedLevel(n);
			}
		}
		communicate()->file("log")->stream() << format("%-29s\t:\t%7.5f\t(%d/%d)\n") % moveName
			% movePtr->getAcceptanceRatio() % movePtr->numAccepted
			% movePtr->numAttempted;
		communicate()->file("log")->stream() << endl;
	}
	communicate()->file("log")->stream() << "---------- End Acceptance Data -----------------" << endl;

	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << endl;

	/* Output the estimator statistics to the log file */
	communicate()->file("log")->stream() << "---------- Begin Estimator Data ----------------" << endl;
	communicate()->file("log")->stream() << endl;
	for (estimator_vector::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->file("log")->stream() << format("%-29s\t:\t%16d\t%16d\n") % estimatorPtr->getName()
			% estimatorPtr->getNumSampled() % estimatorPtr->getTotNumAccumulated();

	}
	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << "---------- End Estimator Data ------------------" << endl;
}

/**************************************************************************//**
 *  Save the state of the simulation to disk, including all path, worm,
 *  move and estimator data.
******************************************************************************/
void PathIntegralMonteCarlo::saveState() {

	savedState = true;
    
    string stateFileName;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++){
        
        stateFileName = "state";
        if(pIdx > 0){
            stateFileName += str(format("%d") % (pIdx+1));
        }
        
        /* Prepare the state file for writing */
        communicate()->file(stateFileName.c_str())->reset();

        /* We First write the current total number of world lines */
        communicate()->file(stateFileName.c_str())->stream() << path.getNumParticles() << endl;

        /* Now write the total acceptance information for all moves */
        communicate()->file(stateFileName.c_str())->stream() << format("%16d\t%16d\n")
            % move.front().totAccepted % move.front().totAttempted;

        /* Now record the individual move acceptance information,
         * first for the diagonal, then off-diagonal*/
        for (move_vector::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
            communicate()->file(stateFileName.c_str())->stream() << format("%16d\t%16d\n")
                % movePtr->numAccepted % movePtr->numAttempted;
        }

        /* Output the estimator sampling information */
        for (estimator_vector::iterator estimatorPtr = estimator.begin();
                estimatorPtr != estimator.end(); ++estimatorPtr) {
            communicate()->file(stateFileName.c_str())->stream() << format("%16d\t%16d\n")
                % estimatorPtr->getTotNumAccumulated()
                % estimatorPtr->getNumSampled();
        }

        /* Now we output the actual path and worldline data */
        communicate()->file(stateFileName.c_str())->stream() << setprecision(16)<< path.beads << endl;
        communicate()->file(stateFileName.c_str())->stream() << path.nextLink << endl;
        communicate()->file(stateFileName.c_str())->stream() << path.prevLink << endl;

        /* Output the worm data */
        communicate()->file(stateFileName.c_str())->stream() << path.worm.beads << endl;

        /* Save the state of the random number generator */
        uint32 randomState[random.SAVE];
        random.save(randomState);
        for (int i = 0; i < random.SAVE; i++)
            communicate()->file(stateFileName.c_str())->stream() << randomState[i] << " ";
        communicate()->file(stateFileName.c_str())->stream() << endl;

        /* Rename and copy the file. */
        communicate()->file(stateFileName.c_str())->rename();
    }
}

/**************************************************************************//**
*  Save the state of the simulation to disk, including all path, worm,
*  move and estimator data.
******************************************************************************/
void PathIntegralMonteCarlo::saveStateString() {

    stringstream stateStrStrm;
    
    for( uint32 pIdx=0; pIdx<Npaths; pIdx++){
        
        stateStrStrm.str("");
        
        /* We First write the current total number of world lines */
        stateStrStrm << pathPtrVec[pIdx].getNumParticles() << endl;
        
        /* Now write the total acceptance information for all moves */
        stateStrStrm << format("%16d\t%16d\n")
            % movePtrVec[pIdx].front().totAccepted % movePtrVec[pIdx].front().totAttempted;
        
        /* Now record the individual move acceptance information,
         * first for the diagonal, then off-diagonal*/
        for (move_vector::iterator movePtr = movePtrVec[pIdx].begin(); movePtr != movePtrVec[pIdx].end(); ++movePtr) {
            stateStrStrm << format("%16d\t%16d\n")
                % movePtr->numAccepted % movePtr->numAttempted;
        }
        
        /* Output the estimator sampling information */
        for (estimator_vector::iterator estimatorPtr = estimatorPtrVec[pIdx].begin();
                estimatorPtr != estimatorPtrVec[pIdx].end(); ++estimatorPtr) {
            stateStrStrm << format("%16d\t%16d\n") % estimatorPtr->getTotNumAccumulated()
                % estimatorPtr->getNumSampled();
        }
        
        /* Now we output the actual path and worldline data */
        stateStrStrm << setprecision(16) << pathPtrVec[pIdx].beads << endl;
        stateStrStrm << pathPtrVec[pIdx].nextLink << endl;
        stateStrStrm << pathPtrVec[pIdx].prevLink << endl;
        
        /* Output the worm data */
        stateStrStrm << pathPtrVec[pIdx].worm.beads << endl;
        
        /* Save the state of the random number generator */
        uint32 randomState[random.SAVE];
        random.save(randomState);
        for (int i = 0; i < random.SAVE; i++)
            stateStrStrm << randomState[i] << " ";
        stateStrStrm << endl;
        
        stateStrings[pIdx] = stateStrStrm.str();
    }

}

/**************************************************************************//**
*  Save the state of the simulation to disk, including all path, worm,
*  move and estimator data.
******************************************************************************/
void PathIntegralMonteCarlo::saveStateFromStr() {
    
	savedState = true;
    string stateFileName;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++){
        
        stateFileName = "state";
        if(pIdx > 0){
            stateFileName += str(format("%d") % (pIdx+1));
        }
        
        /* Prepare the state file for writing */
        communicate()->file(stateFileName.c_str())->reset();
        
        /* Write the stateString to disk */
        communicate()->file(stateFileName.c_str())->stream() << stateStrings[pIdx];
        
        /* Rename and copy the file. */
        communicate()->file(stateFileName.c_str())->rename();
        
    }
    
}

/**************************************************************************//**
 *  Load a classical ground state from file.
******************************************************************************/
void PathIntegralMonteCarlo::loadClassicalState(Array <dVec,2> &tempBeads,
        Array <unsigned int, 2> &tempWormBeads, int numWorldLines) {

    /* We go through each active worldline and create a new classical
     * configuration */
    int ptcl;
    beadLocator beadIndex;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++){
        ptcl = 0;
        for (int n = 0; n < numWorldLines; n++) {
            beadIndex = 0,n;

            /* Check if the bead is on */
            if (tempWormBeads(beadIndex)) {
                
                /* Assign the classical configuration */
                pathPtrVec[pIdx].beads(Range::all(),ptcl) = tempBeads(beadIndex);
                pathPtrVec[pIdx].worm.beads(Range::all(),ptcl) = 1;
                ptcl++;
            }
        }
    }
}

/**************************************************************************//**
 *  Load a quantum ground state from file.
******************************************************************************/
void PathIntegralMonteCarlo::loadQuantumState(Array <dVec,2> &tempBeads, 
        Array <beadLocator,2> &tempNextLink, Array <beadLocator,2> &tempPrevLink,
        int numTimeSlices, int tempNumWorldLines) {

    /* Prevent double counting of worldlines */
    Array <bool, 1> doBead(tempNumWorldLines);

    beadLocator startBead,beadIndex;
    beadLocator newStartBead;
    int ptcl;
    int slice;
    
    for(uint32 pIdx=0; pIdx<Npaths; pIdx++) {

        newStartBead = 0,0;
        ptcl = 0;
        slice = 0;
        doBead = true;
        
        /* Now we iterate through each worldline exactly once */
        for (int n = 0; n < tempNumWorldLines; n++) {

            /* The initial bead to be moved */
            startBead = 0,n;

            /* We make sure we don't try to touch the same worldline twice */
            if (doBead(n)) {

                beadIndex = startBead;

                /* The world line length, we simply advance until we have looped back on 
                 * ourselves. */
                slice = 0;
                newStartBead = slice,ptcl;
                do {
                    /* We turn off any zero-slice beads we have touched */
                    if (beadIndex[0]==0)
                        doBead(beadIndex[1]) = false;

                    pathPtrVec[pIdx].beads(slice % numTimeSlices,ptcl) = tempBeads(beadIndex);
                    pathPtrVec[pIdx].worm.beads(slice % numTimeSlices,ptcl) = 1;

                    beadIndex = tempNextLink(beadIndex);
                    ++slice;

                    /* Do a forward reconnection, provided we are not at the
                     * last bead */
                    if ( ((slice % numTimeSlices) == 0) && !all(beadIndex==startBead)) {
                        pathPtrVec[pIdx].nextLink(numTimeSlices-1,ptcl) = 0,ptcl+1;
                        pathPtrVec[pIdx].prevLink(0,ptcl+1) = numTimeSlices-1,ptcl;
                        ++ptcl;
                    }
                } while (!all(beadIndex==startBead));

                /* Now we have to add the remaining beads and perform the final
                 * reconnection */
                for (int tslice = (slice % numTimeSlices); tslice < numTimeSlices; tslice++) {
                    pathPtrVec[pIdx].beads(tslice,ptcl) = tempBeads(beadIndex);
                    pathPtrVec[pIdx].worm.beads(tslice,ptcl) = 1;
                }
                pathPtrVec[pIdx].nextLink(numTimeSlices-1,ptcl) = newStartBead;
                pathPtrVec[pIdx].prevLink(newStartBead) = numTimeSlices-1,ptcl;
                ++ptcl;

            } // doBead
        } // n
    } // pIdx

    /* Free local memory */
    doBead.free();

}

/**************************************************************************//**
 *  Load a previous state of the simulation from disk.
******************************************************************************/
void PathIntegralMonteCarlo::loadState() {

    string tempString;
    
    string fileInitStr = "init";
    
    for( uint32 pIdx=0; pIdx<Npaths; pIdx++){
        
        if(pIdx>0){
            fileInitStr = str(format("init%d") % (pIdx+1));
        }

        /* Reset the total acceptance information */
        movePtrVec[pIdx].front().resetTotAccept();

        /* Reset all the individual move acceptance information */
        for (move_vector::iterator movePtr = movePtrVec[pIdx].begin(); movePtr != movePtrVec[pIdx].end(); ++movePtr) {
            movePtr->resetAccept();
        }

        /* Reset estimator sampling information */
        for (estimator_vector::iterator estimatorPtr = estimatorPtrVec[pIdx].begin();
                estimatorPtr != estimatorPtrVec[pIdx].end(); ++estimatorPtr) {
            estimatorPtr->restart(0,0);
        }

        /* We first read the former total number of world lines */
        int numWorldLines;
        int numTimeSlices = pathPtrVec[pIdx].numTimeSlices;
        communicate()->file(fileInitStr)->stream() >> numWorldLines;

        /* Now we skip through the input file until we find the beads matrix.  This
         * is signalled by the appearance of an open bracket "(" */
        while (!communicate()->file(fileInitStr)->stream().eof()) {
            if (communicate()->file(fileInitStr)->stream().peek() != '(') 
                getline (communicate()->file(fileInitStr)->stream(), tempString);
            else
                break;
        }

        /* Now we resize all path data members and read them from the init state file */
        pathPtrVec[pIdx].beads.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].nextLink.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].prevLink.resize(numTimeSlices,numWorldLines);
        pathPtrVec[pIdx].worm.beads.resize(numTimeSlices,numWorldLines);

        /* A temporary container for the beads array */
        Array <dVec,2> tempBeads;

        /* Get the worldline configuration */
        communicate()->file(fileInitStr)->stream() >> tempBeads;

        /* The temporary number of time slices */
        int tempNumTimeSlices = tempBeads.rows();

        if (tempNumTimeSlices == numTimeSlices) {
        
            /* Copy over the beads array */
            pathPtrVec[pIdx].beads = tempBeads;

            /* Get the link arrays */
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].nextLink;
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].prevLink;

            /* Repeat for the worm file */
            communicate()->file(fileInitStr)->stream() >> pathPtrVec[pIdx].worm.beads;

        } // locBeads.rows() == numTimeSlices
        else {

            /* Initialize the links */
            firstIndex i1;
            secondIndex i2;
            pathPtrVec[pIdx].prevLink[0] = i1-1;
            pathPtrVec[pIdx].prevLink[1] = i2;
            pathPtrVec[pIdx].nextLink[0] = i1+1;
            pathPtrVec[pIdx].nextLink[1] = i2;
        
            /* Here we implement the initial periodic boundary conditions in 
             * imaginary time */
            pathPtrVec[pIdx].prevLink(0,Range::all())[0] = numTimeSlices-1;
            pathPtrVec[pIdx].nextLink(numTimeSlices-1,Range::all())[0] = 0;

            /* Reset the worm.beads array */
            pathPtrVec[pIdx].worm.beads = 0;

            /* Temporary containers for the links and worm beads */
            Array <beadLocator,2> tempNextLink;
            Array <beadLocator,2> tempPrevLink;
            Array <unsigned int,2> tempWormBeads;

            /* Get the link arrays and worm file */
            communicate()->file(fileInitStr)->stream() >> tempNextLink;
            communicate()->file(fileInitStr)->stream() >> tempPrevLink;
            communicate()->file(fileInitStr)->stream() >> tempWormBeads;

            /* Load a classical (all time slice positions equal) from the input
             * file */
            loadClassicalState(tempBeads,tempWormBeads, numWorldLines);

            /* Load a quantum initial state from a file */
            //if (tempNumTimeSlices < numTimeSlices) {
            //    loadQuantumState(tempBeads,tempNextLink,tempPrevLink,
            //            numTimeSlices,int(sum(tempWormBeads)/tempNumTimeSlices));
            //}

            /* Now we make sure all empty beads are unlinked */
            beadLocator beadIndex;
            for (int slice = 0; slice < numTimeSlices; slice++) {
                for (int ptcl = 0; ptcl < numWorldLines; ptcl++) {
                    beadIndex = slice,ptcl;
                    if (!pathPtrVec[pIdx].worm.beads(beadIndex)) {
                        pathPtrVec[pIdx].nextLink(beadIndex) = XXX;
                        pathPtrVec[pIdx].prevLink(beadIndex) = XXX;
                    }
                }
            }

            /* Free local memory */
            tempPrevLink.free();
            tempNextLink.free();
            tempWormBeads.free();

        } // locBeads.rows() != numTimeSlices

        /* Load the state of the random number generator, only if we are restarting 
         * the simulation */
        if (constants()->restart()) {
            uint32 randomState[random.SAVE];
            for (int i = 0; i < random.SAVE; i++) 
                communicate()->file(fileInitStr)->stream() >> randomState[i];
            random.load(randomState);
        }

        /* Reset the number of on beads */
        pathPtrVec[pIdx].worm.resetNumBeadsOn();

        /* Left pack the input data array */
        pathPtrVec[pIdx].leftPack();

        /* Go through all beads, and make sure they fit inside the simulation cell.
         * At the same time, determine how many active beads there are per slice */
        beadLocator beadIndex;
        pathPtrVec[pIdx].numBeadsAtSlice = 0;
        for (beadIndex[0] = 0; beadIndex[0] < numTimeSlices; ++beadIndex[0]) {
            for (beadIndex[1] = 0; beadIndex[1] < numWorldLines; ++beadIndex[1]) {
                pathPtrVec[pIdx].boxPtr->putInside(path(beadIndex));
                if (pathPtrVec[pIdx].worm.beadOn(beadIndex))
                    ++pathPtrVec[pIdx].numBeadsAtSlice(beadIndex[0]);
            } // particles
        } // slices

        /* Reset the worm parameters */
        pathPtrVec[pIdx].worm.isConfigDiagonal = true;
        pathPtrVec[pIdx].worm.reset();

        /* Resize and update the lookup table arrays */
        pathPtrVec[pIdx].lookup.resizeList(numWorldLines);
        pathPtrVec[pIdx].lookup.updateGrid(path);

        /* Reset the broken/closed worldline vectors */
        pathPtrVec[pIdx].resetBrokenClosedVecs();
        
        /* Close the file */
        communicate()->file(fileInitStr)->close();
        
        /* Free up memory */
        tempBeads.free();
        
    }
}

/**************************************************************************//**
 *  Output the worldline configuration to disk using PDB , suitable
 *  for plotting using vmd. 
 *  
 *  We must post-process the final pdb file and split it up due to connectivity 
 *  changes. 
 *  @see For the PDB specification: 
 *  http://www.wwpdb.org/documentation/format32/v3.2.html
******************************************************************************/
void PathIntegralMonteCarlo::outputPDB() {

	numParticles = path.getNumParticles();
	int numTimeSlices = path.numTimeSlices;

	configNumber++;

	/* We go through all beads, and find the start and end bead for each
	 * worldline, adding them to an array */
	Array <beadLocator,1> startBead,endBead;
	startBead.resize(numParticles);
	endBead.resize(numParticles);

	/* We sort the output by the number of beads in a worldline */
	Array <int,1> wlLength(numParticles);
	wlLength = 0;

	/* This is the index-beadNumber mapping array */
	Array <int,2> beadNum(numTimeSlices,numParticles);
	beadNum = 0;

	int numWorldLines = 0;

	/* Get the list of beads that are active in the simulation */
	Array <bool,2> doBead(numTimeSlices,numParticles);		
	doBead = cast<bool>(path.worm.getBeads());

	/* We go through each particle/worldline */
	int nwl = 0;
	int beadNumber = 0;
	for (int n = 0; n < numParticles; n++) {

		/* The initial bead to be moved */
		startBead(nwl) = 0,n;

		/* We make sure we don't try to touch the same worldline twice */
		if (doBead(startBead(nwl))) {

			/* Check to see if the start bead is on a worm.  If it is, we start
			 * at the worm tail and end at its head. */
			if (path.worm.foundBead(path,startBead(nwl))) {
				startBead(nwl) = path.worm.tail;
				endBead(nwl)   = path.next(path.worm.head);
			}
			/* Otherwise, we loop around until we find the initial bead */
			else
				endBead(nwl) = startBead(nwl);

			/* Mark the beads as touched and increment the number of worldlines */
			beadLocator beadIndex;
			beadIndex = startBead(nwl);
			int length = 1;
			do {
				doBead(beadIndex) = false;
				beadNum(beadIndex) = beadNumber;
				beadNumber++;
				length++;
				beadIndex = path.next(beadIndex);
			} while (!all(beadIndex==endBead(nwl)));

			/* We label each trajectory by the number of particles it contains.
			 * a worm is always given label 0 */
			if ((length % numTimeSlices) == 0)
				wlLength(nwl) = length/numTimeSlices;
			else 
				wlLength(nwl) = 0;

			nwl++;
		} // doBead

	} // n
	numWorldLines = nwl;

	/* Output the PDB header */
	communicate()->file("wl")->stream() << format("REMARK [CONFIG %04d]\n") % configNumber;

	/* Output the unit cell information.  It is always cubic.  Everything is scaled by
	 * an overall factor for better visualization. */

	double scale = 10.0;				
	int i;
	communicate()->file("wl")->stream() << format("%-6s") % "CRYST1";
	for (i = 0; i < NDIM; i++) 
		communicate()->file("wl")->stream() << format("%9.3f") % (scale*path.boxPtr->side[i]);
	while (i < 3) {
		communicate()->file("wl")->stream() << format("%9.3f") % 1.0;
		i++;
	}
	communicate()->file("wl")->stream() << format("%7.2f%7.2f%7.2f %-11s%4d\n") % 90.0 % 90.0 % 90.0 % "P 1" % 1;

	/* We output the atom block */
	beadLocator beadIndex;
	for (int n = 0; n < numWorldLines;  n++) {
		beadIndex = startBead(n);
		do {
			/* We give the zero-time-slice bead a special name */
			if (beadIndex[0] == 0) {
				communicate()->file("wl")->stream() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
					% beadNum(beadIndex) % "H0" % wlLength(n) % " ";
			}
			else {
				communicate()->file("wl")->stream() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
					% beadNum(beadIndex) % "HE" % wlLength(n) % " ";
			}

			/* Output the coordinates in 3D */
			for (i = 0; i < NDIM; i++) {
				communicate()->file("wl")->stream() << format("%8.3f") % (scale*path(beadIndex)[i]);
			}
			while (i < 3) {
				communicate()->file("wl")->stream() << format("%8.3f") % 0.0;
				i++;
			}
			communicate()->file("wl")->stream() << format("%14s\n") % "HE";

			beadIndex = path.next(beadIndex);
		} while (!all(beadIndex==endBead(n)));
	}
	communicate()->file("wl")->stream() <<("TER\n");

	/* Now output the connect block */
	for (int n = 0; n < numWorldLines;  n++) {
		beadIndex = startBead(n);
		do {
			communicate()->file("wl")->stream() << format("%-6s%5d") % "CONECT" % beadNum(beadIndex);
			beadLocator prevIndex,nextIndex;
			prevIndex = path.prev(beadIndex);
			nextIndex = path.next(beadIndex);
			
			/* We make sure that we don't connect beads linked by periodic bondary 
			 * conditions */

			/* Check the previous bead */
			if (path.worm.beadOn(prevIndex)) {
				dVec sep;
				sep = path(beadIndex) - path(prevIndex);
				if (dot(sep,sep) < path.boxPtr->rcut2)
					communicate()->file("wl")->stream() << format("%5d") % beadNum(prevIndex);
			}

			/* Now the next bead */
			if (path.worm.beadOn(nextIndex)) {
				dVec sep;
				sep = path(nextIndex) - path(beadIndex);
				if (dot(sep,sep) < path.boxPtr->rcut2)
					communicate()->file("wl")->stream() << format("%5d") % beadNum(nextIndex);
			}
			communicate()->file("wl")->stream() << endl;

			beadIndex = path.next(beadIndex);
		} while (!all(beadIndex==endBead(n)));
	}
	communicate()->file("wl")->stream() <<("END\n");

	/* Free up memory */
	startBead.free();
	endBead.free();
	wlLength.free();
	beadNum.free();
	doBead.free();
}

/**************************************************************************//**
 *  
 *  
******************************************************************************/
void PathIntegralMonteCarlo::printWormState() {

	/* We make a list of all the beads contained in the worm */
	Array <beadLocator,1> wormBeads;	// Used for debugging
	wormBeads.resize(path.worm.length+1);
	wormBeads = XXX;

	/* Output the worldline configuration */
	communicate()->file("debug")->stream() << " (" << path.getTrueNumParticles() << ")" << endl;
	communicate()->file("debug")->stream() << "head " << path.worm.head[0] << " " << path.worm.head[1]
        << " tail " << path.worm.tail[0] << " " << path.worm.tail[1]
        << " length " << path.worm.length 
        << " gap " << path.worm.gap << endl;

	if (!path.worm.isConfigDiagonal) {
		beadLocator beadIndex;
		beadIndex = path.worm.tail;
		int n = 0;
		do {
			wormBeads(n) = beadIndex;
			beadIndex = path.next(beadIndex);
			++n;
		} while(!all(beadIndex==path.next(path.worm.head)));
	}

	path.printWormConfig(wormBeads);
	path.printLinks<fstream>(communicate()->file("debug")->stream());
	wormBeads.free();
}

