/**
 * @file pimc.cpp
 * @author Adrian Del Maestro
 * 
 * @brief PathIntegralMonteCarlo class implementation.
 */

#include "pimc.h"
#include "path.h"
#include "action.h"
#include "lookuptable.h"

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
PathIntegralMonteCarlo::PathIntegralMonteCarlo (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, const double maxR, const bool _startWithState) :
	random(_random),
	path(_path),
	centerOfMass(_path,_actionPtr,_random),
	staging(_path,_actionPtr,_random),
	open(_path,_actionPtr,_random),
	close(_path,_actionPtr,_random),
	insert(_path,_actionPtr,_random),
	remove(_path,_actionPtr,_random),
	advanceHead(_path,_actionPtr,_random),
	recedeHead(_path,_actionPtr,_random),
	advanceTail(_path,_actionPtr,_random),
	recedeTail(_path,_actionPtr,_random),
	swapHead(_path,_actionPtr,_random),
	swapTail(_path,_actionPtr,_random),
	diagonalFraction(_path),
	energy(_path,_actionPtr),
	numberParticles(_path),
	particlePositions(_path,0),
	superfluidFraction(_path),
	permutationCycle(_path),
	oneBodyDensityMatrix(_path,_actionPtr,_random),
	pairCorrelation(_path,_actionPtr),
	radialDensity(_path),
	wormProperties(_path),
	numberDistribution(_path),
	cylEnergy(_path,_actionPtr,maxR),
	cylNumberParticles(_path,maxR),
	cylSuperFluidFraction(_path,maxR),
	cylPairCorrelation(_path,_actionPtr,maxR),
	cylOneBodyDensityMatrix(_path,_actionPtr,_random,maxR),
	cylNumberDistribution(_path,maxR),
	cylRadialPotential(_path,_actionPtr,_random,maxR)
{
	/* Are we starting from a saved state? */
	startWithState = _startWithState;

	/* We keep simulating until we have stored a certain number of measurements */
	numStoredBins = 0;

    /* Initialize the number of sweeps and updates per sweep to zero */
    numSteps = 0;
    numUpdates = 0;

	/* Calculate the number of sweeps to make sure we touch every bead */
	numImagTimeSweeps = int(ceil((1.0*constants()->numTimeSlices()/(1.0*constants()->Mbar()))));

	/* These counters are used in the equilibration process */
	numConfig      = 0;
	numDiagonal    = 0;
	numCoM         = 0;
	numCoMAccepted = 0;

	fixedNumLevels = false;

	/* Have we saved a state? */
	savedState  = constants()->restart();

	/* Intialize the config number to zero */
	configNumber = 0;

	/* Fill up the move array */
	move.push_back(&centerOfMass);
	move.push_back(&staging);
	move.push_back(&open);
	move.push_back(&close);
	move.push_back(&advanceHead);
	move.push_back(&recedeHead);
	move.push_back(&advanceTail);
	move.push_back(&recedeTail);
	move.push_back(&swapHead);
	move.push_back(&swapTail);
	move.push_back(&insert);
	move.push_back(&remove);

	/* Now we fill up the cumulative attempt move probability arrays */ 

	/* The cumulative diagonal probabilty vector */
	attemptDiagProb.push_back(constants()->openAttemptProb());
	attemptDiagProb.push_back(attemptDiagProb.at(0) + constants()->insertAttemptProb());
	attemptDiagProb.push_back(attemptDiagProb.at(1) + constants()->stagingAttemptProb());
	attemptDiagProb.push_back(attemptDiagProb.at(2) + constants()->comAttemptProb());

    /* Make sure it adds up to 1 */
	PIMC_ASSERT(attemptDiagProb.back()-1.0 < EPS);
    attemptDiagProb.back() = 1.0 + EPS;

	/* The cumulative off diagonal probability vector */
	attemptOffDiagProb.push_back(constants()->closeAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(0) + constants()->advanceHeadAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(1) + constants()->recedeHeadAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(2) + constants()->advanceTailAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(3) + constants()->recedeTailAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(4) + constants()->removeAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(5) + constants()->swapHeadAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(6) + constants()->swapTailAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(7) + constants()->stagingAttemptProb());
	attemptOffDiagProb.push_back(attemptOffDiagProb.at(8) + constants()->comAttemptProb());

    /* Make sure it adds up to 1 */
	PIMC_ASSERT(attemptOffDiagProb.back()-1.0 < EPS);
    attemptOffDiagProb.back() = 1.0 + EPS;

	/* Add all the estimators. The energy estimator has to be
	 * the first one added. */
	estimator.push_back(&energy);
	estimator.push_back(&numberParticles);
	estimator.push_back(&diagonalFraction);
	estimator.push_back(&particlePositions);
	estimator.push_back(&superfluidFraction);
	estimator.push_back(&permutationCycle);
	estimator.push_back(&pairCorrelation);
	estimator.push_back(&numberDistribution);
	estimator.push_back(&oneBodyDensityMatrix);
	estimator.push_back(&wormProperties);

	/* We only consider the radial density estimator for dimensions > 1 */
	if (NDIM > 1) 
		estimator.push_back(&radialDensity);

	/* We only include the cylinder estimators if we have a tube */
	if (constants()->extPotentialType().find("tube") != string::npos) {
		estimator.push_back(&cylEnergy);
		estimator.push_back(&cylSuperFluidFraction);
		estimator.push_back(&cylNumberParticles);
		estimator.push_back(&cylPairCorrelation);
		estimator.push_back(&cylOneBodyDensityMatrix);
		estimator.push_back(&cylNumberDistribution);
		estimator.push_back(&cylRadialPotential);
	}

	/* If we are restarting, we load the last saved state from disk, otherwise
	 * output the estimator headers. */
	if (constants()->restart())
		loadState();
	else {

		/* If we are starting from an initial configuration, load it. */
		if (startWithState)
			loadState();

		/* Setup all the estimator headers */
		for (vector<EstimatorBase*>::iterator estimatorPtr = estimator.begin();
				estimatorPtr != estimator.end(); ++estimatorPtr)
			(*estimatorPtr)->outputHeader();
	}
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
PathIntegralMonteCarlo::~PathIntegralMonteCarlo () {
}

/**************************************************************************//**
 *  Goes through all possible moves, sampling the partition function.
******************************************************************************/
string PathIntegralMonteCarlo::runMoves(const double x, const int sweep) { 

    success = false;
	string moveName = "NONE";

	/* First we try the moves that can operate on the diagonal ensemble */
	if (path.worm.isConfigDiagonal) {

		if (x < attemptDiagProb.at(0)) {
			success = open.attemptMove();
			moveName = open.name;
		}
		else if (x < attemptDiagProb.at(1)) {
			success = insert.attemptMove();
			moveName = insert.name;
		}
		else if (x < attemptDiagProb.at(2)) {
			success = staging.attemptMove();
			moveName = staging.name;
		}
		else if (x < attemptDiagProb.at(3)) {
			if ((sweep % numParticles)== 0) 
				success = centerOfMass.attemptMove();
			moveName = centerOfMass.name;
		}
		else {
			cout << "Cumulative Diagonal Probability" << endl;
			cout << format("x = %16.14E\n") % x;
		}
	}
	/* Now we try all other moves */
	else{
		if (x < attemptOffDiagProb.at(0)) { 
			success = close.attemptMove();
			moveName = close.name;
		}
		else if (x < attemptOffDiagProb.at(1)) {
			success = advanceHead.attemptMove();
			moveName = advanceHead.name;
		}
		else if (x < attemptOffDiagProb.at(2)) {
			success = recedeHead.attemptMove();
			moveName = recedeHead.name;
		}
		else if (x < attemptOffDiagProb.at(3)) {
			success = advanceTail.attemptMove();
			moveName = advanceTail.name;
		}
		else if (x < attemptOffDiagProb.at(4)) {
			success = recedeTail.attemptMove();
			moveName = recedeTail.name;
		}
		else if (x < attemptOffDiagProb.at(5)) {
			success = remove.attemptMove();
			moveName = remove.name;
		}
		else if (x < attemptOffDiagProb.at(6)) {
			success = swapHead.attemptMove();
			moveName = swapHead.name;
		}
		else if (x < attemptOffDiagProb.at(7)) {
			success = swapTail.attemptMove();
			moveName = swapTail.name;
		}
		else if (x < attemptOffDiagProb.at(8)) {
			success = staging.attemptMove();
			moveName = staging.name;
		}
		else if (x < attemptOffDiagProb.at(9)) {
			if ((sweep % numParticles)== 0) 
				success = centerOfMass.attemptMove();
			moveName = centerOfMass.name;
		}
		else {
			cout << "Cumulative Off-Diagonal Probability!" << endl;
			cout << format("x = %16.14E\n") % x;
		}
	} 
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
void PathIntegralMonteCarlo::equilStep(const uint32 iStep, const bool relaxC0) {

	double x;

	/* How far are we into the equilibration? */
	double equilFrac = (1.0*iStep) / (1.0*constants()->numEqSteps());

    numParticles = path.getTrueNumParticles();

	/* For the first 1/3 of equilibration steps, we only do diagonal moves */
	if (equilFrac < 1.0/3.0 && !startWithState) {

        for (int n = 0; n < numParticles; n++) {
            
			x = random.rand();

            /* Here we do the diagonal pre-equilibration moves, and allow for 
             * optimization of simulation parameters */
            if (x < constants()->comAttemptProb()) {
                centerOfMass.attemptMove();
                numCoM += 1;

                /* We check how many CoM moves we have tried.  Every 200 moves, we see if we need
                 * to adjust Delta, provided we are in the pre-equilibration diagonal state. */
                if ( (numCoM >= 200) 
                        && (constants()->Delta() < 0.5*blitz::min(path.boxPtr->side)) 
                        && (constants()->Delta() > 0.01*blitz::min(path.boxPtr->side)) ) {

                    int numCoMAttempted = numCoM*path.getNumParticles();
                    numCoMAccepted = centerOfMass.getNumAccepted() - numCoMAccepted;
                    double CoMRatio = 1.0*numCoMAccepted / (1.0*numCoMAttempted);
                    if (CoMRatio < 0.2)
                        constants()->shiftDelta(-0.6);
                    else if (CoMRatio < 0.3)
                        constants()->shiftDelta(-0.4);
                    else if (CoMRatio < 0.4)
                        constants()->shiftDelta(-0.2);
                    else if (CoMRatio > 0.6)
                        constants()->shiftDelta(0.2);
                    else if (CoMRatio > 0.7)
                        constants()->shiftDelta(0.4);
                    else if (CoMRatio > 0.8)
                        constants()->shiftDelta(0.6);

                    /* Reset the accepted counter */
                    numCoMAccepted = centerOfMass.getNumAccepted();
                    numCoM = 0;
                }
            } // CoM Move
            else {
                /* Attemp the staging Move */
                for (int sweep = 0; sweep < numImagTimeSweeps; sweep++) 
                    staging.attemptMove();

            } //staging move
        } // numParticles

	} // equilFrac < 1/3 and !startWithState

    /* Start the 2/3 off diagonal portion of the equilibration */
	else {

        /* How many updates should we perform? */
        numUpdates = int(ceil(path.worm.getNumBeadsOn()/(1.0*constants()->Mbar())));

        /* Increment the number of off-diagonal steps */
        numSteps++;

        for (int n = 0; n < numUpdates; n++) {

            /* Generate random number and run through all moves */
            x = random.rand();
            runMoves(x,n);

            /* We accumulate the number of diagonal configurations */
            numConfig++;
            if (path.worm.isConfigDiagonal) 
                ++numDiagonal;

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
    if ((iStep > 0) && (iStep % 200) == 0)
        saveState();
}

/**************************************************************************//**
 *  Equilibration.
 * 
 *  The equilibration method, where we perform fully diagonal moves 1/3 of the
 *  time, finding an optimal, COM step, then for another 1/3 we perform 
 *  all moves adjusting C0 to yield a diagonal fraction near 0.75 and and for 
 *  the final 1/3 we just equilibrate.
******************************************************************************/
void PathIntegralMonteCarlo::equilStepOffDiagonal(const uint32 iStep, const bool relaxC0) {

	double x;

	/* How far are we into the equilibration? */
	double equilFrac = (1.0*iStep) / (1.0*constants()->numEqSteps());

    /* Increment the number of steps */
    numSteps++;
    
    /* How many updates should we perform? */
    numUpdates = int(ceil(path.worm.getNumBeadsOn()/(1.0*constants()->Mbar())));
    numParticles = path.getTrueNumParticles();

    for (int n = 0; n < numUpdates; n++) {
		
        /* Generate random number and run through all moves */
        x = random.rand();
        runMoves(x,n);

        /* We accumulate the number of diagonal configurations */
        numConfig++;
        if (path.worm.isConfigDiagonal) 
            ++numDiagonal;

        /* Every 200 steps, we check the diagonal fraction and update the
         * worm constant for the first 1/3 of equilibration */
        if ( (relaxC0) && (numSteps > 200) && (equilFrac < 1.0/3.0) ) {

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

    } // sweeps loop	

    /* Save a state to disk every 200 equilibrium steps */
    if ((iStep > 0) && (iStep % 200) == 0)
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

	uint32 binSize = 100;		// The number of MC steps in an output bin
	string moveName;

    numUpdates = int(ceil(path.worm.getNumBeadsOn()/(1.0*constants()->Mbar())));
    numParticles = path.getTrueNumParticles();

	/* We run through all moves, making sure that we could have touched each bead at least once */
    for (int n = 0; n < numUpdates ; n++) {

        double x = random.rand();
        moveName = runMoves(x,n);

        /* We measure the one body density matrix here to ensure adequate statistics */
        //oneBodyDensityMatrix.sample();

        //if (constants()->extPotentialType().find("tube") != string::npos)
        //		cylOneBodyDensityMatrix.sample();

	} // n

	/* Perform all measurements */
	for (vector<EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr)
		(*estimatorPtr)->sample();

	/* Every binSize measurements, we output averages to disk and record the 
	 * state of the simulation on disk.  */
	if (estimator.front()->getNumAccumulated() >= binSize) {
		for (vector<EstimatorBase*>::iterator estimatorPtr = estimator.begin();
				estimatorPtr != estimator.end(); ++estimatorPtr) {
			if ((*estimatorPtr)->getNumAccumulated() >= binSize)
				(*estimatorPtr)->output();
		}
		saveState();
		++numStoredBins;
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
		% move.front()->getTotAcceptanceRatio();
	communicate()->file("log")->stream() << endl;

	/* Ouptut all the move acceptance information to disk */
	string moveName;
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		moveName = (*movePtr)->name;

		/* We only output levels for those moves which have a variable size */
		if ( (moveName != staging.name) && (moveName != centerOfMass.name)
				&& (moveName != swapHead.name) && (moveName != swapTail.name) ) {
			for (int n = 0; n <= constants()->b(); n++) {
				communicate()->file("log")->stream() << format("%-12s Level %-10d\t:\t%7.5f\t(%d/%d)\n") 
					% moveName % n % (*movePtr)->getAcceptanceRatioLevel(n) 
					% (*movePtr)->numAcceptedLevel(n) % (*movePtr)->numAttemptedLevel(n);
			}
		}
		communicate()->file("log")->stream() << format("%-29s\t:\t%7.5f\t(%d/%d)\n") % moveName
			% (*movePtr)->getAcceptanceRatio() % (*movePtr)->numAccepted
			% (*movePtr)->numAttempted;
		communicate()->file("log")->stream() << endl;
	}
	communicate()->file("log")->stream() << "---------- End Acceptance Data -----------------" << endl;

	communicate()->file("log")->stream() << endl;
	communicate()->file("log")->stream() << endl;

	/* Output the estimator statistics to the log file */
	communicate()->file("log")->stream() << "---------- Begin Estimator Data ----------------" << endl;
	communicate()->file("log")->stream() << endl;
	for (vector <EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->file("log")->stream() << format("%-29s\t:\t%16d\t%16d\n") % (*estimatorPtr)->getName()
			% (*estimatorPtr)->getNumSampled() % (*estimatorPtr)->getTotNumAccumulated();

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

	/* Prepare the state file for writing */
	communicate()->file("state")->reset();

	/* We First write the current total number of world lines */
	communicate()->file("state")->stream() << path.getNumParticles() << endl;

	/* Now write the total acceptance information for all moves */
	communicate()->file("state")->stream() << format("%16d\t%16d\n") 
		% move.front()->totAccepted % move.front()->totAttempted;

	/* Now record the individual move acceptance information,
	 * first for the diagonal, then off-diagonal*/
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		communicate()->file("state")->stream() << format("%16d\t%16d\n") 
			% (*movePtr)->numAccepted % (*movePtr)->numAttempted;
	}

	/* Output the estimator sampling information */
	for (vector <EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->file("state")->stream() << format("%16d\t%16d\n") 
			% (*estimatorPtr)->getTotNumAccumulated() 
			% (*estimatorPtr)->getNumSampled();
	}

	/* Now we output the actual path and worldline data */
	communicate()->file("state")->stream() << path.beads << endl;
	communicate()->file("state")->stream() << path.nextLink << endl;
	communicate()->file("state")->stream() << path.prevLink << endl;

	/* Output the worm data */
	communicate()->file("state")->stream() << path.worm.beads << endl;

	/* Save the state of the random number generator */
	uint32 randomState[random.SAVE];
	random.save(randomState);
	for (int i = 0; i < random.SAVE; i++) 
		communicate()->file("state")->stream() << randomState[i] << " ";
	communicate()->file("state")->stream() << endl;

	/* Rename and copy the file. */
	communicate()->file("state")->rename();
//	communicate()->file("state")->stream().close();
}

/**************************************************************************//**
 *  Load a classical ground state from file.
******************************************************************************/
void PathIntegralMonteCarlo::loadClassicalState(Array <dVec,2> &tempBeads,
        Array <unsigned int, 2> &tempWormBeads, int numWorldLines) {

    /* We go through each active worldline and create a new classical
     * configuration */
    int ptcl = 0;
    beadLocator beadIndex;
    for (int n = 0; n < numWorldLines; n++) {
        beadIndex = 0,n;

        /* Check if the bead is on */
        if (tempWormBeads(beadIndex)) {
            
            /* Assign the classical configuration */
            path.beads(Range::all(),ptcl) = tempBeads(beadIndex);
            path.worm.beads(Range::all(),ptcl) = 1;
            ptcl++;
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
    doBead = true;

    beadLocator startBead,beadIndex;
    beadLocator newStartBead;
    newStartBead = 0,0;
    int ptcl = 0;
    int slice = 0;

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

                path.beads(slice % numTimeSlices,ptcl) = tempBeads(beadIndex);
                path.worm.beads(slice % numTimeSlices,ptcl) = 1;

                beadIndex = tempNextLink(beadIndex);
                ++slice;

                /* Do a forward reconnection, provided we are not at the
                 * last bead */
                if ( ((slice % numTimeSlices) == 0) && !all(beadIndex==startBead)) {
                    path.nextLink(numTimeSlices-1,ptcl) = 0,ptcl+1;
                    path.prevLink(0,ptcl+1) = numTimeSlices-1,ptcl;
                    ++ptcl;
                }
            } while (!all(beadIndex==startBead));

            /* Now we have to add the remaining beads and perform the final
             * reconnection */
            for (int tslice = (slice % numTimeSlices); tslice < numTimeSlices; tslice++) {
                path.beads(tslice,ptcl) = tempBeads(beadIndex);
                path.worm.beads(tslice,ptcl) = 1;
            }
            path.nextLink(numTimeSlices-1,ptcl) = newStartBead;
            path.prevLink(newStartBead) = numTimeSlices-1,ptcl;
            ++ptcl;

        } // doBead
    } // n

    /* Free local memory */
    doBead.free();
}

/**************************************************************************//**
 *  Load a previous state of the simulation from disk.
******************************************************************************/
void PathIntegralMonteCarlo::loadState() {

	uint32 temp;

	/* We first read the former total number of world lines */
	int numWorldLines;
	int numTimeSlices = path.numTimeSlices;

	communicate()->file("init")->stream() >> numWorldLines;

	/* Reset the total acceptance information */
	communicate()->file("init")->stream() >> temp >> temp;
	move.front()->resetTotAccept();

	/* Reset all the individual move acceptance information */
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		communicate()->file("init")->stream() >> temp >> temp;
		(*movePtr)->resetAccept();
	}

	/* Reset previous estimator sampling information */
	for (vector<EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->file("init")->stream() >> temp >> temp;
		(*estimatorPtr)->restart(0,0);
	}

	/* Now we resize all path data members and read them from the init state file */
	path.beads.resize(numTimeSlices,numWorldLines);
	path.nextLink.resize(numTimeSlices,numWorldLines);
	path.prevLink.resize(numTimeSlices,numWorldLines);
    path.worm.beads.resize(numTimeSlices,numWorldLines);

    /* A temporary container for the beads array */
    Array <dVec,2> tempBeads;

	/* Get the worldline configuration */
	communicate()->file("init")->stream() >> tempBeads;

    /* The temporary number of time slices */
    int tempNumTimeSlices = tempBeads.rows();

    if (tempNumTimeSlices == numTimeSlices) {
    
        /* Copy over the beads array */
        path.beads = tempBeads;

        /* Get the link arrays */
        communicate()->file("init")->stream() >> path.nextLink;
        communicate()->file("init")->stream() >> path.prevLink;

        /* Repeat for the worm file */
        communicate()->file("init")->stream() >> path.worm.beads;

    } // locBeads.rows() == numTimeSlices
    else {

        /* Initialize the links */
        firstIndex i1;
        secondIndex i2;
        path.prevLink[0] = i1-1;
        path.prevLink[1] = i2;
        path.nextLink[0] = i1+1;
        path.nextLink[1] = i2;
	
        /* Here we implement the initial periodic boundary conditions in 
         * imaginary time */
        path.prevLink(0,Range::all())[0] = numTimeSlices-1;
        path.nextLink(numTimeSlices-1,Range::all())[0] = 0;

        /* Reset the worm.beads array */
        path.worm.beads = 0;

        /* Temporary containers for the links and worm beads */
        Array <beadLocator,2> tempNextLink;
        Array <beadLocator,2> tempPrevLink;
        Array <unsigned int,2> tempWormBeads;

        /* Get the link arrays and worm file */
        communicate()->file("init")->stream() >> tempNextLink;
        communicate()->file("init")->stream() >> tempPrevLink;
        communicate()->file("init")->stream() >> tempWormBeads;

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
                if (!path.worm.beads(beadIndex)) {
                    path.nextLink(beadIndex) = XXX;
                    path.prevLink(beadIndex) = XXX;
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
            communicate()->file("init")->stream() >> randomState[i];
        random.load(randomState);
    }

    /* Reset the number of on beads */
    path.worm.resetNumBeadsOn();

    /* Left pack the input data array */
    path.leftPack();

	/* Go through all beads, and make sure they fit inside the simulation cell.
	 * At the same time, determine how many active beads there are per slice */
	beadLocator beadIndex;
	path.numBeadsAtSlice = 0;
	for (beadIndex[0] = 0; beadIndex[0] < numTimeSlices; ++beadIndex[0]) {
		for (beadIndex[1] = 0; beadIndex[1] < numWorldLines; ++beadIndex[1]) {
			path.boxPtr->putInside(path(beadIndex));
			if (path.worm.beadOn(beadIndex))
				++path.numBeadsAtSlice(beadIndex[0]);
		} // particles
	} // slices

	/* Reset the worm parameters */
	path.worm.isConfigDiagonal = true;
	path.worm.reset();

	/* Resize and update the lookup table arrays */
	path.lookup.resizeList(numWorldLines);
	path.lookup.updateGrid(path);

	/* Close the file */
	communicate()->file("init")->close();

    /* Free up memory */
    tempBeads.free();
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

	int numParticles = path.getNumParticles();
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
