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
	diagonalFraction(_path,1),
	energy(_path,_actionPtr,1),
	numberParticles(_path,1),
	superfluidFraction(_path,1),
	permutationCycle(_path,1),
	oneBodyDensityMatrix(_path,_actionPtr,_random,15),
	pairCorrelation(_path,_actionPtr,1),
	radialDensity(_path,1),
	wormProperties(_path,5),
	numberDistribution(_path,1),
	cylEnergy(_path,_actionPtr,maxR,1),
	cylNumberParticles(_path,maxR,1),
	cylSuperFluidFraction(_path,maxR,1),
	cylPairCorrelation(_path,_actionPtr,maxR,1),
	cylOneBodyDensityMatrix(_path,_actionPtr,_random,maxR,15),
	cylNumberDistribution(_path,maxR,1),
	cylRadialPotential(_path,_actionPtr,_random,maxR,1)
{
	/* Are we starting from a saved state? */
	startWithState = _startWithState;

	/* We keep simulating until we have stored a certain number of measurements */
	numStoredBins = 0;

	/* Calculate the number of sweeps to make sure we touch every bead */
	numSweeps = int(ceil((1.0*constants()->numTimeSlices()/(1.0*constants()->Mbar()))));

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
	PIMC_ASSERT(attemptDiagProb.back()-1.0 < EPS);

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

	PIMC_ASSERT(attemptOffDiagProb.back()-1.0 < EPS);

	/* Add all the estimators. The energy estimator has to be
	 * the first one added. */
	estimator.push_back(&energy);
	estimator.push_back(&numberParticles);
	estimator.push_back(&diagonalFraction);
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

	bool success = false;		// Was a move successful?
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
			if (sweep == 0) 
				success = centerOfMass.attemptMove();
			moveName = centerOfMass.name;
		}
		else {
			cout << "Problem with Cumulative Diagonal Probability!" << endl;
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
			if (sweep == 0)
				success = centerOfMass.attemptMove();
			moveName = centerOfMass.name;
		}
		else {
			cout << "Problem with Cumulative Off-Diagonal Probability!" << endl;
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

	/* For the first 1/3 of equilibration steps, we only do diagonal moves */
	if (equilFrac < 1.0/3.0 && !startWithState) {

		/* We loop over all particles */
		for (int n = 0; n < constants()->initialNumParticles(); n++) {

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
				for (int sweep = 0; sweep < numSweeps; sweep++) 
						staging.attemptMove();

			} //staging move

		} // initialNumParticles

	} // equilFrac < 1/3 and !startWithState
	else {

		/* Start the 2/3 off diagonal portion of the equilibration */
		
		/* Save a diagonal configuration to disk */
		if (!savedState) 
			saveState();

		for (int sweep = 0; sweep < numSweeps; sweep++) {

			for (int n = 0; n < constants()->initialNumParticles(); n++) {

				/* Generate random number and run through all moves */
				x = random.rand();
				runMoves(x,sweep);

				/* We accumulate the number of diagonal configurations */
				numConfig++;
				if (path.worm.isConfigDiagonal) 
					++numDiagonal;

				/* Every 200*checkNum steps, we check the diagonal fraction and update the
				 * worm constant for the second 1/3 of equilibration */
				int checkNum = constants()->initialNumParticles()*numSweeps;
				if ( (constants()->C0() > 0.0) && (numConfig > 200*checkNum) && 
						(equilFrac < 2.0/3.0) ) {

					double diagFrac = 1.0*numDiagonal / (1.0*numConfig);

					cout << format("%4.2f\t%8.5f\t") % diagFrac % constants()->C0();

					if (relaxC0) {
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

					} // relaxC0

					cout << format("%8.5f\t%5d\t%8.6f\n") % constants()->C0() 
						% path.getTrueNumParticles() 
						% (1.0*path.getTrueNumParticles()/path.boxPtr->volume);

					/* Reset the counters */
					numDiagonal = 0;
					numConfig = 0;

				}

			} // particle loop

		} // sweeps loop

	} // 2/3 off-diagonal equilibration 
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

	uint32 binSize = 200;		// The number of MC steps in an output bin
	string moveName;

	/* We run through all moves, making sure that we could have touched each bead at least once */
	for (int sweep = 0; sweep < numSweeps; sweep++) {

		for (int n = 0; n < constants()->initialNumParticles(); n++) {

			double x = random.rand();
			moveName = runMoves(x,sweep);

			/* We measure the one body density matrix here to ensure adequate statistics */
			oneBodyDensityMatrix.sample();

			if (constants()->extPotentialType().find("tube") != string::npos)
				cylOneBodyDensityMatrix.sample();

		} // n

	} // sweep

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
	communicate()->logFile() << endl;
	communicate()->logFile() << endl;

	communicate()->logFile() << "---------- Begin Acceptance Data ---------------" << endl;
	communicate()->logFile() << endl;
	communicate()->logFile() << format("%-29s\t:\t%7.5f\n") % "Total Rate" 
		% move.front()->getTotAcceptanceRatio();
	communicate()->logFile() << endl;

	/* Ouptut all the move acceptance information to disk */
	string moveName;
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		moveName = (*movePtr)->name;

		/* We only output levels for those moves which have a variable size */
		if ( (moveName != staging.name) && (moveName != centerOfMass.name)
				&& (moveName != swapHead.name) && (moveName != swapTail.name) ) {
			for (int n = 0; n <= constants()->b(); n++) {
				communicate()->logFile() << format("%-12s Level %-10d\t:\t%7.5f\t(%d/%d)\n") 
					% moveName % n % (*movePtr)->getAcceptanceRatioLevel(n) 
					% (*movePtr)->numAcceptedLevel(n) % (*movePtr)->numAttemptedLevel(n);
			}
		}
		communicate()->logFile() << format("%-29s\t:\t%7.5f\t(%d/%d)\n") % moveName
			% (*movePtr)->getAcceptanceRatio() % (*movePtr)->numAccepted
			% (*movePtr)->numAttempted;
		communicate()->logFile() << endl;
	}
	communicate()->logFile() << "---------- End Acceptance Data -----------------" << endl;

	communicate()->logFile() << endl;
	communicate()->logFile() << endl;

	/* Output the estimator statistics to the log file */
	communicate()->logFile() << "---------- Begin Estimator Data ----------------" << endl;
	communicate()->logFile() << endl;
	for (vector <EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->logFile() << format("%-29s\t:\t%16d\t%16d\n") % (*estimatorPtr)->name
			% (*estimatorPtr)->getNumSampled() % (*estimatorPtr)->getTotNumAccumulated();

	}
	communicate()->logFile() << endl;
	communicate()->logFile() << "---------- End Estimator Data ------------------" << endl;

}

/**************************************************************************//**
 *  Save the state of the simulation to disk, including all path, worm,
 *  move and estimator data.
******************************************************************************/
void PathIntegralMonteCarlo::saveState() {

	savedState = true;

	/* Prepare the state file for writing */
	communicate()->resetStateFile(ios::out|ios::trunc);

	/* We First write the current total number of world lines */
	communicate()->stateFile() << path.getNumParticles() << endl;

	/* Now write the total acceptance information for all moves */
	communicate()->stateFile() << format("%16d\t%16d\n") 
		% move.front()->totAccepted % move.front()->totAttempted;

	/* Now record the individual move acceptance information,
	 * first for the diagonal, then off-diagonal*/
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		communicate()->stateFile() << format("%16d\t%16d\n") 
			% (*movePtr)->numAccepted % (*movePtr)->numAttempted;
	}

	/* Output the estimator sampling information */
	for (vector <EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->stateFile() << format("%16d\t%16d\n") 
			% (*estimatorPtr)->getTotNumAccumulated() 
			% (*estimatorPtr)->getNumSampled();
	}

	/* Now we output the actual path and worldline data */
	communicate()->stateFile() << path.beads << endl;
	communicate()->stateFile() << path.nextLink << endl;
	communicate()->stateFile() << path.prevLink << endl;

	/* Output the worm data */
	communicate()->stateFile() << path.worm.beads << endl;

	/* Save the state of the random number generator */
	uint32 randomState[random.SAVE];
	random.save(randomState);
	for (int i = 0; i < random.SAVE; i++) 
		communicate()->stateFile() << randomState[i] << " ";
	communicate()->stateFile() << endl;

	/* Close the file */
	communicate()->stateFile().close();
}

/**************************************************************************//**
 *  Load a previous state of the simulation from disk.
******************************************************************************/
void PathIntegralMonteCarlo::loadState() {

	uint32 temp;

	/* We first read the former total number of world lines */
	int numWorldLines;
	int numTimeSlices = path.numTimeSlices;

	communicate()->initFile() >> numWorldLines;

	/* Reset the total acceptance information */
	communicate()->initFile() >> temp >> temp;
	move.front()->resetTotAccept();

	/* Reset all the individual move acceptance information */
	for (vector<MoveBase*>::iterator movePtr = move.begin(); movePtr != move.end(); ++movePtr) {
		communicate()->initFile() >> temp >> temp;
		(*movePtr)->resetAccept();
	}

	/* Reset previous estimator sampling information */
	for (vector<EstimatorBase*>::iterator estimatorPtr = estimator.begin();
			estimatorPtr != estimator.end(); ++estimatorPtr) {
		communicate()->initFile() >> temp >> temp;
		(*estimatorPtr)->restart(0,0);
	}

	/* Now we resize all path data members and read them from the init state file */
	path.beads.resize(numTimeSlices,numWorldLines);
	path.nextLink.resize(numTimeSlices,numWorldLines);
	path.prevLink.resize(numTimeSlices,numWorldLines);

	/* Get the worldline configuration */
	communicate()->initFile() >> path.beads;

	/* Get the link arrays */
	communicate()->initFile() >> path.nextLink;
	communicate()->initFile() >> path.prevLink;

	/* Repeat for the worm file */
	path.worm.beads.resize(numTimeSlices,numWorldLines);
	communicate()->initFile() >> path.worm.beads;

	/* Load the state of the random number generator, only if we are restarting 
	 * the simulation */
	if (constants()->restart()) {
		uint32 randomState[random.SAVE];
		for (int i = 0; i < random.SAVE; i++) 
			communicate()->initFile() >> randomState[i];
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
	communicate()->initFile().close();
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
	communicate()->wlFile() << format("REMARK [CONFIG %04d]\n") % configNumber;

	/* Output the unit cell information.  It is always cubic.  Everything is scaled by
	 * an overall factor for better visualization. */

	double scale = 10.0;				
	int i;
	communicate()->wlFile() << format("%-6s") % "CRYST1";
	for (i = 0; i < NDIM; i++) 
		communicate()->wlFile() << format("%9.3f") % (scale*path.boxPtr->side[i]);
	while (i < 3) {
		communicate()->wlFile() << format("%9.3f") % 1.0;
		i++;
	}
	communicate()->wlFile() << format("%7.2f%7.2f%7.2f %-11s%4d\n") % 90.0 % 90.0 % 90.0 % "P 1" % 1;

	/* We output the atom block */
	beadLocator beadIndex;
	for (int n = 0; n < numWorldLines;  n++) {
		beadIndex = startBead(n);
		do {
			/* We give the zero-time-slice bead a special name */
			if (beadIndex[0] == 0) {
				communicate()->wlFile() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
					% beadNum(beadIndex) % "H0" % wlLength(n) % " ";
			}
			else {
				communicate()->wlFile() << format("%-6s%5d %-4s %03d %9s") % "ATOM" 
					% beadNum(beadIndex) % "HE" % wlLength(n) % " ";
			}

			/* Output the coordinates in 3D */
			for (i = 0; i < NDIM; i++) {
				communicate()->wlFile() << format("%8.3f") % (scale*path(beadIndex)[i]);
			}
			while (i < 3) {
				communicate()->wlFile() << format("%8.3f") % 0.0;
				i++;
			}
			communicate()->wlFile() << format("%14s\n") % "HE";

			beadIndex = path.next(beadIndex);
		} while (!all(beadIndex==endBead(n)));
	}
	communicate()->wlFile() <<("TER\n");

	/* Now output the connect block */
	for (int n = 0; n < numWorldLines;  n++) {
		beadIndex = startBead(n);
		do {
			communicate()->wlFile() << format("%-6s%5d") % "CONECT" % beadNum(beadIndex);
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
					communicate()->wlFile() << format("%5d") % beadNum(prevIndex);
			}

			/* Now the next bead */
			if (path.worm.beadOn(nextIndex)) {
				dVec sep;
				sep = path(nextIndex) - path(beadIndex);
				if (dot(sep,sep) < path.boxPtr->rcut2)
					communicate()->wlFile() << format("%5d") % beadNum(nextIndex);
			}
			communicate()->wlFile() << endl;

			beadIndex = path.next(beadIndex);
		} while (!all(beadIndex==endBead(n)));
	}
	communicate()->wlFile() <<("END\n");

	/* Free up memory */
	startBead.free();
	endBead.free();
	wlLength.free();
	beadNum.free();
	doBead.free();
}
