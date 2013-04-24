/**
 * @file move.cpp
 * @author Adrian Del Maestro
 *
 * @brief Move class impementaions.
 */

#include "move.h"
#include "path.h"
#include "action.h"
#include "lookuptable.h"
#include "communicator.h"

uint32 MoveBase::totAttempted = 0;
uint32 MoveBase::totAccepted = 0;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// MOVE BASE CLASS -----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
MoveBase::MoveBase (Path &_path, ActionBase *_actionPtr, MTRand &_random, 
        string _name, ensemble _operateOnConfig) :
    name(_name),
    operateOnConfig(_operateOnConfig),
	path(_path), 
	actionPtr(_actionPtr), 
	random(_random) {

 	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;
	success = false;

	sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
MoveBase::~MoveBase() {
	originalPos.free();
}


///@cond DEBUG
/*************************************************************************//**
 *  Print the current worm configuration after a move if DEBUG_WORM is
 *  activated.
******************************************************************************/
inline void MoveBase::printMoveState(string state) {
#ifdef DEBUG_WORM
	/* We make a list of all the beads contained in the worm */
	Array <beadLocator,1> wormBeads;	// Used for debugging
	wormBeads.resize(path.worm.length+1);
	wormBeads = XXX;

	/* Output the worldline configuration */
	communicate()->file("debug")->stream() << "Move State: " << state 
		<< " (" << path.getTrueNumParticles() << ")" << endl;
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
#endif
}

/*************************************************************************//**
 *  We perform a check to make sure that there are no problems with a move,
 *  this is only performed if DEBUG_MOVE is defined.
******************************************************************************/
inline void MoveBase::checkMove(int callNum, double diffA) {
#ifdef DEBUG_MOVE
	/* The first time we call, we calculate the kinetic and potential action */
	if (callNum == 0) {
		oldV = actionPtr->potentialAction();
		oldK = actionPtr->kineticAction();
	}

	/* The second time, we compute the updated values and compare. This would
	 * be used after a keep move */
	if (callNum == 1) {
		newV = actionPtr->potentialAction();
		newK = actionPtr->kineticAction();
		double diffV = newV - oldV;
		if (abs(diffV-diffA) > EPS) {
			communicate()->file("debug")->stream() << format("%-16s%16.6e\t%16.6e\t%16.6e\n") % name 
				% diffV % diffA % (diffV - diffA);
			cout << name << " KEEP " << diffA << endl;
			exit(EXIT_FAILURE);
		}
	}

	/* The third time would be used after an undo move */
	if (callNum == 2) {
		newV = actionPtr->potentialAction();
		newK = actionPtr->kineticAction();
		double diffV = newV - oldV;
		double diffK = newK - oldK;
		if ( (abs(diffV) > EPS) || abs(diffK) > EPS) {
			communicate()->file("debug")->stream() << format("%-16s%16.6e\t%16.6e\n") % name
				% diffV % diffK;
			cout << name << " UNDO " << diffV << " " << diffK << endl;
			exit(EXIT_FAILURE);
		}
	}

	/* This call is used to debug the distributions in the changes in potential
	 * and kinetic energy per particle. diffA is the number of beads touched
	 * in the move here.*/
	if (callNum == -1) {
		newV = actionPtr->potentialAction();
		newK = actionPtr->kineticAction();
		double diffV = newV - oldV;
		communicate()->file("debug")->stream() << format("%-16s%16.6e\t%16.6e\n") % name 
			% ((newK-oldK)/diffA) % ((diffV)/diffA);
	}
#endif
}
///@endcond DEBUG

/*************************************************************************//**
 *  For all move types, if the move is accepted, we simply increment our
 *  accept counters 
******************************************************************************/
void MoveBase::keepMove() {
	numAccepted++;
	totAccepted++;

	/* Restore the shift level for the time step to 1 */
	actionPtr->setShift(1);

	success = true;
}

/*************************************************************************//**
 * Returns a new staging position which will exactly sample the kinetic
 * action. 
 *
 * @param neighborIndex The index of the bead to be updated's neighbor
 * @param endIndex The index of the final bead in the stage
 * @param stageLength The length of the stage
 * @param k The position along the stage
 * @return A NDIM-vector which holds a new random position.
******************************************************************************/
dVec MoveBase::newStagingPosition(const beadLocator &neighborIndex, const beadLocator &endIndex,
		const int stageLength, const int k) {

	PIMC_ASSERT(path.worm.beadOn(neighborIndex));

    /* The rescaled value of lambda used for staging */
    double f1 = 1.0 * (stageLength - k - 1);
    double f2 = 1.0 / (1.0*(stageLength - k));
    double sqrtLambdaKTau = sqrt2LambdaTau * sqrt(f1 * f2);

	/* We find the new 'midpoint' position which exactly samples the kinetic 
	 * density matrix */
	neighborPos = path(neighborIndex);
	newRanPos = path(endIndex) - neighborPos;
	path.boxPtr->putInBC(newRanPos);
	newRanPos *= f2;
	newRanPos += neighborPos;

	/* This is the random kick around that midpoint */
	for (int i = 0; i < NDIM; i++)
		newRanPos[i] = random.randNorm(newRanPos[i],sqrtLambdaKTau);

    path.boxPtr->putInside(newRanPos);

    return newRanPos;
}

/*************************************************************************//**
 * Generates a new position, which exactly samples the free particle
 * density matrix.
 *
 * Compute a position which is selected from a guassian distribution
 * with a mean at a neighboring position, and variance equal to 
 * 2 * lambda * tau.  
 * @param neighborIndex the beadLocator for a neighboring bead
 * @return A randomly generated position which exactly samples 1/2 the
 * kinetic action.
******************************************************************************/
dVec MoveBase::newFreeParticlePosition(const beadLocator &neighborIndex) {

	PIMC_ASSERT(path.worm.beadOn(neighborIndex));

	/* The Gaussian distributed random position */
	for (int i = 0; i < NDIM; i++)
		newRanPos[i] = random.randNorm(path(neighborIndex)[i],sqrt2LambdaTau);

    path.boxPtr->putInside(newRanPos);

    return newRanPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CENTEROFMASS MOVE CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
CenterOfMassMove::CenterOfMassMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name,  ensemble _operateOnConfig) :
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zera */
	numAccepted = numAttempted = numToMove = 0;

	/* We setup the original position array which will store the shift
	 * of a single particle at a time so it can be undone later */
	originalPos.resize(1);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CenterOfMassMove::~CenterOfMassMove() {
}

/*************************************************************************//**
 *  Performs a Center of Mass move.
 * 
 *  For each active particle, we move the center of mass of its wordline
 *  by a fixed random amount.
******************************************************************************/
bool CenterOfMassMove::attemptMove() {

	success = false;

	/* We don't bother performing a center of mass move, unless we have at 
	 * least one bead on slice 0 */
	if (path.numBeadsAtSlice(0) == 0)
		return false;

	/* The initial bead to be moved */
	startBead = 0,random.randInt(path.numBeadsAtSlice(0)-1);

	checkMove(0,0.0);

	/* Increment the number of center of mass moves and the total number
	 * of moves */
	numAttempted++;
	totAttempted++;

	/* Check to see if the start bead is on a worm.  If it is, we start
	 * at the worm tail and end at its head. */
	if (path.worm.foundBead(path,startBead)) {
		startBead = path.worm.tail;
		endBead   = path.next(path.worm.head);
	}
	/* Otherwise, we loop around until we find the initial bead */
	else
		endBead = startBead;

	/* The random shift that will be applied to the center of mass*/
	for (int i = 0; i < NDIM; i++) 
		originalPos(0)[i] = constants()->Delta()*(-0.5 + random.rand());

	/* Here we test to see if any beads are placed outside the box (if we
	 * don't have periodic boundary conditions).  If that is the case,
	 * we don't continue */
	beadLocator beadIndex;
	dVec pos;
	if (int(sum(path.boxPtr->periodic)) != NDIM) {
		beadIndex = startBead;
		do {
			pos = path(beadIndex) + originalPos(0);
			path.boxPtr->putInBC(pos);
			for (int i = 0; i < NDIM; i++) {
				if ((pos[i] < -0.5*path.boxPtr->side[i]) || (pos[i] >= 0.5*path.boxPtr->side[i]))
					return false;
			}
			beadIndex = path.next(beadIndex);
		} while (!all(beadIndex==endBead));
	}

	/* Reset the action */
	newAction = oldAction = 0.0;

	/* We go through the entire worldline, accumulating the old potential action,
	 * shifting the center of mass, then accumulating the new potential action. */ 
	beadIndex = startBead;
	int wlLength = 0;
	do {
		oldAction += actionPtr->potentialAction(beadIndex);
		pos = path(beadIndex) + originalPos(0);
		path.boxPtr->putInBC(pos);
		path.updateBead(beadIndex,pos);

		newAction += actionPtr->potentialAction(beadIndex);
		beadIndex = path.next(beadIndex);
		wlLength++;
	} while (!all(beadIndex==endBead));

	/* The metropolis acceptance step */
	if (random.rand() < min(exp(-(newAction - oldAction)),1.0))  {
		keepMove();
		checkMove(1,newAction-oldAction);
	}
	else {
		undoMove();
		checkMove(2,0.0);
	}

	return success;
}

/*************************************************************************//**
 *  Undo the move by restoring the original particle positions.
******************************************************************************/
void CenterOfMassMove::undoMove() {

	dVec pos;
	beadLocator beadIndex;
	beadIndex = startBead;
	do {
		pos = path(beadIndex) - originalPos(0);
		path.boxPtr->putInBC(pos);
		path.updateBead(beadIndex,pos);
		beadIndex = path.next(beadIndex);
	} while (!all(beadIndex==endBead));

	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// STAGING MOVE CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
StagingMove::StagingMove (Path &_path, ActionBase *_actionPtr, 
        MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zera */
	numAccepted = numAttempted = numToMove = 0;

	/* Resize the original position array */
	originalPos.resize(constants()->Mbar()-1);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
StagingMove::~StagingMove() {
    // empty destructor
}

/*************************************************************************//**
 *  Perform the staging move.
 * 
 *  The staging move is non-local in imaginary time and exactly samples the
 *  kinetic density matrix.
******************************************************************************/
bool StagingMove::attemptMove() {

	success = false;

	/* We don't bother staging if we only have a worm, as this will be taken care
	 * of with other moves.  It is also possible that we could get stuck in an 
	 * infinite loop here.*/
	if (path.getTrueNumParticles()==0)
		return false;

	/* Randomly select the start bead of the stage */
	startBead[0] = random.randInt(path.numTimeSlices-1);
	startBead[1] = random.randInt(path.numBeadsAtSlice(startBead[0])-1);

	/* Now we have to make sure that we are moving an active trajectory, 
	 * otherwise we exit immediatly */
	beadLocator beadIndex;
	beadIndex = startBead;
	for (int k = 0; k < (constants()->Mbar()); k++) {
		if (!path.worm.beadOn(beadIndex) || all(beadIndex==path.worm.head))
			return false;
		beadIndex = path.next(beadIndex);
	}
	endBead = beadIndex;

	/* If we haven't found the worm head, and all beads are on, try to perform the move */

	/* Increment the attempted counters */
	numAttempted++;
	totAttempted++;

	newAction = oldAction = 0.0;

	beadIndex = startBead;
	int k = 0;
	dVec pos;
	do {
		beadIndex = path.next(beadIndex);

		/* Save the old position, and compute the old potential action */
		originalPos(k) = path(beadIndex);
		oldAction += actionPtr->potentialAction(beadIndex);

		/* Generate a new position, and compute the new potential action */
		path.updateBead(beadIndex,
				newStagingPosition(path.prev(beadIndex),endBead,constants()->Mbar(),k));
		newAction += actionPtr->potentialAction(beadIndex);
		++k;
	} while (!all(beadIndex==path.prev(endBead)));

	/* The actual Metropolis test */
	if (random.rand() < min(exp(-(newAction-oldAction)),1.0))
		keepMove();
	else 
		undoMove();
	
	return success;
}

/*************************************************************************//**
 *  Undo the staging move, restoring the original position of each
 *  particle that has been tampered with.
******************************************************************************/
void StagingMove::undoMove() {

	int k = 0;
	beadLocator beadIndex;
	beadIndex = startBead;
	do {
		beadIndex = path.next(beadIndex);
		path.updateBead(beadIndex,originalPos(k));
		++k;
	} while (!all(beadIndex==path.prev(endBead)));

	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// OPEN MOVE CLASS -----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 * Constructor.
******************************************************************************/
OpenMove::OpenMove (Path &_path, ActionBase *_actionPtr, MTRand &_random, 
        string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
OpenMove::~OpenMove() {
}

/*************************************************************************//**
 * Perform an open move.
 *
 * Here we actually attempt to open up the world line configuration with
 * the resulting creation of a worm.  We select a particle and then time 
 * slice at random, and provided the current configuration is diagonal
 * (worm free) we just remove a portion of the particle's worldline.
******************************************************************************/
bool OpenMove::attemptMove() {

	success = false;

	/* We only perform a move if our configuration is diagonal */
	if (path.worm.isConfigDiagonal) {

		/* Get the length of the proposed gap to open up*/
		gapLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int (ceil(log(1.0*gapLength) / log(2.0)-EPS));

		/* Randomly select the head bead, and make sure it is turned on */
		/* THIS IS EXTREMELY IMPORTANT FOR DETAILED BALANCE */
		headBead[0] = random.randInt(path.numTimeSlices-1);
		headBead[1] = random.randInt(path.numBeadsAtSlice(headBead[0])-1);

		/* Find the tail bead */
		tailBead = path.next(headBead,gapLength);

		/* Determine how far apart they are */
		dVec sep;
		sep = path.getSeparation(headBead,tailBead);

		/* We make sure that the proposed worm is not too costly */
		if ( !path.worm.tooCostly(sep,gapLength) ) {

			/* We use the 'true' number of particles here because we are still diagonal, 
			 * so it corresponds to the number of worldlines*/
			double norm = (constants()->C() * constants()->Mbar() * path.worm.getNumBeadsOn())
 					/ actionPtr->rho0(headBead,tailBead,gapLength);

			/* We rescale to take into account different attempt probabilities */
			norm *= constants()->attemptProb("close")/constants()->attemptProb("open");

			/* Weight for ensemble */
			norm *= actionPtr->ensembleWeight(-gapLength+1);

			double muShift = gapLength*constants()->mu()*constants()->tau();

			/* Increment the number of open moves and the total number of moves */
			numAttempted++;
			totAttempted++;
			numAttemptedLevel(numLevels)++;

			/* The temporary head and tail are special beads */
			path.worm.special1 = headBead;
			path.worm.special2 = tailBead;

			/* We now compute the potential energy of the beads that would be
			 * removed if the worm is inserted */
			oldAction = 0.0;
			beadLocator beadIndex;
			beadIndex = headBead;
			do {
				oldAction += actionPtr->potentialAction(beadIndex);
				beadIndex = path.next(beadIndex);
			} while (!all(beadIndex==path.next(tailBead)));

			/* Now perform the metropolis acceptance test based on removing a chunk of
			 * worldline. */
			if (random.rand() < min(norm*exp(oldAction - muShift),1.0)) 
				keepMove();
			else 
				undoMove();

		} // too costly

	} // if config diagonal?

	return success;
}

/*************************************************************************//**
 * We preserve the move, update the acceptance ratios and the location
 * of the worm, and must delete the beads and links in the gap.
******************************************************************************/
void OpenMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;
	
	/* Remove the beads and links from the gap */
	beadLocator beadIndex;
	beadIndex = path.next(headBead);
	while (!all(beadIndex==tailBead)) {
		beadIndex = path.delBeadGetNext(beadIndex);
	} 

	/* Update all the properties of the worm */
	path.worm.update(path,headBead,tailBead);

	/* The configuration with a worm is not diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Opened up a worm.");
	success = true;
}

/*************************************************************************//**
 * Undo the open move, restoring the beads and links.
******************************************************************************/
void OpenMove::undoMove() {

	/* Reset the worm parameters */
	path.worm.reset();
 	path.worm.isConfigDiagonal = true;
	
	printMoveState("Failed to open up a worm.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CLOSE MOVE CLASS ----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
CloseMove::CloseMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CloseMove::~CloseMove() {
	oldBeadOn.free();
}

/*************************************************************************//**
 * Perform a close move.
 * 
 * We attempt to close up the world line configuration if a worm 
 * is already present. This consists of both filling in the beadOn 
 * array as well as generating new positions for the particle beads. 
 * After a sucessful close, we update the number of particles. 
******************************************************************************/
bool CloseMove::attemptMove() {

	success = false;
	/* We first make sure we are in an off-diagonal configuration, and that that
	 * gap is neither too large or too small and that the worm cost is reasonable.
	 * Otherwise, we simply exit the move */
	if ( path.worm.isConfigDiagonal || (path.worm.gap > constants()->Mbar()) 
			|| (path.worm.gap == 0) || path.worm.tooCostly() )
		return false;

	/* Otherwise, proceed with the close move */
	numLevels = int (ceil(log(1.0*path.worm.gap) / log(2.0)-EPS));

	/* Increment the number of close moves and the total number of moves */
	numAttempted++;
	numAttemptedLevel(numLevels)++;
	totAttempted++;

	/* Get the head and new 'tail' slices for the worldline length
	 * to be closed.  These beads are left untouched */
	headBead = path.worm.head;
	tailBead = path.worm.tail;

	/* Compute the part of the acceptance probability that does not
	 * depend on the change in potential energy */
	double norm = actionPtr->rho0(path.worm.head,path.worm.tail,path.worm.gap) /  
		(constants()->C() * constants()->Mbar() * (path.worm.getNumBeadsOn() + path.worm.gap - 1));

	/* We rescale to take into account different attempt probabilities */
	norm *= constants()->attemptProb("open")/constants()->attemptProb("close");

	/* Weight for ensemble */
	norm *= actionPtr->ensembleWeight(path.worm.gap-1);

	/* The change in the number sector */
	double muShift = path.worm.gap*constants()->mu()*constants()->tau();

	/* Compute the potential action for the new trajectory, inclucing a piece coming from the
	 * head and tail. */
	beadLocator beadIndex;
	beadIndex = path.worm.head;
	newAction = actionPtr->potentialAction(beadIndex);
	for (int k = 0; k < (path.worm.gap-1); k++) {
		beadIndex = path.addNextBead(beadIndex,
				newStagingPosition(beadIndex,path.worm.tail,path.worm.gap,k));
		newAction += actionPtr->potentialAction(beadIndex);
	}
	path.next(beadIndex) = path.worm.tail;
	path.prev(path.worm.tail) = beadIndex;
	newAction += actionPtr->potentialAction(path.worm.tail);

	/* Perform the metropolis test */
	if ( random.rand() < min(norm*exp(-newAction + muShift),1.0) ) 
		keepMove();
	else 
		undoMove();

	return success;
}

/*************************************************************************//**
 * We preserve the move, update the acceptance ratios and fully close the
 * worm.
******************************************************************************/
void CloseMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	numAcceptedLevel(numLevels)++;
	totAccepted++;

	/* Update all the properties of the closed worm and our newly
	 * diagonal configuration. */
	path.worm.reset();
 	path.worm.isConfigDiagonal = true;
	
	printMoveState("Closed up a worm.");
	success = true;
}

/*************************************************************************//**
 * We undo the close move, restoring the worm and gap.
******************************************************************************/
void CloseMove::undoMove() {

	/* Delete all the beads that were added. */
	beadLocator beadIndex;
	beadIndex = path.next(path.worm.head);
	while (!all(beadIndex==path.worm.tail))
		beadIndex = path.delBeadGetNext(beadIndex);

	path.next(path.worm.head) = XXX;
	path.prev(path.worm.tail) = XXX;

	/* Make sure we register the off-diagonal configuration */
 	path.worm.isConfigDiagonal = false;
	
	printMoveState("Failed to close up a worm.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// INSERT MOVE CLASS ---------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
InsertMove::InsertMove (Path &_path, ActionBase *_actionPtr, MTRand &_random, 
        string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
InsertMove::~InsertMove() {
}

/*************************************************************************//**
 * Perform an insert move.
 * 
 * Attempt to insert a worm, which acts like a new worldline without
 * periodic boundary conditions.  It is only possible if we are already
 * in a diagonal configuration.  We have to be careful to properly grow
 * our bead and link arrays. The actual number of particles doesn't increase,
 * just the number of active worldlines.
******************************************************************************/
bool InsertMove::attemptMove() {

    /* Get the length of the proposed worm to insert */
    wormLength = 1 + random.randInt(constants()->Mbar()-1);
    numLevels = int(ceil(log(1.0*wormLength) / log(2.0)-EPS));

    /* Increment the number of insert moves and the total number of moves */
    numAttempted++;
    totAttempted++;
    numAttemptedLevel(numLevels)++;

    /* Move normalization factor */
    double norm = constants()->C() * constants()->Mbar() * path.numTimeSlices 
        * path.boxPtr->volume;
    double muShift = wormLength*constants()->tau()*constants()->mu();
    
    /* We rescale to take into account different attempt probabilities */
    norm *= constants()->attemptProb("remove") / constants()->attemptProb("insert");

    /* Weight for ensemble */
    norm *= actionPtr->ensembleWeight(wormLength);

    double actionShift = (log(norm) + muShift)/wormLength;

    /* We pick a random tail slice, and add a new bead */
    int slice = random.randInt(constants()->numTimeSlices()-1);

    /* Choose a random position inside the box for the proposed tail */
    tailBead = path.addBead(slice,path.boxPtr->randPosition(random));
    path.worm.special2 = tailBead;

    double deltaAction = 0.0;
    double PNorm = 1.0;
    double P;

    /* Generate the action for the proposed worm */
    beadLocator beadIndex;
    beadIndex = tailBead;
    deltaAction += actionPtr->potentialAction(beadIndex) - 0.5*actionShift;
    P = min(exp(-deltaAction)/PNorm,1.0);

    /* We perform a metropolis test on the tail bead */
    if ( random.rand() >= P ) {
        undoMove();
        return success;
    }
    PNorm *= P;

    /* Now go through the non head/tail beads */
    for (int k = 1; k < wormLength; k++) {

        beadIndex = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));
        deltaAction += actionPtr->potentialAction(beadIndex) - actionShift;
        P = min(exp(-deltaAction)/PNorm,1.0);

        /* We perform a metropolis test on the single bead */
        if ( random.rand() >= P ) {
            undoMove();
            return success;
        }
        PNorm *= P;
    }
    headBead = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));
    path.worm.special1 = headBead;
    deltaAction += actionPtr->potentialAction(headBead) - 0.5*actionShift;

    /* Perform a final Metropolis test for inserting the full worm*/
    if ( random.rand() < (exp(-deltaAction)/PNorm) )
        keepMove();
    else
        undoMove();

	return success;
}

/*************************************************************************//**
 * Perform an insert move.
 * 
 * Attempt to insert a worm, which acts like a new worldline without
 * periodic boundary conditions.  It is only possible if we are already
 * in a diagonal configuration.  We have to be careful to properly grow
 * our bead and link arrays. The actual number of particles doesn't increase,
 * just the number of active worldlines.
******************************************************************************/
bool InsertMove::attemptMove1() {

	/* We first make sure we are in a diagonal configuration */
	if (path.worm.isConfigDiagonal) {

		/* Get the length of the proposed worm to insert */
		wormLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int(ceil(log(1.0*wormLength) / log(2.0)-EPS));

		/* Increment the number of insert moves and the total number of moves */
		numAttempted++;
		totAttempted++;
		numAttemptedLevel(numLevels)++;

		/* Move normalization factor */
		double norm = constants()->C() * constants()->Mbar() * path.numTimeSlices 
			* path.boxPtr->volume;
		double muShift = wormLength*constants()->tau()*constants()->mu();
		
		/* We rescale to take into account different attempt probabilities */
		norm *= constants()->attemptProb("remove") / constants()->attemptProb("insert");

		/* Weight for ensemble */
		norm *= actionPtr->ensembleWeight(wormLength);

		/* We pick a random tail slice, and add a new bead */
		int slice = random.randInt(constants()->numTimeSlices()-1);

		/* Choose a random position inside the box for the proposed tail */
		tailBead = path.addBead(slice,path.boxPtr->randPosition(random));
		path.worm.special2 = tailBead;

		/* Generate the action for the proposed worm */
		beadLocator beadIndex;
		beadIndex = tailBead;
		newAction = actionPtr->potentialAction(beadIndex);
		for (int k = 1; k < wormLength; k++) {
			beadIndex = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));
			newAction += actionPtr->potentialAction(beadIndex);
		}
		headBead = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));
		path.worm.special1 = headBead;
		newAction += actionPtr->potentialAction(headBead);

		/* Perform the metropolis test */
		if ( random.rand() < min(norm*exp(-newAction + muShift),1.0) ) 
			keepMove();
		else
			undoMove();

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 *  We preserve the insert move, update the acceptance ratios and the worm
 *  location.
******************************************************************************/
void InsertMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Update all the properties of the inserted worm */
	path.worm.update(path,headBead,tailBead);

	/* The configuration with a worm is off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Inserted a worm.");
	success = true;
}

/*************************************************************************//**
 * We undo the insert move, restoring the original worldlines.
******************************************************************************/
void InsertMove::undoMove() {

	/* To undo an insert, we simply remove all the beads that we have
	 * added and reset the worm state.*/
	beadLocator beadIndex;
	beadIndex = tailBead;
	do {
		beadIndex = path.delBeadGetNext(beadIndex);
	} while (!all(beadIndex==XXX));

	path.worm.reset();

	/* Reset the configuration to diagonal */
 	path.worm.isConfigDiagonal = true;

	printMoveState("Failed to insert a worm.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// REMOVE MOVE CLASS ---------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
RemoveMove::RemoveMove (Path &_path, ActionBase *_actionPtr, MTRand &_random, 
        string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RemoveMove::~RemoveMove() {
}

/*************************************************************************//**
 * Perform a remove move.
 * 
 * Attempt to remove a worm thus restoring a diagonal configuration.
 * It is only possible if we are already have an off-diagonal configuration.  
 * We have to be careful to properly shrink our bead and link arrays.  Again
 * the number of true particles doesn't change here, just the number of
 * active worldlines.
******************************************************************************/
bool RemoveMove::attemptMove() {

	/* We first make sure we are in an off-diagonal configuration, and the worm isn't
	 * too short or long, also that we don't remove our last particle */
    if ( (path.worm.length <= constants()->Mbar()) && (path.worm.length >= 1)
			&& (path.getTrueNumParticles() > 1) ) {

		numLevels = int(ceil(log(1.0*path.worm.length) / log(2.0)-EPS));

		/* Increment the number of insert moves and the total number of moves */
		numAttempted++;
		numAttemptedLevel(numLevels)++;
		totAttempted++;

		/* The normalization constant and action shift due to the chemical potential */
		double norm = 1.0 / (constants()->C() * constants()->Mbar() * path.numTimeSlices 
				* path.boxPtr->volume);
		double muShift = path.worm.length * constants()->mu() * constants()->tau();
		
		/* We rescale to take into account different attempt probabilities */
		norm *= constants()->attemptProb("insert") / constants()->attemptProb("remove");

		/* Weight for ensemble */
		norm *= actionPtr->ensembleWeight(-path.worm.length);

        double actionShift = (-log(norm) + muShift)/path.worm.length;

		oldAction = 0.0;
        double deltaAction = 0.0;
        double PNorm = 1.0;
        double P;

        /* First do the head */
		beadLocator beadIndex;
		beadIndex = path.worm.head;
        double factor = 0.5;
		do {
            deltaAction -= actionPtr->potentialAction(beadIndex) - factor*actionShift;
            P = min(exp(-deltaAction)/PNorm,1.0);

            /* We do a single slice Metropolis test and exit the move if we
             * wouldn't remove the single bead */
            if ( random.rand() >= P ) {
                undoMove();
                return success;
            }
            PNorm *= P;

            factor = 1.0; 
			beadIndex = path.prev(beadIndex);
		} while (!all(beadIndex==path.worm.tail));

        /* Add the part from the tail */
        deltaAction -= actionPtr->potentialAction(path.worm.tail) - 0.5*actionShift;

        if ( random.rand() < (exp(-deltaAction)/PNorm) )
            keepMove();
        else
            undoMove();

	} // is the worm length appropriate to remove

	return success;
}

/*************************************************************************//**
 * Perform a remove move.
 * 
 * Attempt to remove a worm thus restoring a diagonal configuration.
 * It is only possible if we are already have an off-diagonal configuration.  
 * We have to be careful to properly shrink our bead and link arrays.  Again
 * the number of true particles doesn't change here, just the number of
 * active worldlines.
******************************************************************************/
bool RemoveMove::attemptMove1() {

	/* We first make sure we are in an off-diagonal configuration, and the worm isn't
	 * too short or long, also that we don't remove our last particle */
	if ( (!path.worm.isConfigDiagonal) 
			&& (path.worm.length <= constants()->Mbar())
			&& (path.worm.length >= 1)
			&& (path.getTrueNumParticles() > 1) ) {

		numLevels = int(ceil(log(1.0*path.worm.length) / log(2.0)-EPS));

		/* Increment the number of insert moves and the total number of moves */
		numAttempted++;
		numAttemptedLevel(numLevels)++;
		totAttempted++;

		/* The normalization constant and action shift due to the chemical potential */
		double norm = 1.0 / (constants()->C() * constants()->Mbar() * path.numTimeSlices 
				* path.boxPtr->volume);
		
		/* We rescale to take into account different attempt probabilities */
		norm *= constants()->attemptProb("insert") / constants()->attemptProb("remove");

		/* Weight for ensemble */
		norm *= actionPtr->ensembleWeight(-path.worm.length);
		
		double muShift = path.worm.length * constants()->mu() * constants()->tau();

		oldAction = 0.0;

		beadLocator beadIndex;
		beadIndex = path.worm.head;
		do {
			oldAction += actionPtr->potentialAction(beadIndex);
			beadIndex = path.prev(beadIndex);
		} while (!all(beadIndex==path.prev(path.worm.tail)));

		/* Perform the metropolis test */
		if ( random.rand() < min(norm*exp(oldAction - muShift),1.0) )
			keepMove();
		else
			undoMove();

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * We preserve the remove move, update the acceptance ratios and shrink
 * the data sets.
******************************************************************************/
void RemoveMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	numAcceptedLevel(numLevels)++;
	totAccepted++;

	printMoveState("About to remove a worm.");

	/* We delete the worm from our data sets and reset all
	 * worm properties.  */
	beadLocator beadIndex;
	beadIndex = path.worm.head;
	do {
		beadIndex = path.delBeadGetPrev(beadIndex);
	} while (!all(beadIndex==XXX));

	path.worm.reset();

	/* The configuration without a worm is diagonal */
 	path.worm.isConfigDiagonal = true;

	printMoveState("Removed a worm.");
	success = true;
}

/*************************************************************************//**
 *  We undo the remove move, but since we haven't actually touched the 
 *  wordlines, we don't have to do anything here.
******************************************************************************/
void RemoveMove::undoMove() {

	/* Reset the configuration to off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Failed to remove a worm.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ADVANCE HEAD MOVE CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
AdvanceHeadMove::AdvanceHeadMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
AdvanceHeadMove::~AdvanceHeadMove() {
}

/*************************************************************************//**
 * Perform an advance head move.
 * 
 * Attempt to advance a worm head in imaginary time by a random number
 * of slices.  We generate the new bead positions from the free particle
 * density matrix.  It is only possible if we already have an off-diagonal 
 * configuration.  
******************************************************************************/
bool AdvanceHeadMove::attemptMove() {

	success = false;
	/* We first make sure we are in an off diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Get the length of the proposed advancement by computing
		 * a random number of levels to be used in the bisection algorithm. */
		advanceLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int (ceil(log(1.0*advanceLength) / log(2.0)-EPS));

		/* Increment the number of advance moves and the total number of moves */
		numAttempted++;
		numAttemptedLevel(numLevels)++;
		totAttempted++;

		double muShift = advanceLength*constants()->tau()*constants()->mu();
		double norm = constants()->attemptProb("recede head") / 
			constants()->attemptProb("advance head");

		/* Weight for ensemble */
		norm *= actionPtr->ensembleWeight(advanceLength);

		/* Make the old head a special bead, and undefine the head */
		path.worm.special1 = path.worm.head;
		path.worm.head = XXX;

		/* Generate the new path, and compute its action, assigning the new head */
		beadLocator beadIndex;
		beadIndex = path.worm.special1;
		newAction = actionPtr->potentialAction(beadIndex);
		for (int k = 0; k < (advanceLength-1); k++) {
			beadIndex = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));
			newAction += actionPtr->potentialAction(beadIndex);
		}
		headBead = path.addNextBead(beadIndex,newFreeParticlePosition(beadIndex));

		/* Assign the new head and compute its action */
		path.worm.head = headBead;
		newAction += actionPtr->potentialAction(headBead);

		/* Perform the metropolis test */
		if ( random.rand() < min(norm*exp(-newAction + muShift),1.0))
			keepMove();
		else 
			undoMove();

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * We preserve the advance head move, update the acceptance ratios and the worm
 * coordinates.
******************************************************************************/
void AdvanceHeadMove::keepMove() {
	
	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Update all the properties of the inserted worm */
	path.worm.update(path,headBead,path.worm.tail);

	/* The configuration with a worm is off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Advanced a worm.");
	success = true;
}

/*************************************************************************//**
 * We undo the advance move, turning off the re-initialized beads and making
 * sure to restore the orginal worm.
******************************************************************************/
void AdvanceHeadMove::undoMove() {

	/* Copy back the head bead */
	path.worm.head = path.worm.special1;

	/* We remove all the beads and links that have been added. */
	beadLocator beadIndex;
	beadIndex = path.next(path.worm.head);
	do {
		beadIndex = path.delBeadGetNext(beadIndex);
	} while (!all(beadIndex==XXX));

	/* Reset the configuration to off-diagonal */
 	path.worm.isConfigDiagonal = false;

	/* Unset the special marker */
	path.worm.special1 = XXX;

	printMoveState("Failed to advance a worm.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ADVANCE TAIL MOVE CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
AdvanceTailMove::AdvanceTailMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
AdvanceTailMove::~AdvanceTailMove() {
}

/*************************************************************************//**
 * Perform an advance tail move.
 * 
 * Here we attempt to advance the tail of a worm in imaginary time by a 
 * random number of slices. This is accomplished by removing beads and the
 * result is a shorter worm.  It is only possible if we already have an 
 * off-diagonal configuration.  
******************************************************************************/
bool AdvanceTailMove::attemptMove() {

	success = false;

	/* We first make sure we are in an off-diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Get the number of time slices we will try to shrink the worm by */
		advanceLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int (ceil(log(1.0*advanceLength) / log(2.0)-EPS));

		/* We now make sure that this is shorter than the length of the worm, 
		 * otherewise we reject the move immediatly */
		if (advanceLength < path.worm.length) {

			/* The proposed new tail */
			tailBead = path.next(path.worm.tail,advanceLength);

			/* The action shift due to the chemical potential */
			double muShift = advanceLength*constants()->tau()*constants()->mu();

			double norm = constants()->attemptProb("recede tail") / 
				constants()->attemptProb("advance tail");

			/* Weight for ensemble */
			norm *= actionPtr->ensembleWeight(-advanceLength);

			/* Increment the number of recede moves and the total number of moves */
			numAttempted++;
			numAttemptedLevel(numLevels)++;
			totAttempted++;

			/* Mark the proposed tail as special */
			path.worm.special1 = tailBead;
			oldAction = 0.0;

			beadLocator beadIndex;
			beadIndex = path.worm.tail;
			do {
				oldAction += actionPtr->potentialAction(beadIndex);
				beadIndex = path.next(beadIndex);
			} while (!all(beadIndex==path.next(tailBead)));

			/* Perform the metropolis test */
			if ( random.rand() < min(norm*exp(oldAction - muShift),1.0) ) 
				keepMove();
			else 
				undoMove();

		} // advanceLength

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * We preserve the advance tail move by deleting the appropriate links and
 * beads and updating the properties of the worm.
******************************************************************************/
void AdvanceTailMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Delete beads and links */
	beadLocator beadIndex;
	beadIndex = path.prev(tailBead);
	do {
		beadIndex = path.delBeadGetPrev(beadIndex);
	} while (!all(beadIndex==XXX));

	/* Update all the changed properties of the inserted worm */
	path.worm.update(path,path.worm.head,tailBead);

	/* The configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Advanced a worm tail.");
	success = true;
}

/*************************************************************************//**
 * We undo the advance tail move, but since we haven't actually touched the 
 * wordlines, we don't have to do anything here.
******************************************************************************/
void AdvanceTailMove::undoMove() {

	/* Make sure the configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	/* Unset the special marker */
	path.worm.special1 = XXX;

	printMoveState("Failed to advance a worm tail.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RECEDE HEAD MOVE CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
RecedeHeadMove::RecedeHeadMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RecedeHeadMove::~RecedeHeadMove() {
}

/*************************************************************************//**
 * Perform a recede head move.
 * 
 * Attempt to propagate a worm backwards in imaginary time by
 * randomly selecting a number of links then attempting to remove them.
 * The number of true particles doesn't change here.
******************************************************************************/
bool RecedeHeadMove::attemptMove() {

	success = false;

	/* We first make sure we are in an off-diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Get the number of time slices we will try to shrink the worm by */
		recedeLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int (ceil(log(1.0*recedeLength) / log(2.0)-EPS));

		/* We now make sure that this is shorter than the length of the worm, 
		 * otherewise we reject the move immediatly */
		if (recedeLength < path.worm.length) {

			/* The proposed new head */
			headBead = path.prev(path.worm.head,recedeLength);

			/* The action shift due to the chemical potential */
			double muShift = recedeLength*constants()->tau()*constants()->mu();

			double norm = constants()->attemptProb("advance head") / 
				constants()->attemptProb("recede head");

			/* Weight for ensemble */
			norm *= actionPtr->ensembleWeight(-recedeLength);

			/* Increment the number of recede moves and the total number of moves */
			numAttempted++;
			numAttemptedLevel(numLevels)++;
			totAttempted++;

			/* Set the proposed head as special */
			path.worm.special1 = headBead;

			oldAction = 0.0;

			beadLocator beadIndex;
			beadIndex = path.worm.head;
			do {
				oldAction += actionPtr->potentialAction(beadIndex);
				beadIndex = path.prev(beadIndex);
			} while (!all(beadIndex==path.prev(headBead)));

			/* Perform the metropolis test */
			if ( random.rand() < min(norm*exp(oldAction - muShift),1.0) ) 
				keepMove();
			else 
				undoMove();

		} // recedeLength

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * We preserve the recede move by deleting the appropriate links and
 * particles and updating the properties of the worm.
******************************************************************************/
void RecedeHeadMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Delete beads and links */
	beadLocator beadIndex;
	beadIndex = path.next(headBead);
	do {
		beadIndex = path.delBeadGetNext(beadIndex);
	} while (!all(beadIndex==XXX));

	/* Update all the changed properties of the inserted worm */
	path.worm.update(path,headBead,path.worm.tail);

	/* The configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Receded a worm head.");
	success = true;
}

/*************************************************************************//**
 * We undo the recede move, but since we haven't actually touched the 
 * wordlines, we don't have to do anything here.
******************************************************************************/
void RecedeHeadMove::undoMove() {

	/* Make sure the configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	/* Unset the special mark */
	path.worm.special1 = XXX;

	printMoveState("Failed to recede a worm head.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RECEDE TAIL MOVE CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
RecedeTailMove::RecedeTailMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RecedeTailMove::~RecedeTailMove() {
}

/*************************************************************************//**
 * Perform a recede tail move.
 * 
 * Attempt to propagate a worm tail backwards in imaginary time by
 * randomly selecting a number of links then attempting to generate new positions
 * which exactly sample the free particle density matrix.
******************************************************************************/
bool RecedeTailMove::attemptMove() {

	success = false;
	/* We first make sure we are in an off diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Get the length of the proposed advancement by computing
		 * a random number of levels to be used in the bisection algorithm. */
		recedeLength = 1 + random.randInt(constants()->Mbar()-1);
		numLevels = int (ceil(log(1.0*recedeLength) / log(2.0)-EPS));

		/* Increment the number of advance moves and the total number of moves */
		numAttempted++;
		numAttemptedLevel(numLevels)++;
		totAttempted++;

		double muShift = recedeLength*constants()->tau()*constants()->mu();
		double norm = constants()->attemptProb("advance tail") / 
			constants()->attemptProb("recede tail");

		/* Weight for ensemble */
		norm *= actionPtr->ensembleWeight(recedeLength);

		/* Make the current tail special, and undefine the tail */
		path.worm.special1 = path.worm.tail;
		path.worm.tail = XXX;

		/* Compute the action for the proposed path, and assign the new tail */
		beadLocator beadIndex;
		beadIndex = path.worm.special1;
		newAction = actionPtr->potentialAction(beadIndex);
		for (int k = 0; k < (recedeLength-1); k++) {
			beadIndex = path.addPrevBead(beadIndex,newFreeParticlePosition(beadIndex));
			newAction += actionPtr->potentialAction(beadIndex);
		}
		tailBead = path.addPrevBead(beadIndex,newFreeParticlePosition(beadIndex));

		/* Assign the new tail bead and compute its action */
		path.worm.tail = tailBead;
		newAction += actionPtr->potentialAction(tailBead);

		/* Perform the metropolis test */
		if ( random.rand() < min(norm*exp(-newAction + muShift),1.0)) 
			keepMove();
		else 
			undoMove();


	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * We preserve the recede tail move, update the acceptance ratios and the worm
 * coordinates.
******************************************************************************/
void RecedeTailMove::keepMove() {
	
	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Update all the properties of the inserted worm */
	path.worm.update(path,path.worm.head,tailBead);

	/* The configuration with a worm is off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Receded a worm tail.");
	success = true;
}

/*************************************************************************//**
 * We undo the recede tail move, turning off the re-initialized beads and making
 * sure to restore the orginal worm.
******************************************************************************/
void RecedeTailMove::undoMove() {

	/* Copy back the tail bead */
	path.worm.tail = path.worm.special1;

	/* We remove all the beads and links that have been added. */
	beadLocator beadIndex;
	beadIndex = path.prev(path.worm.tail);
	do {
		beadIndex = path.delBeadGetPrev(beadIndex);
	} while (!all(beadIndex==XXX));
	path.prev(path.worm.tail) = XXX;

	/* Reset the configuration to off-diagonal */
 	path.worm.isConfigDiagonal = false;

	/* Unset the special mark */
	path.worm.special1 = XXX;

	printMoveState("Failed to recede a worm tail.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SWAP MOVE BASE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
SwapMoveBase::SwapMoveBase (Path &_path, ActionBase *_actionPtr, MTRand &_random,
        string _name, ensemble _operateOnConfig) : 
    MoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
SwapMoveBase::~SwapMoveBase() {
}

/*************************************************************************//**
 * Get the normalization constant for a swap move.
 * 
 * We compute the normalization constant used in both the pivot selection
 * probability as well as the overall acceptance probabilty.  
 * @see Eq. (2.23) of PRE 74, 036701 (2006).
******************************************************************************/
double SwapMoveBase::getNorm(const beadLocator &beadIndex) {

	double Sigma = 0.0;
	double crho0;

	/* We sum up the free particle density matrices for each bead
	 * in the list */
	Sigma = actionPtr->rho0(beadIndex,path.lookup.fullBeadList(0),swapLength);
	cumulant.at(0) = Sigma;
	for (int n = 1; n < path.lookup.fullNumBeads; n++) {
		crho0 = actionPtr->rho0(beadIndex,path.lookup.fullBeadList(n),swapLength);
		Sigma += crho0;
		cumulant.at(n) = cumulant.at(n-1) + crho0;
	}

	/* Normalize the cumulant */
	for (int n = 0; n < path.lookup.fullNumBeads; n++)
        cumulant.at(n) /= Sigma;


	return Sigma;
}

/*************************************************************************//**
 * Select the pivot bead for a swap move.
 * 
 * Here we select a pivot bead from a list with the probability given by 
 * Eq. (2.22) of PRE 74, 036701 (2006).  We use the trick in Ceperly's 
 * lecture notes where we evaluate the cumulative distribution function
 * then generate a uniform random variable and find where it lies in 
 * the ordered CDF list.
******************************************************************************/
beadLocator SwapMoveBase::selectPivotBead() {
	
	/* Generate a uniform deviate, and figure out where it fits in our cumulant
	 * array using a binary search */
	double u = random.rand();
	//Array<double,1>::iterator lb;
	//lb = lowerBound(cumulant.begin(),cumulant.end(),u);

    int index = std::lower_bound(cumulant.begin(),cumulant.end(),u)
            - cumulant.begin();
	/* Copy over the pivot bead and return it */
	beadLocator pivotBead;
	//pivotBead = path.lookup.fullBeadList(lb.position()[0]);
	pivotBead = path.lookup.fullBeadList(index);
	
	return pivotBead;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SWAP HEAD MOVE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
SwapHeadMove::SwapHeadMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    SwapMoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;

	/* Update the sizes of the original position array */
	originalPos.resize(constants()->Mbar()-1);
	originalPos = 0.0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
SwapHeadMove::~SwapHeadMove() {
}

/*************************************************************************//**
 * Perform a swap head move.
 * 
 * Attempt to perform a swap move that samples particle permuations due to the
 * indistinguishability of bosons by reattaching the worm head to another
 * worldline..
******************************************************************************/
bool SwapHeadMove::attemptMove() {

	success = false;

	/* We first make sure we are in an off-diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Initialize */
		pivot = XXX;
		swap = XXX;

		/* Now we figure out how many beads will be involved with the swap bisection. */
		swapLength = constants()->Mbar();
		numLevels = int (ceil(log(1.0*swapLength) / log(2.0)-EPS));

		/* Now form the list of beads which are in the neighbourhood of
		 * the head, but at the advanced time slice */
		int pivotSlice = path.worm.head[0] + swapLength;
		if (pivotSlice >= constants()->numTimeSlices())
			pivotSlice -= constants()->numTimeSlices();

		/* Update the full interaction list */
		path.lookup.updateFullInteractionList(path.worm.head,pivotSlice);

		/* We can only try to make a move if we have at least one bead to swap with */
		if (path.lookup.fullNumBeads > 0) {

			cumulant.resize(path.lookup.fullNumBeads);

			/* We compute the normalization factors using the head bead */
			SigmaHead = getNorm(path.worm.head);

			/* Get the pivot bead */
			pivot = selectPivotBead();

			/* Now we try to find the swap bead.  If we find the worm tail, we immediatly
			 * exit the move */
			beadLocator beadIndex;
			beadIndex = pivot;
			for (int k = 0; k < swapLength; k++) {
				if (all(beadIndex==path.worm.tail))
					return false;
				beadIndex = path.prev(beadIndex);
			}
			swap = beadIndex;

			/* We only continue if the swap is not the tail, and the swap and pivot
			 * grid boxes coincide. */
			if ( !all(path.worm.tail==swap) && path.lookup.gridNeighbors(pivot,swap) ) {

				/* Increment the number of swap moves and the total number of moves */
				numAttempted++;
				totAttempted++;
				numAttemptedLevel(numLevels)++;

				/* Now we create a new list of beads (using the same beadList)
				 * which contains all beads in neighborhood of the swap bead
				 * grid box, but at a time slice advanced by Mbar. We only need to
				 * do this if head and swap are in different grid boxes. */ 
				if (!path.lookup.gridShare(path.worm.head,swap)) {
					path.lookup.updateFullInteractionList(swap,pivotSlice);
					cumulant.resize(path.lookup.fullNumBeads);
				}

				/* Get the normalization factor for the new list */
				SigmaSwap = getNorm(swap);

				/* We now perform a pre-metropolis step on the selected bead. If this 
				 * is not accepted, it is extremely likely that the potential change
				 * will have any effect, so we don't bother with it. */
				if (random.rand() < min(SigmaHead/SigmaSwap,1.0)) {

					/* Mark the special beads */
					path.worm.special1 = swap;
					path.worm.special2 = pivot;

					int k = 0;
					oldAction = newAction = 0.0;

					/* First we compute the old potential action */
					beadIndex = swap;
					do {
						/* Store the original positions */
						if (!all(beadIndex==swap) && !all(beadIndex==pivot)) {
							originalPos(k) = path(beadIndex);
							++k;
						}
						oldAction += actionPtr->potentialAction(beadIndex);
						beadIndex = path.next(beadIndex);
					} while (!all(beadIndex==path.next(pivot)));

					/* Because the staging algorithm requires that we move foward
					 * and backward in imaginary time (for periodic boundary conditions)
					 * we perform the relinking now, and change back to the original 
					 * linking only if we reject the move */

					/* Store temporary bead locators for all the linkages that
					 * will be changed */
					nextSwap = path.next(swap);

					/* Update the links.  Here we exploit the non-constant return reference 
					 * semantics of the next/prev methods */
					path.next(path.worm.head) = nextSwap;
					path.next(swap)           = XXX; 
					path.prev(nextSwap)       = path.worm.head;

					/* Change the former head to a special bead, and assign the new head */
					path.worm.special1 = path.worm.head;
					path.worm.head = swap;

					/* Now we propose a new trajectory and compute its action */
					beadIndex = path.worm.special1;
					k = 0;
					do {
						if (!all(beadIndex==path.worm.special1) && !all(beadIndex==pivot)) {
							path.updateBead(beadIndex,
									newStagingPosition(path.prev(beadIndex),pivot,swapLength,k));
							++k;
						}
						newAction += actionPtr->potentialAction(beadIndex);
						beadIndex = path.next(beadIndex);
					} while (!all(beadIndex==path.next(pivot)));


					if (random.rand() < min(exp(-(newAction - oldAction)),1.0))
						keepMove();
					else
						undoMove();

				} // Pnorm 

			} // if we didn't find the tail and grid boxes coincide

		} // numListBeads = 0

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * Preserve the swap head move.
******************************************************************************/
void SwapHeadMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Now we update the properties of our worm.  Head is replaced with
	 * swap and we must update the length and gap */
	path.worm.update(path,swap,path.worm.tail);

	printMoveState("Performed a head swap.");
	success = true;
}

/*************************************************************************//**
 * If we reject the swap head move, all we need to do is restore the 
 * original positions and undo the new linkages.
******************************************************************************/
void SwapHeadMove::undoMove() {

	/* Reassign the head bead */
	path.worm.head = path.worm.special1;

	/* Return all links to their un-swapped values */
	path.next(path.worm.head) = XXX;
	path.next(swap) = nextSwap;
	path.prev(nextSwap) = swap;

	/* Return the path to its original coordinates */
	beadLocator beadIndex;
	beadIndex = path.next(swap);
	int k = 0;
	do {
		path.updateBead(beadIndex,originalPos(k));
		++k;
		beadIndex = path.next(beadIndex);
	} while (!all(beadIndex==pivot));

	/* Make sure the configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	/* Unset the special beads */
	path.worm.special1 = XXX;
	path.worm.special2 = XXX;

	printMoveState("Failed to perform a head swap.");
	success = false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SWAP TAIL MOVE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
******************************************************************************/
SwapTailMove::SwapTailMove (Path &_path, ActionBase *_actionPtr, 
		MTRand &_random, string _name, ensemble _operateOnConfig) : 
    SwapMoveBase(_path,_actionPtr,_random,_name,_operateOnConfig) {

	/* Initialize private data to zero */
	numAccepted = numAttempted = numToMove = 0;

	/* Initialize the acceptance by level counters */
	numAcceptedLevel.resize(constants()->b()+1);
	numAttemptedLevel.resize(constants()->b()+1);
	numAcceptedLevel  = 0;
	numAttemptedLevel = 0;

	/* Update the sizes of the original position array */
	originalPos.resize(constants()->Mbar()-1);
	originalPos = 0.0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
SwapTailMove::~SwapTailMove() {
}

/*************************************************************************//**
 * Perform a swap tail move.
 * 
 * Try to perform a swap tail move that samples particle permuations due to the
 * indistinguishability of bosons by reattaching the tail to a worldline.
******************************************************************************/
bool SwapTailMove::attemptMove() {

	success = false;

	/* We first make sure we are in an off-diagonal configuration */
	if (!path.worm.isConfigDiagonal) {

		/* Initialize */
		pivot = XXX;
		swap = XXX;

		/* Now we figure out how many beads will be involved with the swap bisection. */
		swapLength = constants()->Mbar();
		numLevels = int (ceil(log(1.0*swapLength) / log(2.0)-EPS));

		/* Now form the list of beads which are in the neighbourhood of
		 * the tail, but at the regressed time slice */
		int pivotSlice = path.worm.tail[0] - swapLength;
		if (pivotSlice < 0)
			pivotSlice += constants()->numTimeSlices();

		path.lookup.updateFullInteractionList(path.worm.tail,pivotSlice);

		/* We can only try to make a move if we have at least one bead to swap with */
		if (path.lookup.fullNumBeads > 0) {

			cumulant.resize(path.lookup.fullNumBeads);

			/* We compute the normalization factors using the tail bead */
			SigmaTail = getNorm(path.worm.tail);

			/* Get the pivot bead */
			pivot = selectPivotBead();

			/* Now we try to find the swap bead.  If we find the worm head, we immediatly
			 * exit the move */
			beadLocator beadIndex;
			beadIndex = pivot;
			for (int k = 0; k < swapLength; k++) {
				if (all(beadIndex==path.worm.head))
					return false;
				beadIndex = path.next(beadIndex);
			}
			swap = beadIndex;

			/* We only continue if we don't find the head, and the pivot and swap grid
			 * boxes coincide, otherwise we reject the move. */
			if ( !all(path.worm.head==swap) && path.lookup.gridNeighbors(pivot,swap) ) {

				/* Increment the number of swap moves and the total number of moves */
				numAttempted++;
				totAttempted++;
				numAttemptedLevel(numLevels)++;

				/* Now we create a new list of beads (using the same beadList)
				 * which contains all beads in neighborhood of the swap bead
				 * grid box, but at a time slice advanced by Mbar. We only 
				 * do this if tail and swap are in different grid boxes. */ 
				if (!path.lookup.gridShare(path.worm.tail,swap)) {
					path.lookup.updateFullInteractionList(swap,pivotSlice);
					cumulant.resize(path.lookup.fullNumBeads);
				}

				/* Get the normalization factor for the new list */
				SigmaSwap = getNorm(swap);

				/* We now perform a pre-metropolis step on the selected bead. If this 
				 * is not accepted, it is extremely unlikely that the potential change
				 * will have any effect, so we don't bother with it. */
				if (random.rand() < min(SigmaTail/SigmaSwap,1.0)) {

					/* Mark the swap and pivot as special */
					path.worm.special1 = swap;
					path.worm.special2 = pivot;

					int k = 0;
					oldAction = newAction = 0.0;

					/* First we compute the old potential action */
//					beadIndex = pivot;
//					do {
//						if (!all(beadIndex==swap) && !all(beadIndex==pivot)) {
//							originalPos(k) = path(beadIndex);
//							++k;
//						}
//						oldAction += actionPtr->potentialAction(beadIndex);
//						beadIndex = path.next(beadIndex);
//					} while (!all(beadIndex==path.next(swap)));

					beadIndex = swap;
					do {
						if (!all(beadIndex==swap) && !all(beadIndex==pivot)) {
							originalPos(k) = path(beadIndex);
							++k;
						}
						oldAction += actionPtr->potentialAction(beadIndex);
						beadIndex = path.prev(beadIndex);
					} while (!all(beadIndex==path.prev(pivot)));

					/* Because the bisection algorithm requires that we move foward
					 * and backward in imaginary time (for periodic boundary conditions)
					 * we perform the relinking now, and change back to the original 
					 * linking only if we reject the move */

					/* Store temporary bead locators for all the linkages that
					 * will be changed */
					prevSwap = path.prev(swap);

					/* Update the links.  Here we exploit the non-constant return reference 
					 * semantics of the next/prev methods */
					path.prev(path.worm.tail) = prevSwap;
					path.prev(swap)           = XXX;
					path.next(prevSwap)       = path.worm.tail;

					/* Change the former tail to a special bead, and assign the new tail */
					path.worm.special1 = path.worm.tail;
					path.worm.tail = swap;

					/* Now we propose a new trajectory and compute its action */
//					beadIndex = pivot;
//					k = 0;
//					do {
//						if (!all(beadIndex==path.worm.special1) && !all(beadIndex==pivot)) {
//							path.updateBead(beadIndex,
//									newStagingPosition(path.prev(beadIndex),path.worm.special1,swapLength,k));
//							++k;
//						}
//						newAction += actionPtr->potentialAction(beadIndex);
//						beadIndex = path.next(beadIndex);
//					} while (!all(beadIndex==path.next(path.worm.special1)));

					k = 0;
					beadIndex = path.worm.special1;
					do {
						if (!all(beadIndex==path.worm.special1) && !all(beadIndex==pivot)) {
							path.updateBead(beadIndex,
									newStagingPosition(path.next(beadIndex),pivot,swapLength,k));
							++k;
						}
						newAction += actionPtr->potentialAction(beadIndex);
						beadIndex = path.prev(beadIndex);
					} while (!all(beadIndex==path.prev(pivot)));

					if (random.rand() < min(exp(-(newAction - oldAction)),1.0))
						keepMove();
					else
						undoMove();

				} // Pnorm 

			} // if we didn't find the tail and the grid boxes coincide

		} // numListBeads = 0

	} // isConfigDiagonal

	return success;
}

/*************************************************************************//**
 * Preserve the swap tail move.
******************************************************************************/
void SwapTailMove::keepMove() {

	/* update the acceptance counters */
	numAccepted++;
	totAccepted++;
	numAcceptedLevel(numLevels)++;

	/* Now we update the properties of our worm.  Tail is replaced with
	 * swap and we must update the length and gap */
	path.worm.update(path,path.worm.head,swap);

	/* Restore the shift level for the time step to 1 */
	actionPtr->setShift(1);

	printMoveState("Performed a tail swap.");
	success = true;
}

/*************************************************************************//**
 * If we reject the swap tail move, all we need to do is restore the 
 * original positions and undo the new linkages.
******************************************************************************/
void SwapTailMove::undoMove() {

	/* Copy back the tail bead */
	path.worm.tail = path.worm.special1;

	/* Return all links to their un-swapped values */
	path.prev(path.worm.tail) = XXX;
	path.prev(swap)           = prevSwap;
	path.next(prevSwap)       = swap;

	beadLocator beadIndex;
//	beadIndex = path.next(pivot);
//	int k = 0;
//	do {
//		path.updateBead(beadIndex,originalPos(k));
//		++k;
//		beadIndex = path.next(beadIndex);
//	} while (!all(beadIndex==swap));

	beadIndex = path.prev(swap);
	int k = 0;
	do {
		path.updateBead(beadIndex,originalPos(k));
		++k;
		beadIndex = path.prev(beadIndex);
	} while (!all(beadIndex==pivot));

	/* Unset the special beads */
	path.worm.special1 = XXX;
	path.worm.special2 = XXX;

	/* Make sure the configuration is still off-diagonal */
 	path.worm.isConfigDiagonal = false;

	printMoveState("Failed to perform a tail swap.");
	success = false;
}
