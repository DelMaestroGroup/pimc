/**
 * @file action.cpp
 * @author Adrian Del Maestro
 *
 * @brief Action implementation.
 */

#include "action.h"
#include "path.h"
#include "potential.h"
#include "lookuptable.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ACTION BASE CLASS --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Setup the path data members and base action factors.
******************************************************************************/
ActionBase::ActionBase (const Path &_path, Potential *_potentialPtr) :
	path(_path), potentialPtr(_potentialPtr) {
		/* The default tau scale is 1 */
		shift = 1;

		/* By default, there is no even/odd difference and the correction
		 * coming from gradients of the potential is zero */
		pFactor = 1.0;
		cFactor = 0.0;
		eFactor = 0.0;

		/* Needed for canonical ensemble weighting */
		canonical = constants()->canonical();
		numBeads0  = constants()->initialNumParticles()*constants()->numTimeSlices();
		deltaNumBeads2 = 1.0*constants()->deltaNumParticles()*constants()->numTimeSlices();
		deltaNumBeads2 *= deltaNumBeads2;
	}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
ActionBase::~ActionBase() {
}

/**************************************************************************//**
 *  The free-particle density matrix for two particles with imaginary
 *  time separation M.
******************************************************************************/
double ActionBase::rho0(const dVec &r1, const dVec &r2, int M) {
	dVec vel;
	vel = r2 - r1;
	path.boxPtr->putInBC(vel);
	double rho0Norm = pow(4.0*M_PI*constants()->lambda()*M*constants()->tau(),-0.5*NDIM);
	return ( rho0Norm * exp(-dot(vel,vel) / (4.0*constants()->lambda()*M*constants()->tau()) ) );
}

/**************************************************************************//**
 *  The free-particle density matrix for two beadLocators with imaginary
 *  time separation M.
******************************************************************************/
double ActionBase::rho0(const beadLocator &bead1, const beadLocator &bead2, int M) {
	dVec vel;
	vel = path.getSeparation(bead1,bead2);
	double rho0Norm = pow(4.0*M_PI*constants()->lambda()*M*constants()->tau(),-0.5*NDIM);
	return ( rho0Norm * exp(-dot(vel,vel) / (4.0*constants()->lambda()*M*constants()->tau()) ) );
}

/**************************************************************************//**
 *  The ensemble weighting factor.
 *
 *  For the grand canonical ensemble we return 1, for the canonical ensemble
 *  we return exp(-[(nb'-nb_0)^2 - (nb-nb_0)^2]/dN^2) where
 *  nb = number of beads
******************************************************************************/
double ActionBase::ensembleWeight(const int deltaNumBeads) {
	if (!canonical)
		return 1.0;
	else {
		int numBeads = path.worm.getNumBeadsOn();
		int numBeadsP = numBeads + deltaNumBeads;
		double xp = 1.0*(numBeadsP-numBeads0)*(numBeadsP-numBeads0)/deltaNumBeads2;
		double x = 1.0*(numBeads-numBeads0)*(numBeads-numBeads0)/deltaNumBeads2;
		return exp(-xp+x);
	}
}

/**************************************************************************//**
 *  Return the kinetic action for a single bead. 
 *
 *  The kinetic action involves two contributions as each bead is linked to 
 *  two other ones at the next and previous time slice.
******************************************************************************/
double ActionBase::kineticAction (const beadLocator &beadIndex) {

	double totK = 0.0;
	/* We get the Kinetic action for a single slice, which connects
	 * to two other slices */
	dVec vel;
	vel = path.getVelocity(path.prev(beadIndex));
	totK += dot(vel,vel);
	vel = path.getVelocity(beadIndex);
	totK += dot(vel,vel);

	return totK / (4.0*constants()->lambda()*tau());
}

/**************************************************************************//**
 *  Return the kinetic action for wlLength time slices starting at the 
 *  bead given by beadIndex.  
 *
 *  The total kinetic action for a string of connected beads all on the same
 *  worldline.
******************************************************************************/
double ActionBase::kineticAction (const beadLocator &beadIndex, int wlLength) {

	double totK = 0.0;

	/* Make a local copy of beadIndex */
	beadLocator beadID;
	beadID = beadIndex;

	/* Starting from the initial bead, move wlLength slices in imaginary time
	 * and accumulate the kinetic action */
	dVec vel;
	for (int n = 0; n < wlLength; n++) {
		vel = path.getVelocity(beadID);
		totK += dot(vel, vel);
		beadID = path.next(beadID);
	} 

	return totK / (4.0*constants()->lambda()*tau());
}

/**************************************************************************//**
 *  Return the total kinetic action.
 *
 *  The total kinetic action evaluated for all particles and all time slices.
******************************************************************************/
double ActionBase::kineticAction () {

	double totK = 0.0;

	/* Calculate the kinetic energy.  Even though there
	 * may be multiple mixing and swaps, it doesn't matter as we always
	 * just advance one time step at a time, as taken care of through the
	 * linking arrays.  This has been checked! */
	beadLocator beadIndex;
	for (int slice = 0; slice < path.numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
			beadIndex = slice,ptcl;
			totK += kineticAction(beadIndex);
		}
	}

	return totK;
}

/**************************************************************************//**
 *  Return the potential action for all time slices and all particles.
 *
 *  Computes the total potential energy by summing over all particles and time
 *  slices.  
******************************************************************************/
double ActionBase::potentialAction () {

	double totV = 0.0;

	/* Combine the potential and a possible correction */
	for (int slice = 0; slice < path.numTimeSlices; slice++) {
		eo = (slice % 2);
		totV += pFactor[eo] * potentialPtr->V(slice);
		
		/* We only add the correction if we are not using the primitive action, and
		 * if the correction pre-factor is not zero */
		if ( (getActionName() != "Base Action") &&  (cFactor[eo] > EPS) )
			totV += pFactor[eo] * cFactor[eo] * tau() * tau() * constants()->lambda() 
				* potentialPtr->gradVSquared(slice);
	}

	return ( tau() * totV );
}

/**************************************************************************//**
 *  Return the potential action for a single bead indexed with beadIndex.  
 *
 *  We have to use the beadIndex form of gradVSquared here to ensure that bead
 *  indexed action changes are equal to total potential action changes. We use
 *  a pre-processor directive to determine whether or not to use the
 *  nearest neighbor lookup table.
******************************************************************************/
double ActionBase::potentialAction (const beadLocator &beadIndex) {

	double totV = 0.0;
	double totF = 0.0;
	eo = (beadIndex[0] % 2);

	/* We use pre-processor directives to determine whether or not we are 
	 * using the NN lookup table for the action */
#ifndef NN_TABLE
	/* We first compute the base potential action piece */
	totV = potentialPtr->V(beadIndex);

	/* If we are not using the base action, and we are at the base level
	 * of any bisection (shift=1) then we must use the full action */
	if ( (shift == 1) && (cFactor[eo] > EPS) )
		totF = cFactor[eo] * tau() * tau() * constants()->lambda() 
			* potentialPtr->gradVSquared(beadIndex);
#else
	/* We first compute the base potential action piece */
	totV = potentialPtr->Vnn(beadIndex);

	/* If we are not using the base action, and we are at the base level
	 * of any bisection (shift=1) then we must use the full action */
	if ( (shift == 1) && (cFactor[eo] > EPS) )
		totF = cFactor[eo] * tau() * tau() * constants()->lambda() 
			* potentialPtr->gradVnnSquared(beadIndex);
#endif

	/* We only include the force term if it is small */
	if (pFactor[eo]*tau()*abs(totF) < 1.0)
		totV += totF;

	return ( pFactor[eo] * tau() * totV );
}

/**************************************************************************//**
 *  Return the potential action for a worldline of length wlLength, 
 *  starting at a bead indexed by beadIndex.
******************************************************************************/
double ActionBase::potentialAction (const beadLocator &beadIndex, 
		int wlLength) {

	double totV = 0.0;

	beadLocator beadID;
	beadID = beadIndex;

	/* Compute the total potential part */
	for (int n = 0; n <= wlLength; n++) {
		eo = (beadID[0] % 2);

		totV += pFactor[eo] * potentialPtr->V(beadID);
		
		/* We only add the correction if we are not using the primitive action */
		if ((getActionName() != "Base Action") && (cFactor[eo] > EPS))
			totV += pFactor[eo] * cFactor[eo] * tau() * tau() * constants()->lambda() 
				* potentialPtr->gradVSquared(beadID);
		
		beadID = path.next(beadID);
	}

	return (tau() * totV);
}

/**************************************************************************//**
 *  Return the total action, for all particles and all time slices.
******************************************************************************/
double ActionBase::totalAction () {

	return ( kineticAction() + potentialAction());
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LI-BROUGHTON ACTION CLASS -------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
LiBroughtonAction::LiBroughtonAction(const Path &_path, 
		Potential *_potentialPtr) : ActionBase(_path,_potentialPtr)
{
	/* Set the action shift */
	shift = 1;

	/* We have no even/odd slice difference here */
	pFactor = 1.0;
	cFactor = 1.0 / 12.0;
	eFactor = 1.0 / 24.0;
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
LiBroughtonAction::~LiBroughtonAction() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GSF ACTION CLASS ----------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
GSFAction::GSFAction(const Path &_path, Potential *_potentialPtr) : 
	ActionBase(_path,_potentialPtr)
{
	/* Set the action shift */
	shift = 1;

	/* Set the value of alpha */ 
	//alpha = 1.0 / 3.0;
	alpha = 0.0;

	/* There are different pre-factors for even/odd slices */
//	pFactor[0] = 4.0/3.0;
//	pFactor[1] = 2.0/3.0;
//
//	cFactor[0] = (1.0 - alpha)/6.0;
//	cFactor[1] = alpha/3.0;
//
//	eFactor[0] = (1.0 - alpha)/9.0;
//	eFactor[1] = alpha/9.0;
	
	/* There are different pre-factors for even/odd slices, we use the 
	 * Prokof'ev version here.*/
	pFactor[0] = 2.0/3.0;
	pFactor[1] = 4.0/3.0;

	cFactor[0] = 0.0;
	cFactor[1] = 1.0/6.0;

	eFactor[0] = 0.0;
	eFactor[1] = 1.0/9.0;
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
GSFAction::~GSFAction() {}

