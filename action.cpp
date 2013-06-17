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
 *  Setup the path data members and canonical re-weighting factors.
******************************************************************************/
ActionBase::ActionBase (const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, 
        bool _local, string _name) :
	local(_local), 
    externalPtr(_externalPtr),
    interactionPtr(_interactionPtr),
    name(_name),
    lookup(_lookup),
    path(_path)
{

    /* The default tau scale is 1 */
    shift = 1;

    /* Initialize the separation histogram */
    sepHist.resize(NPCFSEP);
    sepHist = 0;
    cylSepHist.resize(NPCFSEP);
    cylSepHist = 0;
    dSep = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

    /* Needed for canonical ensemble weighting */
    canonical = constants()->canonical();
    numBeads0  = constants()->initialNumParticles()*constants()->numTimeSlices();
    deltaNumBeads2 = 1.0*constants()->deltaNumParticles()*constants()->numTimeSlices();
    deltaNumBeads2 *= deltaNumBeads2;
    window = constants()->window();
    windowWidth = constants()->windowWidth();
    gaussianEnsemble = constants()->gaussianEnsemble();
    gaussianEnsembleSD = constants()->gaussianEnsembleSD();
}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
ActionBase::~ActionBase() {
    sepHist.free();
    cylSepHist.free();
}

/*************************************************************************//**
 *  Update the separation histogram. 
 *
 *  We multiply by the periodic array, as we will only measure separations 
 *  along spatial dimensions with periodic boundary conditions.
******************************************************************************/
inline void ActionBase::updateSepHist(const dVec &_sep) {
	dVec psep;
	psep = path.boxPtr->periodic*_sep;
	int nR = int(sqrt(dot(psep,psep))/dSep);
	if (nR < NPCFSEP)
		++sepHist(nR);
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
   /*if (!canonical)
		return 1.0;
	else {
		int numBeads = path.worm.getNumBeadsOn();
		int numBeadsP = numBeads + deltaNumBeads;
		//double xp = 1.0*(numBeadsP-numBeads0)*(numBeadsP-numBeads0)/deltaNumBeads2;
		//double x = 1.0*(numBeads-numBeads0)*(numBeads-numBeads0)/deltaNumBeads2;
		//return exp(-xp+x)
    }*/
    
    // ensembleWieght returns 1.0 unless window or ensemble weight is used
    if ( window || gaussianEnsemble ){
        int numBeads = path.worm.getNumBeadsOn();
        int numBeadsP = numBeads + deltaNumBeads;
        
        // Check if we are within bead window
        if (window){
            int beadWindowWidth = path.numTimeSlices*windowWidth;
            if ( ( numBeadsP > (numBeads0+beadWindowWidth) )||
                (numBeadsP < (numBeads0-beadWindowWidth)) )
                return 0.0;
        }
        // If we are within window check gaussian weight
        if (gaussianEnsemble){
            double sigmaBeads = path.numTimeSlices*gaussianEnsembleSD;
            double xp = 1.0*(numBeadsP-numBeads0)*(numBeadsP-numBeads0)/
                                (2.0*sigmaBeads*sigmaBeads);
            double x = 1.0*(numBeads-numBeads0)*(numBeads-numBeads0)/
                                (2.0*sigmaBeads*sigmaBeads);
            return exp(-xp+x);
        } 
        else
            return 1.0;
    } 
    else
        return 1.0; // Return 1.0 if no enembleWeight is used
    
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
 *  Return the potential action for a path.
 *
 *  Given a starting and ending bead, compute the potential action for the
 *  path.
 *
 *  NOTE: if startBead == endBead this returns potentialAction(startBead)
 *        I'm not sure if this is the desired behavior yet.
 *
 *  @param startBead the starting bead
 *  @param endBead the end bead on the path
 *  @return the total potential action for a given path
******************************************************************************/
double ActionBase::potentialAction (const beadLocator &startBead, 
        const beadLocator &endBead) {

    double totU = 0.0;
    beadLocator beadIndex;
    beadIndex = startBead;
    do {
        totU += potentialAction(beadIndex);
        beadIndex = path.next(beadIndex);
	} while (!all(beadIndex==path.next(endBead)));

    return totU;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LOCAL ACTION BASE CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Setup the path data members and local action factors.
******************************************************************************/
LocalAction::LocalAction (const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, string _name) :
    ActionBase(_path,_lookup,_externalPtr,_interactionPtr,true,_name) 
{
    /* By default, there is no even/odd difference and the correction
     * coming from gradients of the potential is zero */
    VFactor= 1.0;
    gradVFactor = 0.0;
    eFactor = 0.0;

    setFactor();

}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
LocalAction::~LocalAction() {
}

/**************************************************************************//**
 *  Setup the local potential action and force correction factors.
 *
 *  We pre-compute constant numbers that don't change during the
 *  simulation and are used to compute the action
******************************************************************************/
void LocalAction::setFactor() {
    potentialFactor = constants()->tau()*VFactor;
    fFactor = constants()->tau() * constants()->tau() * constants()->tau() * 
        constants()->lambda() * gradVFactor;
}

/**************************************************************************//**
 *  Return the potential action for all time slices and all particles.
 *
 *  Computes the total potential energy by summing over all particles and time
 *  slices.  We must be very careful with the factor of 1/2 or the
 *  fact that the sum is over i<j for particles to avoid double counting.
 *
 *  !!NB!! It is important to use V and gradVSquared here, as this is used for 
 *  debugging  purposes with check move and must match up with the single bead
 *  calculation.
******************************************************************************/
double LocalAction::potentialAction () {

    double totU = 0.0;

	/* Combine the potential and a possible correction */
	for (int slice = 0; slice < path.numTimeSlices; slice++) {
		eo = (slice % 2);
        totU += VFactor[eo]*tau()*V(slice);
		
		/* We only add the correction if it is finite */
         if ( gradVFactor[eo] > EPS ) 
             totU += gradVFactor[eo] * tau() * tau() * tau() * constants()->lambda() 
                 * gradVSquared(slice);
    }

	return ( totU );
}

/**************************************************************************//**
 *  Return the potential action for a single bead indexed with beadIndex.  
 *
 *  We have to use the beadIndex form of gradVSquared here to ensure that bead
 *  indexed action changes are equal to total potential action changes. We use
 *  a pre-processor directive to determine whether or not to use the
 *  nearest neighbor lookup table.
 *
 *  !!NB!! If we are using checkMove to debug moves you *MUST* always include
 *  the correction term.
******************************************************************************/
double LocalAction::potentialAction (const beadLocator &beadIndex) {

	double bareU = 0.0;
	double corU = 0.0;
	eo = (beadIndex[0] % 2);

 	/* We first compute the base potential action piece */
 	bareU = VFactor[eo]*tau()*Vnn(beadIndex);
 
 	/* If we have a finite correction and are not at the base level
     * of bisection (shift=1) we compute the full correction */
 	 if ( (shift == 1) && (gradVFactor[eo] > EPS) )
 	 	corU = gradVFactor[eo] * tau() * tau() * tau() * constants()->lambda() 
 	 		* gradVnnSquared(beadIndex);
     
     return (bareU + corU);
}

/**************************************************************************//**
 *  Return the bare potential action for a single bead indexed with beadIndex.  
 *
 *  This action corresponds to the primitive approxiamtion and may be used
 *  when attempting updates that use single slice rejection.
******************************************************************************/
double LocalAction::barePotentialAction (const beadLocator &beadIndex) {
    return VFactor[beadIndex[0]%2]*tau()*Vnn(beadIndex);
}

/**************************************************************************//**
 *  Return the potential action correction for a single bead.
 *
 *  Provided that the correction is small, we return its value.
******************************************************************************/
double LocalAction::potentialActionCorrection (const beadLocator &beadIndex) {

	eo = (beadIndex[0] % 2);
    /* If we have a finite correction */
	if (gradVFactor[eo] > EPS) {

        /* We need to update the lookup table here */
		lookup.updateInteractionList(path,beadIndex);

        /* Compute the correction */
       return (gradVFactor[eo] * tau() * tau() * tau() * constants()->lambda() 
               * gradVnnSquared(beadIndex));
    }
    else 
        return 0.0;

}

/**************************************************************************//**
 *  Return the potential action correction for a path.
 *
 *  Given a starting and ending bead, compute the potential action correction
 *  for the path.
 *
 *  @param startBead the starting bead
 *  @param endBead the end bead on the path
******************************************************************************/
double LocalAction::potentialActionCorrection (const beadLocator &startBead, 
        const beadLocator &endBead) {

    double corU = 0.0;
    beadLocator beadIndex;
    beadIndex = startBead;
    do {
        corU += potentialActionCorrection(beadIndex);
        beadIndex = path.next(beadIndex);
	} while (!all(beadIndex==path.next(endBead)));

    return corU;
}

/**************************************************************************//**
 *  The derivative of the potential action wrt tau on all links starting on
 *  slice.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
******************************************************************************/
double LocalAction::derivPotentialActionTau (int slice) {

    double dU = 0.0;
    eo = (slice % 2);

 	/* We first compute the base potential action piece */
 	dU = VFactor[eo]*V(slice);
 
 	/* As long as there is a finite correction, we include it */
 	if ( gradVFactor[eo] > EPS )
 		dU += 3.0 * gradVFactor[eo] * tau() * tau() * constants()->lambda() 
 			* gradVSquared(slice);
    return (dU);

}

/**************************************************************************//**
 *  The derivative of the potential action wrt lambda for all links starting
 *  on slice.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
******************************************************************************/
double LocalAction::derivPotentialActionLambda (int slice) {
    eo = (slice % 2);

 	/* As long as thers is a finite correction, we include it */
 	if ( gradVFactor[eo] > EPS )
 		return gradVFactor[eo] * tau() * tau() * tau() * gradVSquared(slice);
    else
        return 0.0;
}

/**************************************************************************//**
 *  The derivative of the potential action wrt tau on all links starting on
 *  slice within a cutoff radius.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
 *  @param maxR the cutoff radius
******************************************************************************/
double LocalAction::derivPotentialActionTau (int slice, double maxR) {

    double dU = 0.0;
    eo = (slice % 2);

 	/* We first compute the base potential action piece */
 	dU = VFactor[eo]*V(slice,maxR);
 
 	/* As long as there is a finite correction, we include it */
 	if ( gradVFactor[eo] > EPS )
 		dU += 3.0 * gradVFactor[eo] * tau() * tau() * constants()->lambda() 
 			* gradVSquared(slice,maxR);
    return (dU);

}

/**************************************************************************//**
 *  The derivative of the potential action wrt lambda for all links starting
 *  on slice within a cutoff radius.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
 *  @param maxR the cutoff radius
******************************************************************************/
double LocalAction::derivPotentialActionLambda (int slice, double maxR) {
    eo = (slice % 2);

 	/* As long as thers is a finite correction, we include it */
 	if ( gradVFactor[eo] > EPS )
 		return gradVFactor[eo] * tau() * tau() * tau() * gradVSquared(slice,maxR);
    else
        return 0.0;
}

/*************************************************************************//**
 *  Returns the total value of the potential energy, including both the
 *  external potential, and that due to interactions for a single bead.
******************************************************************************/
double LocalAction::V(const beadLocator &bead1) {

	/* We only need to calculate the potential if the bead is on */
	if (path.worm.beadOn(bead1)) {

		bead2[0] = bead1[0];
		
		/* Calculate the total external potential */
		double totVext = externalPtr->V(path(bead1));

		/* Now calculate the total interation potential, neglecting self-interactions */
		double totVint = 0.0;

		/* Get the state of bead 1 */
		beadState state1 = path.worm.getState(bead1);
		
		for (bead2[1]= 0; bead2[1] < path.numBeadsAtSlice(bead1[0]); bead2[1]++) {

			/* Skip self interactions */
			if ( bead2[1] != bead1[1] ) {

				/* get the separation between the two particles */
				sep = path.getSeparation(bead2,bead1);

				/* Now add the interaction potential */
				totVint += path.worm.factor(state1,bead2) * interactionPtr->V(sep);
			} // bead2 != bead1 

		} // for bead2
		return ( totVext + totVint );

	} // if bead1On
	else
		return 0.0;
}

/*************************************************************************//**
 *  Returns the total value of the potential energy, including both the
 *  external potential, and that due to interactions for all particles in
 *  a single time slice.  
 *
 *  This is really only used for either debugging or during the calculation 
 *  of the potential energy. As such, we update the separation histogram here.
******************************************************************************/
double LocalAction::V(const int slice) {

	double totVint = 0.0;
	double totVext = 0.0;

	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

	int numParticles = path.numBeadsAtSlice(slice);

	/* Initialize the separation histogram */
	sepHist = 0;

	/* Calculate the total potential, including external and interaction
	 * effects*/
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

            /* Get the state of bead 1 */
            beadState state1 = path.worm.getState(bead1);

			/* Evaluate the external potential */
			totVext += externalPtr->V(path(bead1));

			/* The loop over all other particles, to find the total interaction
			 * potential */
			for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
				sep = path.getSeparation(bead2,bead1);
				updateSepHist(sep);
				totVint += path.worm.factor(state1,bead2) * interactionPtr->V(sep);
			} // bead2

	} // bead1
	return ( totVext + totVint );
}

/*************************************************************************//**
 *  Returns the total value of the potential energy, including both the
 *  external potential, and that due to interactions for all particles in
 *  a single time slice within a certain cuttoff radius of the center of 
 *  a cylinder.  
 *
 *  This is really only used for either debugging or during the calculation 
 *  of the potential energy. As such, we update the separation histogram here.
******************************************************************************/
double LocalAction::V(const int slice, const double maxR) {

	double totVint = 0.0;
	double totVext = 0.0;
	dVec r1,r2;

	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

	int numParticles = path.numBeadsAtSlice(slice);

	/* Initialize the separation histogram */
	cylSepHist = 0;

	/* Calculate the total potential, including external and interaction
	 * effects*/
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

		r1 = path(bead1);

		if (r1[0]*r1[0] + r1[1]*r1[1] < maxR*maxR) {

			/* Evaluate the external potential */
			totVext += externalPtr->V(path(bead1));

			/* The loop over all other particles, to find the total interaction
			 * potential */
			for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
				r2 = path(bead2);
				sep = path.getSeparation(bead2,bead1);
				if (r2[0]*r2[0] + r2[1]*r2[1] < maxR*maxR) {
					int nR = int(abs(sep[2])/dSep);
					if (nR < NPCFSEP)
						++cylSepHist(nR);
				}
				totVint += interactionPtr->V(sep);
			} // bead2

		} // maxR
	} // bead1

	return ( totVext + totVint );
}

/*************************************************************************//**
 *  Returns the total value of the potential energy, including both the
 *  external potential, and that due to interactions for a single bead with
 *  all interactions confined to a single time slice using a nearest neighbor
 *  grid lookup table.
******************************************************************************/
double LocalAction::Vnn(const beadLocator &bead1) {

	double totVint = 0.0;
	double totVext = 0.0;

	/* We only continue if bead1 is turned on */
	if (path.worm.beadOn(bead1)) {

		/* Fill up th nearest neighbor list */
		lookup.updateInteractionList(path,bead1);

		/* Evaluate the external potential */
		totVext = externalPtr->V(path(bead1));

		/* Get the state of bead 1 */
		beadState state1 = path.worm.getState(bead1);

		/* Sum the interaction potential over all NN beads */
		for (int n = 0; n < lookup.numBeads; n++) {
			totVint += path.worm.factor(state1,lookup.beadList(n)) 
				* interactionPtr->V(lookup.beadSep(n));
		}
	}
	return ( totVext + totVint );
}


/*************************************************************************//**
 *  Returns the total value of the potential energy, including both the
 *  external potential, and that due to interactions for all particles in
 *  a single time slice using a lookup table which only considers particles
 *  within a sphere of some cutoff radius.
******************************************************************************/
double LocalAction::Vnn(const int slice) {

	Array <bool,1> doParticles(path.numBeadsAtSlice(slice));
	doParticles = true;

	double totVint = 0.0;
	double totVext = 0.0;

	iVec gIndex,nngIndex;			// The grid box of a particle
	TinyVector<int,NDIM+1> nnIndex;	// The nearest neighbor boxes of a particle
	TinyVector<int,NDIM+2> hI1,hI2;	// The hash indices

	dVec pos;						// The position of a particle

	beadLocator bead1; 		// Interacting beads
	bead1[0] = slice;

	for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(slice); bead1[1]++) {

		doParticles(bead1[1]) = false;

		/* Accumulate the external potential */
		pos = path(bead1);
		totVext += externalPtr->V(pos);

		/* Get the interaction list */
		lookup.updateInteractionList(path,bead1);

		/* Get the state of bead 1 */
		beadState state1 = path.worm.getState(bead1);

		/* Sum the interaction potential over all NN beads */
		for (int n = 0; n < lookup.numBeads; n++) {
			bead2 = lookup.beadList(n);
			if (doParticles(bead2[1])) {
				sep = path.getSeparation(bead2,bead1);
				totVint += path.worm.factor(state1,bead2) * interactionPtr->V(sep);
			}
		} // n

	} // bead1

	return ( totVext + totVint );
}

/**************************************************************************//**
 *  Return the gradient of the full potential squared for a single bead.  This
 *  includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::gradVSquared(const beadLocator &bead1) {

	double totF2 = 0.0;		// The total force squared

	if (path.worm.beadOn(bead1)) {

		/* All interacting beads are on the same slice. */
		bead2[0] = bead3[0] = bead1[0];

		/* The 'forces' and particle separation */
		dVec Fext1,Fext2;
		dVec Fint1,Fint2,Fint3;

		int numParticles = path.numBeadsAtSlice(bead1[0]);
	
		/* Get the gradient squared part for the external potential*/
		Fext1 = externalPtr->gradV(path(bead1));

		/* We loop through all beads and compute the forces between beads
		 * 1 and 2 */
		Fint1 = 0.0;
		for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

			if (!all(bead1==bead2)) {

				sep = path.getSeparation(bead2,bead1);
				Fint2 = interactionPtr->gradV(sep);
				Fint1 -= Fint2;
				Fext2 = externalPtr->gradV(path(bead2));

				/* There is a single term that depends on this additional interaction
				 * between beads 2 and 3.  This is where all the time is taken */
				Fint3 = 0.0;
				for (bead3[1] = 0; bead3[1] < numParticles; bead3[1]++) {
					if ( !all(bead3==bead2) && !all(bead3==bead1) ) {
						sep = path.getSeparation(bead2,bead3);
						Fint3 += interactionPtr->gradV(sep);
					}
				} // for bead3

				totF2 += dot(Fint2,Fint2) + 2.0*dot(Fext2,Fint2) + 2.0*dot(Fint2,Fint3);

			} //bead1 != bead2

		} // for bead2

		totF2 += dot(Fext1,Fext1) + 2.0 * dot(Fext1,Fint1) + dot(Fint1,Fint1);

	} // bead1 On

	return totF2;
}

/**************************************************************************//**
 *  Return the gradient of the full potential squared for all beads at a single
 *  time slice. 
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::gradVSquared(const int slice) {

	double totF2 = 0.0;

	int numParticles = path.numBeadsAtSlice(slice);

	/* The two interacting particles */
	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

	dVec F;					// The 'force'

	/* We loop over the first bead */
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

		F = 0.0;
		/* Sum up potential for all other active beads in the system */
		for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

			/* Avoid self interactions */
			if (!all(bead1==bead2)) {

				/* The interaction component of the force */
				F += interactionPtr->gradV(path.getSeparation(bead1,bead2));
			} 

		} // end bead2

		/* Now add the external component */
		F += externalPtr->gradV(path(bead1));

		totF2 += dot(F,F);
	} // end bead1

	return totF2;
}

/**************************************************************************//**
 *  Return the gradient of the full potential squared for all beads at a single
 *  time slice restricted inside some cut-off radius.
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::gradVSquared(const int slice, const double maxR) {

	double totF2 = 0.0;

	int numParticles = path.numBeadsAtSlice(slice);

	/* The two interacting particles */
	beadLocator bead1;
	bead1[0] = bead2[0] = slice;
	dVec r1;

	dVec F;					// The 'force'

	/* We loop over the first bead */
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

		r1 = path(bead1);
		if (r1[0]*r1[0] + r1[1]*r1[1] < maxR*maxR) {

			F = 0.0;
			/* Sum up potential for all other active beads in the system */
			for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {
				/* Avoid self interactions */
				if (!all(bead1==bead2)) {

					/* The interaction component of the force */
					F += interactionPtr->gradV(path.getSeparation(bead1,bead2));
				} 
			} // end bead2

			/* Now add the external component */
			F += externalPtr->gradV(path(bead1));

			totF2 += dot(F,F);
		} // maxR
	} // end bead1

	return totF2;
}

/**************************************************************************//**
 *  Return the gradient of the full potential squared for a single bead.  
 *
 *  This includes both the external and interaction potentials using the nearest
 *  neighbor lookup table.
******************************************************************************/
double LocalAction::gradVnnSquared(const beadLocator &bead1) {

	/* This should always be called after Vnn, such that the lookup table 
	 * interaction list has been updated!!! */

	double totF2 = 0.0;		// The total force squared

	if (path.worm.beadOn(bead1)) {

		/* The 'forces' and particle separation */
		dVec Fext1,Fext2;
		dVec Fint1,Fint2,Fint3;

		/* Get the gradient squared part for the external potential*/
		Fext1 = externalPtr->gradV(path(bead1));

		/* We first loop over bead2's interacting with bead1 via the nn lookup table */
		Fint1 = 0.0;
		for (int n = 0; n < lookup.numBeads; n++) {
			bead2 = lookup.beadList(n);

			/* Eliminate self and null interactions */
			if ( !all(bead1 == bead2) ) {

				/* Get the separation between beads 1 and 2 and compute the terms in the
				 * gradient squared */
				sep = lookup.beadSep(n); 
				Fint2 = interactionPtr->gradV(sep);
				Fint1 -= Fint2;
				Fext2 = externalPtr->gradV(path(bead2));

				/* We now loop over bead3, this is the time-intensive part of the calculation */
				Fint3 = 0.0;
				for (int m = 0; m < lookup.numBeads; m++) {
					bead3 = lookup.beadList(m);

					/* Eliminate self-interactions */
					if ( !all(bead3==bead2) && !all(bead3==bead1) ) {
						sep = path.getSeparation(bead2,bead3);
						Fint3 += interactionPtr->gradV(sep);
					}

				} // end bead3

				totF2 += dot(Fint2,Fint2) + 2.0*dot(Fext2,Fint2) + 2.0*dot(Fint2,Fint3);

			} // bead2 != bead1

		} // end bead2

		totF2 += dot(Fext1,Fext1) + 2.0 * dot(Fext1,Fint1) + dot(Fint1,Fint1);

	} // bead1 on

	return totF2;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PRIMITIVE ACTION CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
PrimitiveAction::PrimitiveAction(const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, string _name) :
    LocalAction(_path,_lookup,_externalPtr,_interactionPtr,_name) 
{
	/* Set the action shift */
	shift = 1;

	/* We have no even/odd slice difference here */
	VFactor = 1.0;
	gradVFactor = 0.0;
	eFactor = 0.0;

    setFactor();
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
PrimitiveAction::~PrimitiveAction() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LI-BROUGHTON ACTION CLASS -------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
LiBroughtonAction::LiBroughtonAction(const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, string _name) :
    LocalAction(_path,_lookup,_externalPtr,_interactionPtr,_name) 
{
	/* Set the action shift */
	shift = 1;

	/* We have no even/odd slice difference here */
	VFactor = 1.0;
	gradVFactor = 1.0 / 12.0;
	eFactor = 1.0 / 24.0;

    setFactor();
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
GSFAction::GSFAction(const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, string _name) :
    LocalAction(_path,_lookup,_externalPtr,_interactionPtr,_name) 
{
	/* Set the action shift */
	shift = 1;

	/* Set the value of alpha */ 
	//alpha = 1.0 / 3.0;
	alpha = 0.0;

	/* There are different pre-factors for even/odd slices */
//	VFactor[0] = 4.0/3.0;
//	VFactor[1] = 2.0/3.0;
//
//	gradVFactor[0] = 2.0*(1.0 - alpha)/9.0;
//	gradVFactor[1] = 2.0*alpha/9.0;
//
//	eFactor[0] = (1.0 - alpha)/9.0;
//	eFactor[1] = alpha/9.0;
	
	/* There are different pre-factors for even/odd slices, we use the 
	 * Prokof'ev version here. I have taken alpha = 0 here. */
	VFactor[0] = 2.0/3.0;
	VFactor[1] = 4.0/3.0;

	gradVFactor[0] = 0.0;
	gradVFactor[1] = 2.0/9.0;

	eFactor[0] = 0.0;
	eFactor[1] = 1.0/9.0;

    setFactor();
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
GSFAction::~GSFAction() {}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// NON LOCAL ACTION BASE CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Setup the path data members for non-local actions.
******************************************************************************/
NonLocalAction::NonLocalAction (const Path &_path, LookupTable &_lookup, 
		PotentialBase *_externalPtr, PotentialBase *_interactionPtr, string _name) :
    ActionBase(_path,_lookup,_externalPtr,_interactionPtr,false,_name) 
{
//
}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
NonLocalAction::~NonLocalAction() {
}

/**************************************************************************//**
 *  Return the potential action for all time slices and all particles.
 *
 *  Computes the total potential energy by summing over all particles and time
 *  slices.  We must be very careful with the factor of 1/2 or the
 *  fact that the sum is over i<j for particles to avoid double counting.
******************************************************************************/
double NonLocalAction::potentialAction () {

    double totU = 0.0;
	for (int slice = 0; slice < path.numTimeSlices; slice++) 
        totU += U(slice);

	return ( totU );
}

/**************************************************************************//**
 *  Return the potential action for a single bead indexed with beadIndex.  
 *
 *  @param beadIndex the bead to compute the advanced link action for.
******************************************************************************/
double NonLocalAction::potentialAction (const beadLocator &bead1) {

    /* We only need to calculate the potential action if the two neighboring
     * beads are on */
	if (!path.worm.beadOn(bead1))
        return 0.0;
   
    /* Since the non-local action lives on a link, both bead1 and nextBead1
     * must exist */
    beadLocator nextBead1,nextBead2;
    nextBead1 = path.next(bead1);

    /* Make sure nextBead1 is a real bead and that it is active */
    if ( (all(nextBead1==XXX)) || (!path.worm.beadOn(nextBead1)) )
        return 0.0;

    bead2[0] = bead1[0];

    /* Needed to compute the effective action */
    double lambdaTau = constants()->lambda() * constants()->tau();

    /* Now calculate the total effective interation potential, neglecting self-interactions */
    double totU = 0.0;

    for (bead2[1]= 0; bead2[1] < path.numBeadsAtSlice(bead1[0]); bead2[1]++) {

        /* Skip self interactions */
        if ( bead2[1] != bead1[1] ) {

            /* get the separation between the two particles at the first time
             * slice */
            sep = path.getSeparation(bead1,bead2);

            /* Now we find the separation at the advanced time step */
            nextBead2 = path.next(bead2);

            /* If the imaginary time neighbor exists, compute the effective
             * potential */
            if ( (!all(nextBead2==XXX)) && (path.worm.beadOn(nextBead2)) ) {
                sep2 = path.getSeparation(nextBead1,nextBead2);
                totU += interactionPtr->V(sep,sep2,lambdaTau);
            }
        } // bead2 != bead1 

    } // for bead2

	return ( totU );
}

/**************************************************************************//**
 *  Return the potential action for all time slices and all particles.
 *
 *  Computes the total potential energy by summing over all particles and time
 *  slices.  
******************************************************************************/
double NonLocalAction::U(int slice) {

    double totU = 0.0;

	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

    beadLocator nextBead1,nextBead2;

	int numParticles = path.numBeadsAtSlice(slice);

	/* Initialize the separation histogram */
	sepHist = 0;

	/* Calculate the total potential, including external and interaction
	 * effects*/
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        /* Get the advanced neightbor of bead1 */
        nextBead1 = path.next(bead1);

        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->V(sep,sep2,constants()->lambda()*constants()->tau());
        } // bead2

	} // bead1
	return ( totU );
}

/**************************************************************************//**
 *  The derivative of the potential action wrt tau on all links starting on
 *  slice.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
******************************************************************************/
double NonLocalAction::derivPotentialActionTau (int slice) {

    double totU = 0.0;

	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

    beadLocator nextBead1,nextBead2;

	int numParticles = path.numBeadsAtSlice(slice);

	/* Initialize the separation histogram */
	sepHist = 0;

	/* Calculate the total potential, including external and interaction
	 * effects*/
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        /* Get the advanced neightbor of bead1 */
        nextBead1 = path.next(bead1);

        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);
            updateSepHist(sep);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->dVdtau(sep,sep2,constants()->lambda(),constants()->tau());
        } // bead2

	} // bead1
	return ( totU );
}

/**************************************************************************//**
 *  The derivative of the potential action wrt lambda on all links starting on
 *  slice.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
******************************************************************************/
double NonLocalAction::derivPotentialActionLambda (int slice) {

    double totU = 0.0;

	beadLocator bead1;
	bead1[0] = bead2[0] = slice;

    beadLocator nextBead1,nextBead2;

	int numParticles = path.numBeadsAtSlice(slice);

	/* Initialize the separation histogram */
	sepHist = 0;

	/* Calculate the total potential, including external and interaction
	 * effects*/
	for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        /* Get the advanced neightbor of bead1 */
        nextBead1 = path.next(bead1);

        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->dVdlambda(sep,sep2,constants()->lambda(),constants()->tau());
        } // bead2

	} // bead1
	return ( totU );
}
