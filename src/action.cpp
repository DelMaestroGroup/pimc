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
#include "wavefunction.h"

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
        WaveFunctionBase *_waveFunctionPtr, bool _local, std::string _name, 
        double _endFactor, int _period) :
    local(_local), 
    period(_period),
    externalPtr(_externalPtr),
    interactionPtr(_interactionPtr),
    name(_name),
    lookup(_lookup),
    path(_path),
    waveFunctionPtr(_waveFunctionPtr),
    endFactor(_endFactor)
{

    /* The default tau scale is 1 */
    shift = 1;

    /* Initialize the separation histogram */
    sepHist.resize(NPCFSEP);
    sepHist = 0;
    cylSepHist.resize(NPCFSEP);
    cylSepHist = 0;
    dSep = 0.5*sqrt(NDIM)*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);
    dPerSep = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

    /* Needed for canonical ensemble weighting */
    canonical = constants()->canonical();
    numBeads0  = constants()->initialNumParticles()*constants()->numTimeSlices();
    window = constants()->window();
    windowWidth = constants()->windowWidth();
    gaussianEnsemble = constants()->gaussianEnsemble();
    gaussianEnsembleSD = constants()->gaussianEnsembleSD();
}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
ActionBase::~ActionBase() {
}

/*************************************************************************//**
 *  Update the separation histogram. 
 *
 *  We multiply by the periodic array, as we will only measure separations 
 *  along spatial dimensions with periodic boundary conditions.
******************************************************************************/
inline void ActionBase::updateSepHist(const dVec &_sep) {
    dVec psep;
    /* Temporarly computing full radial g(r) */
    /* psep = path.boxPtr->periodic*_sep; */
    psep = _sep;
    int nR = int(sqrt(dot(psep,psep))/dSep);
    if (nR >= 0 && nR < NPCFSEP)
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
double ActionBase::rho0(const beadLocator &bead1, const beadLocator &bead4, int M) {
    dVec vel;
    vel = path.getSeparation(bead1,bead4);
    double rho0Norm = pow(4.0*M_PI*constants()->lambda()*M*constants()->tau(),-0.5*NDIM);
    return ( rho0Norm * exp(-dot(vel,vel) / (4.0*constants()->lambda()*M*constants()->tau()) ) );
}

/**************************************************************************//**
 *  The free-particle density matrix for a given spatial and temporal
 *  separation.
 *
 *  @param vel The spatial separtion R - R'
 *  @param M The imaginary time separation
 *  @return The free density matrix
******************************************************************************/
double ActionBase::rho0(const dVec &vel, const int M) {
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
 *  The total kinetic action for a std::string of connected beads all on the same
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
    /* Calculate the kinetic energy.  Even though there
     * may be multiple mixing and swaps, it doesn't matter as we always
     * just advance one time step at a time, as taken care of through the
     * linking arrays.  This has been checked! */

    double totK = 0.0;
    beadLocator beadIndex;
    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        beadIndex[0] = slice;
        int numBeads = path.numBeadsAtSlice(slice);
        for (int ptcl = 0; ptcl < numBeads; ptcl++) {
            beadIndex[1] = ptcl;
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
    } while (!(beadIndex == path.next(endBead)));
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
        PotentialBase *_externalPtr, PotentialBase *_interactionPtr, 
        WaveFunctionBase *_waveFunctionPtr, const std::array<double,2> &_VFactor, 
        const std::array<double,2> & _gradVFactor, bool _local, std::string _name,
        double _endFactor, int _period) :
    ActionBase(_path,_lookup,_externalPtr,_interactionPtr,_waveFunctionPtr,
            _local,_name,_endFactor,_period), 
    VFactor(_VFactor),
    gradVFactor(_gradVFactor)
{
    shift = 1;
}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
LocalAction::~LocalAction() {
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
        totU += VFactor[eo]*tau()*sum(V(slice));
        
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
 *
 *  @param beadIndex the bead that we are computing the potential action for
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
     
#if PIGS
     /* We tack on a trial wave function and boundary piece if necessary */  
     if ( (beadIndex[0] == 0) || (beadIndex[0] == (constants()->numTimeSlices()-1)) ) {
         bareU *= 0.5*endFactor;
         bareU -= log(waveFunctionPtr->PsiTrial(beadIndex[0]));
     }
#endif

     return (bareU + corU);
}

/**************************************************************************//**
 *  Return the bare potential action for a single bead indexed with beadIndex.  
 *
 *  This action corresponds to the primitive approximation and may be used
 *  when attempting updates that use single slice rejection.
******************************************************************************/
double LocalAction::barePotentialAction (const beadLocator &beadIndex) {

    double bareU = VFactor[beadIndex[0]%2]*tau()*Vnn(beadIndex);

#if PIGS
     /* We tack on a trial wave function and boundary piece if necessary */  
     if ( (beadIndex[0] == 0) || (beadIndex[0]== (constants()->numTimeSlices()-1)) ) {
         bareU *= 0.5*endFactor;
         bareU -= log(waveFunctionPtr->PsiTrial(beadIndex[0]));
     }
#endif

    return bareU;
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
    } while (!(beadIndex == path.next(endBead)));

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
    dU = VFactor[eo]*sum(V(slice));
 
    /* As long as there is a finite correction, we include it */
    if ( gradVFactor[eo] > EPS )
        dU += 3.0 * gradVFactor[eo] * tau() * tau() * constants()->lambda() 
            * gradVSquared(slice);
    return (dU);

}

/**************************************************************************//**
 *  The second derivative of the potential action wrt tau on all links 
 *  starting on slice.
 *
 *  It is essential to have these slice overloaded values as we need to be
 *  careful of the factor of 1/2 or i<j in the full potential action.
 *
 *  @param slice the imaginary timeslice of the first link
******************************************************************************/
double LocalAction::secondderivPotentialActionTau (int slice) {

    double d2U = 0.0;
    eo = (slice % 2);

    /* As long as there is a finite correction, we include it */
    if ( gradVFactor[eo] > EPS )
        d2U += 6.0 * gradVFactor[eo] * tau() * constants()->lambda() 
            * gradVSquared(slice);
    return (d2U);

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

    /* As long as there is a finite correction, we include it */
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

        /* Get the state of bead 1 */
        beadState state1 = path.worm.getState(bead1);
        
        /* Calculate the total external potential */
        double totVext = path.worm.factor(state1)*externalPtr->V(path(bead1));

        /* Now calculate the total interation potential, neglecting self-interactions */
        double totVint = 0.0;

        
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
 *  Returns the external and interaction potential energy, for all particles on
 *  a single time slice.  
 *
 *  This is really only used for either debugging or during the calculation 
 *  of the potential energy. As such, we update the separation histogram here.
******************************************************************************/
std::array<double,2> LocalAction::V(const int slice) {

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
            totVext += path.worm.factor(state1)*externalPtr->V(path(bead1));

            /* The loop over all other particles, to find the total interaction
             * potential */
            for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
                sep = path.getSeparation(bead2,bead1);
                updateSepHist(sep);
                totVint += path.worm.factor(state1,bead2) * interactionPtr->V(sep);
            } // bead2

    } // bead1

    /* Separate the external and interaction parts */ 
    return {totVext,totVint};
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

    double r1sq,r2sq;

    beadLocator bead1;
    bead1[0] = bead2[0] = slice;

    int numParticles = path.numBeadsAtSlice(slice);

    /* Initialize the separation histogram */
    cylSepHist = 0;

    /* Calculate the total potential, including external and interaction
     * effects*/
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        r1 = path(bead1);

        r1sq = 0.0;
        for (int i = 0; i < NDIM-1; i++)
            r1sq += r1[i]*r1[i];

        if (r1sq < maxR*maxR) {

            /* Evaluate the external potential */
            totVext += externalPtr->V(path(bead1));

            /* The loop over all other particles, to find the total interaction
             * potential */
            for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
                r2 = path(bead2);

                r2sq = 0.0;
                for (int i = 0; i < NDIM-1; i++)
                    r2sq += r2[i]*r2[i];

                sep = path.getSeparation(bead2,bead1);
                if (r2sq < maxR*maxR) {
                    int nR = int(abs(sep[2])/dPerSep);
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

        /* Get the state of bead 1 */
        beadState state1 = path.worm.getState(bead1);

        /* Evaluate the external potential */
        totVext = path.worm.factor(state1)*externalPtr->V(path(bead1));

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

    DynamicArray <bool,1> doParticles(path.numBeadsAtSlice(slice));
    doParticles = true;

    double totVint = 0.0;
    double totVext = 0.0;

    //FIXME these variables are unused?
    //iVec gIndex,nngIndex;           // The grid box of a particle
    //std::array<int,NDIM+1> nnIndex; // The nearest neighbor boxes of a particle
    //std::array<int,NDIM+2> hI1,hI2; // The hash indices

    dVec pos;                       // The position of a particle

    beadLocator bead1;      // Interacting beads
    bead1[0] = slice;

    for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(slice); bead1[1]++) {

        doParticles(bead1[1]) = false;

        /* Get the state of bead 1 */
        beadState state1 = path.worm.getState(bead1);

        /* Accumulate the external potential */
        pos = path(bead1);
        totVext += path.worm.factor(state1)*externalPtr->V(pos);

        /* Get the interaction list */
        lookup.updateInteractionList(path,bead1);

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

    double totF2 = 0.0;     // The total force squared

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
        Fint1 = dVec{};
        for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

            if (!(bead1 == bead2)) {

                sep = path.getSeparation(bead2,bead1);
                Fint2 = interactionPtr->gradV(sep);
                Fint1 -= Fint2;
                Fext2 = externalPtr->gradV(path(bead2));

                /* There is a single term that depends on this additional interaction
                 * between beads 2 and 3.  This is where all the time is taken */
                Fint3 = dVec{};
                for (bead3[1] = 0; bead3[1] < numParticles; bead3[1]++) {
                    if ( !(bead3 == bead2) && !(bead3 == bead1) ) {
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

    dVec F;                 // The 'force'

    /* We loop over the first bead */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        F = dVec{};
        /* Sum up potential for all other active beads in the system */
        for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

            /* Avoid self interactions */
            if (!(bead1 == bead2)) {

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
    double r1sq;

    dVec F;                 // The 'force'

    /* We loop over the first bead */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        r1 = path(bead1);
        r1sq = 0.0;
        for (int i = 0; i < NDIM-1; i++)
            r1sq += r1[i]*r1[i];

        if (r1sq < maxR*maxR) {

            F = dVec{};
            /* Sum up potential for all other active beads in the system */
            for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {
                /* Avoid self interactions */
                if (!(bead1 == bead2)) {

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

    double totF2 = 0.0;     // The total force squared

    if (path.worm.beadOn(bead1)) {

        /* The 'forces' and particle separation */
        dVec Fext1,Fext2;
        dVec Fint1,Fint2,Fint3;

        /* Get the gradient squared part for the external potential*/
        Fext1 = externalPtr->gradV(path(bead1));

        /* We first loop over bead2's interacting with bead1 via the nn lookup table */
        Fint1 = dVec{};
        for (int n = 0; n < lookup.numBeads; n++) {
            bead2 = lookup.beadList(n);

            /* Eliminate self and null interactions */
            if ( !(bead1 == bead2) ) {

                /* Get the separation between beads 1 and 2 and compute the terms in the
                 * gradient squared */
                sep = lookup.beadSep(n); 
                Fint2 = interactionPtr->gradV(sep);
                Fint1 -= Fint2;
                Fext2 = externalPtr->gradV(path(bead2));

                /* We now loop over bead3, this is the time-intensive part of the calculation */
                Fint3 = dVec{};
                for (int m = 0; m < lookup.numBeads; m++) {
                    bead3 = lookup.beadList(m);

                    /* Eliminate self-interactions */
                    if ( !(bead3 == bead2) && !(bead3 == bead1) ) {
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

/**************************************************************************//**
 *  Return the gradient of the full potential for all beads at a single
 *  time slice. 
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
dVec LocalAction::gradientV(const int slice) {

    int numParticles = path.numBeadsAtSlice(slice);

    /* The two interacting particles */
    beadLocator bead1;
    bead1[0] = bead2[0] = slice;

    dVec gV;                    // The 'force'

    /* We loop over the first bead */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        gV = dVec{};
        /* Sum up potential for all other active beads in the system */
        for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

            /* Avoid self interactions */
            if (!(bead1 == bead2)) {

                /* The interaction component of the force */
                gV += interactionPtr->gradV(path.getSeparation(bead1,bead2));
            } 
        } // end bead2

        /* Now add the external component */
        gV += externalPtr->gradV(path(bead1));
    } // end bead1

    return gV;
}

/*************************************************************************//**
 *  Returns the T-matrix needed to compute gradU.
 *
 *  This is used for the virial energy estimator in the TI or GSF action.
 *  CURRENTLY ONLY RETURNS VALUES FOR INTERACTIONS, NOT EXTERNAL.
******************************************************************************/
dMat LocalAction::tMatrix(const int slice) {

    int numParticles = path.numBeadsAtSlice(slice);

    dMat tMat{}; // tMat(row, col)

    dVec rDiff{};
    double rmag = 0.0;
    double d2V = 0.0;
    double dV = 0.0;
    
    /* two interacting particles */
    beadLocator bead1;
    bead1[0] = bead2[0] = slice;
    
    for (int a=0; a<NDIM; a++){
        for (int b=0; b<NDIM; b++){

            /* loop over first bead */
            for (bead1[1]=0; bead1[1]<numParticles; bead1[1]++){
                
                /* loop over all other beads */
                for (bead2[1]=0; bead2[1]<numParticles; bead2[1]++){

                    /* avoid self-interactions */
                    if (!(bead1 == bead2)) {
                        
                        rDiff = path.getSeparation(bead1, bead2);
                        rmag = sqrt(dot(rDiff,rDiff));
                        d2V = interactionPtr->grad2V(rDiff);
                        d2V += externalPtr->grad2V(path(bead1));
                        dV = sqrt(dot(interactionPtr->gradV(rDiff)
                                    ,interactionPtr->gradV(rDiff)));
                        dV += sqrt(dot(externalPtr->gradV(path(bead1))
                                    ,externalPtr->gradV(path(bead1))));

                        tMat(a,b) += rDiff(a)*rDiff(b)*d2V/(rmag*rmag);
                        if (a != b)
                            tMat(a,b) -= (rDiff(a)*rDiff(b)/pow(rmag,3))*dV;
                        else
                            tMat(a,b) += (1.0/rmag - rDiff(a)*rDiff(b)/pow(rmag,3))*dV;
                    }
                } // end bead2
            } // end bead1
        }
    }

    tMat /= 2.0; // don't want to double count interactions
    
    return ( tMat ); 
}

/**************************************************************************//**
 *  Return the sum over particles at a given time slice of the 
 *  gradient of the action potential for a single bead dotted into the
 *  bead's position.  Used for regular virial energy and thermodynamic
 *  pressure.
 *  
 *  NOTE:  This is the first term.  These were split up because it is
 *      beneficial to be able to return them separately when computing
 *      the specific heat via the centroid virial estimator.
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::rDOTgradUterm1(const int slice) {

    eo = slice % 2;
    int numParticles = path.numBeadsAtSlice(slice);

    /* The two interacting particles */
    beadLocator bead1;
    bead1[0] = bead2[0] = slice;

    dVec gVi, gV;
    double rDotgV = 0.0;

    /* We loop over the first bead */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        gVi = dVec{};
        /* Sum potential of bead1 interacting with all other beads at
         * a given time slice.*/
        for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

            /* Avoid self interactions */
            if (!(bead1 == bead2)) {

                /* The interaction component of the force */
                gVi += interactionPtr->gradV(path.getSeparation(bead1,bead2));
            } 
        } // end bead2
        
        /* Now add the external component */
        gV = (externalPtr->gradV(path(bead1)) + gVi);

        rDotgV += dot(gV, path(bead1));

    } // end bead1

    return (VFactor[eo]*constants()->tau()*rDotgV);
}

/**************************************************************************//**
 *  Return the sum over particles at a given time slice of the 
 *  gradient of the action potential for a single bead dotted into the
 *  bead's position.  Used for regular virial energy and thermodynamic
 *  pressure.
 *  
 *  NOTE:  This is the second term
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::rDOTgradUterm2(const int slice) {

    eo = slice % 2;
    int numParticles = path.numBeadsAtSlice(slice);

    /* The two interacting particles */
    beadLocator bead1;
    bead1[0] = bead2[0] = slice;

    dVec gVi, gVe, gV, g2V;

    double term2 = 0.0;
    if (gradVFactor[eo] > EPS){

        /* constants for tMatrix */
        dVec rDiff{};
        double rmag = 0.0;
        double d2V = 0.0;
        double dV = 0.0;
        double dVe = 0.0;
        double dVi = 0.0;
        double g2Vi = 0.0;
        double g2Ve = 0.0;

        /* We loop over the first bead */
        for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
            gV = dVec{};
            g2V = dVec{};
            dMat tMat{}; // tMat(row, col)
            dVec gVdotT{};

            /* compute external potential derivatives */
            gVe = externalPtr->gradV(path(bead1));
            dVe = sqrt(dot(gVe,gVe));
            g2Ve = externalPtr->grad2V(path(bead1));

            gV += gVe;  // update full slice gradient at bead1

            /* Sum potential of bead1 interacting with all other beads at
             * a given time slice.*/
            for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

                /* separation vector */
                rDiff = path.getSeparation(bead1, bead2);
                rmag = sqrt(dot(rDiff,rDiff));

                /* Avoid self interactions */
                if (!(bead1 == bead2)) {

                    /* Compute interaction potential derivatives */
                    gVi = interactionPtr->gradV(rDiff);
                    dVi = sqrt(dot(gVi,gVi));
                    g2Vi = interactionPtr->grad2V(rDiff);
                    
                    /* total derivatives between bead1 and bead2 at bead1 */
                    dV = dVi + dVe;
                    d2V = g2Vi + g2Ve;

                    /* compute the T-matrix for bead1 interacting with bead2 */
                    for (int a=0; a<NDIM; a++){
                        for (int b=0; b<NDIM; b++){
                            tMat(a,b) += rDiff(a)*rDiff(b)*d2V/(rmag*rmag)
                                - rDiff(a)*rDiff(b)*dV/pow(rmag,3);
                            if (a == b)
                                tMat(a,b) += dV/rmag;
                        }
                    }   // end T-matrix 
                    
                    gV += gVi;  // update full slice gradient at bead1

                }   
            }   // end bead2

	    apply_matrix_vector_product(gVdotT, gV, tMat);
            term2 += dot(gVdotT, path(bead1));

        }   // end bead1

        term2 *= (2.0 * gradVFactor[eo] * pow(tau(),3) * constants()->lambda());
    }

    return (term2);
}

/**************************************************************************//**
 *  Return the sum over particles at a given time slice of the 
 *  gradient of the action potential for a single bead dotted into the
 *  bead's position minus the center of mass of the WL that it belongs to. 
 *  
 *  NOTE:  This is the first term.  These were split up because it is
 *      beneficial to be able to return them separately when computing
 *      the specific heat via the centroid virial estimator.
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::deltadotgradUterm1(const int slice) {

    eo = slice % 2;
    int numParticles = path.numBeadsAtSlice(slice);
    int virialWindow = constants()->virialWindow();

    /* The two interacting particles */
    beadLocator bead1, beadNext, beadPrev, beadNextOld, beadPrevOld;
    bead1[0] = bead2[0] = slice;

    dVec gVi, gVe, gV, delta;
    double rDotgV = 0.0;

    /* We loop over the first bead */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        gVi = dVec{};
        gVe = dVec{};
        gV = dVec{};
        delta = dVec{};
        /* Sum potential of bead1 interacting with all other beads at
         * a given time slice.*/
        for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

            /* Avoid self interactions */
            if (!(bead1 == bead2)) {

                /* The interaction component of the force */
                gVi += interactionPtr->gradV(path.getSeparation(bead1,bead2));
            } 
        } // end bead2
        
        /* Now add the external component */
        gVe += externalPtr->gradV(path(bead1));
        
        /* Compute deviation of bead from COM of worldline, 
         * WITHOUT mirror image conv. */
        dVec runTotMore{};
        dVec runTotLess{};
        dVec COM{};
        dVec pos1 = path(bead1);
        beadNextOld = bead1;
        beadPrevOld = bead1;
        for (int gamma = 0; gamma < virialWindow; gamma++) {
            // move along worldline to next (previous) bead
            beadNext = path.next(bead1, gamma);
            beadPrev = path.prev(bead1, gamma);
            // keep running total of distance from first bead
            // to the current beads of interest 
            runTotMore += path.getSeparation(beadNext, beadNextOld);
            runTotLess += path.getSeparation(beadPrev, beadPrevOld);
            // update center of mass of WL
            COM += (pos1 + runTotMore) + (pos1 + runTotLess);
            // store current bead locations
            beadNextOld = beadNext;
            beadPrevOld = beadPrev;
        }
        COM /= (2.0*virialWindow);

        delta = pos1 - COM; // end delta computation.
        path.boxPtr->putInBC(delta);
        
        gV += (gVe + gVi);

        rDotgV += dot(gV, delta);

    } // end bead1

    return (VFactor[eo]*constants()->tau()*rDotgV);
}

/**************************************************************************//**
 *  Return the sum over particles at a given time slice of the 
 *  gradient of the action potential for a single bead dotted into the
 *  bead's position minus the center of mass of the WL that it belongs to. 
 *  
 *  NOTE:  This is the second term
 *
 *  This includes both the external and interaction potentials.
******************************************************************************/
double LocalAction::deltadotgradUterm2(const int slice) {

    eo = slice % 2;
    int numParticles = path.numBeadsAtSlice(slice);
    int virialWindow = constants()->virialWindow();

    /* The two interacting particles */
    beadLocator bead1, beadNext, beadPrev, beadNextOld, beadPrevOld;
    bead1[0] = bead2[0] = slice;

    dVec gVi, gVe, gV, g2V, delta;

    double term2 = 0.0;
    if (gradVFactor[eo] > EPS){

        /* constants for tMatrix */
        dVec rDiff{};
        double rmag = 0.0;
        double d2V = 0.0;
        double dV = 0.0;
        double dVe = 0.0;
        double dVi = 0.0;
        double g2Vi = 0.0;
        double g2Ve = 0.0;

        /* We loop over the first bead */
        for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
            gV    = dVec{};
            g2V   = dVec{};
            delta = dVec{};
            dMat tMat{}; // tMat(row, col)
            dVec gVdotT{};

            /* compute external potential derivatives */
            gVe = externalPtr->gradV(path(bead1));
            dVe = sqrt(dot(gVe,gVe));
            g2Ve = externalPtr->grad2V(path(bead1));

            gV += gVe;  // update full slice gradient at bead1

            /* Sum potential of bead1 interacting with all other beads at
             * a given time slice.*/
            for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

                /* separation vector */
                rDiff = path.getSeparation(bead1, bead2);
                rmag = sqrt(dot(rDiff,rDiff));

                /* Avoid self interactions */
                if (!(bead1 == bead2)) {

                    /* Compute interaction potential derivatives */
                    gVi = interactionPtr->gradV(rDiff);
                    dVi = sqrt(dot(gVi,gVi));
                    g2Vi = interactionPtr->grad2V(rDiff);
                    
                    /* total derivatives between bead1 and bead2 at bead1 */
                    dV = dVi + dVe;
                    d2V = g2Vi + g2Ve;

                    /* compute the T-matrix for bead1 interacting with bead2 */
                    for (int a=0; a<NDIM; a++){
                        for (int b=0; b<NDIM; b++){
                            tMat(a,b) += rDiff(a)*rDiff(b)*d2V/(rmag*rmag)
                                - rDiff(a)*rDiff(b)*dV/pow(rmag,3);
                            if (a == b)
                                tMat(a,b) += dV/rmag;
                        }
                    }   // end T-matrix 
                    
                    gV += gVi;  // update full slice gradient at bead1

                }   
            }   // end bead2

	    apply_matrix_vector_product(gVdotT, gV, tMat);
            
            /* Compute deviation of bead from COM of worldline, 
             * WITHOUT mirror image conv.*/
            dVec runTotMore = 0.0;
            dVec runTotLess = 0.0;
            dVec COM = 0.0;
            dVec pos1 = path(bead1);
            beadNextOld = bead1;
            beadPrevOld = bead1;
            for (int gamma = 0; gamma < virialWindow; gamma++) {
                // move along worldline to next (previous) bead
                beadNext = path.next(bead1, gamma);
                beadPrev = path.prev(bead1, gamma);
                // keep running total of distance from first bead
                // to the current beads of interest
                runTotMore += path.getSeparation(beadNext, beadNextOld);
                runTotLess += path.getSeparation(beadPrev, beadPrevOld);
                // update center of mass of WL
                COM += (pos1 + runTotMore) + (pos1 + runTotLess);
                // store current bead locations
                beadNextOld = beadNext;
                beadPrevOld = beadPrev;
            }
            COM /= (2.0*virialWindow);
            delta = pos1 - COM; // end delta computation.
           
            path.boxPtr->putInBC(delta);

            term2 += dot(gVdotT, delta);

        }   // end bead1

        term2 *= (2.0 * gradVFactor[eo] * pow(tau(),3) * constants()->lambda());
    }

    return (term2);
}

/*************************************************************************//**
 *  Returns the value of the virial kinetic energy correction term.
 *
 * @see S. Jang, S. Jang and G.A. Voth, J. Chem. Phys. 115, 7832 (2001).
******************************************************************************/
double LocalAction::virialKinCorrection(const int slice) {

    double vKC = 0.0;

    /* Combine the potential and a possible correction */
    eo = (slice % 2);

    /* We only add the correction if it is finite */
    if ( gradVFactor[eo] > EPS ) {
        
        vKC += gradVSquared(slice);
        
        /* scale by constants */
        vKC *= gradVFactor[eo] * pow(tau(),3) * constants()->lambda();
    }
    
    return ( vKC );
}

/*************************************************************************//**
 *  Returns the value of the gradient of the potential action for all
 *  beads at a given time slice.
 *
 *  This should be used for the pressure.
******************************************************************************/
dVec LocalAction::gradU(const int slice) {

    dVec gU2 = 0.0;
    dMat tM = tMatrix(slice);
    dVec gV = gradientV(slice);

    /* Combine the potential and a possible correction */
    eo = (slice % 2);
    dVec gU1 = VFactor[eo]*tau()*gradientV(slice);

    /* We only add the correction if it is finite */
    if ( gradVFactor[eo] > EPS ) {
	apply_matrix_vector_product(gU2, gV, tM);
        /* now scale by constants */
        gU2 *= 2.0 * gradVFactor[eo] * pow(tau(),3) * constants()->lambda();
    }
    
    return ( gU1+gU2 );
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// NON LOCAL ACTION BASE CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Setup the path data members for non-local actions.
******************************************************************************/
NonLocalAction::NonLocalAction (const Path &_path, LookupTable &_lookup, 
        PotentialBase *_externalPtr, PotentialBase *_interactionPtr, 
        WaveFunctionBase *_waveFunctionPtr, bool _local, std::string _name) :
    ActionBase(_path,_lookup,_externalPtr,_interactionPtr,_waveFunctionPtr,
            _local,_name) 
{
    NNbead.resize(constants()->initialNumParticles(),false);
}

/**************************************************************************//**
 *  Empty constructor.
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
        totU += sum(U(slice));

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

    double totU = 0.0;

#if PIGS
    /* We tack on a trial wave function and boundary piece if necessary */
    if ( (bead1[0] == 0) || (bead1[0] == (constants()->numTimeSlices()-1)) ) 
            totU -= log(waveFunctionPtr->PsiTrial(bead1));
#endif

    /* Make sure nextBead1 is a real bead and that it is active */
    if ( (nextBead1 == {XXX, XXX}) || (!path.worm.beadOn(nextBead1)) )
        return totU;

    /* Evaluate the external potential */
    double totVext = externalPtr->V(path(bead1))+externalPtr->V(path(nextBead1));

    /* Get the interaction list */
    lookup.updateInteractionList(path,bead1);
    for (int n = 0; n < lookup.numBeads; n++) {
        bead2 = lookup.beadList(n);
        if(!(path.next(bead2) == {XXX, XXX}))
            NNbead[bead2[1]] = true;
    }
    
    /* Get the next interaction list */
    lookup.updateInteractionList(path,nextBead1);
    for (int n = 0; n < lookup.numBeads; n++) {
        bead2 = path.prev(lookup.beadList(n));
        if(!(bead2 == {XXX, XXX}))
           NNbead[bead2[1]] = true;
    }
    
    bead2[0] = bead1[0];
    for (bead2[1]= 0; bead2[1] < path.numBeadsAtSlice(bead1[0]); bead2[1]++) {
        if(NNbead[bead2[1]]){
            nextBead2 = path.next(bead2);
            sep = path.getSeparation(bead1,bead2);
            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->V(sep,sep2);
            NNbead[bead2[1]] = false;
        }
    }
    
    double totUext = tau()*totVext/(2.0);

    /* bead2[0] = bead1[0]; */ 
    
    /* /1* Now calculate the total effective interation potential, neglecting self-interactions *1/ */
    /* for (bead2[1]= 0; bead2[1] < path.numBeadsAtSlice(bead1[0]); bead2[1]++) { */
        
    /*     /1* Skip self interactions *1/ */
    /*     if ( bead2[1] != bead1[1] ) { */
            
    /*         /1* get the separation between the two particles at the first time */
    /*          * slice *1/ */
    /*         sep = path.getSeparation(bead1,bead2); */
            
    /*         /1* Now we find the separation at the advanced time step *1/ */
    /*         nextBead2 = path.next(bead2); */
            
    /*         /1* If the imaginary time neighbor exists, compute the effective */
    /*          * potential *1/ */
    /*         if ( (!(nextBead2 == {XXX, XXX})) && (path.worm.beadOn(nextBead2)) ) { */
    /*             sep2 = path.getSeparation(nextBead1,nextBead2); */
    /*             totU += interactionPtr->V(sep,sep2,constants()->lambda(),constants()->tau()); */
    /*         } */
    /*     } // bead2 != bead1 */
        
        
    /* } // for bead2 */
    
    return (totU + totUext);
}

/**************************************************************************//**
 *  Return the potential action for all time slices and all particles.
 *
 *  Computes the total potential energy by summing over all particles and time
 *  slices.  
******************************************************************************/
std::array<double,2> NonLocalAction::U(int slice) {

    double totUint = 0.0;
    double totUext = 0.0;

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
            sep = path.getSeparation(bead1,bead2);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totUint += interactionPtr->V(sep,sep2);
        } // bead2

    } // bead1
    return std::array<double,2>(totUext,totUint);
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
            sep = path.getSeparation(bead1,bead2);
            updateSepHist(sep);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->dVdtau(sep,sep2);
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

    /* Calculate the total potential, including external and interaction
     * effects*/
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

        /* Get the advanced neightbor of bead1 */
        nextBead1 = path.next(bead1);

        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead1,bead2);

            /* Get the advanced neighbor of the second bead */
            nextBead2 = path.next(bead2);

            sep2 = path.getSeparation(nextBead1,nextBead2);
            totU += interactionPtr->dVdlambda(sep,sep2);
        } // bead2

    } // bead1
    return ( totU );
}
