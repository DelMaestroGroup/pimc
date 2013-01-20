/**
 * @file estimator.cpp
 * @author Adrian Del Maestro
 *
 * @brief Estimator class implementations.
 */

#include "estimator.h"
#include "path.h"
#include "action.h"
#include "potential.h"
#include "communicator.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ESTIMATOR BASE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 * Initialize and specify the output file.
******************************************************************************/
EstimatorBase::EstimatorBase (const Path &_path, int _frequency, string _label) : 
    path(_path),
    name(),
    label(_label),
    numSampled(0),
	numAccumulated(0),
	totNumAccumulated(0),
    frequency(_frequency),
    diagonal(true),
    endLine(true)
{
    /* Two handy local constants */
	canonical = constants()->canonical();
	numBeads0 = constants()->initialNumParticles()*constants()->numTimeSlices();
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
EstimatorBase::~EstimatorBase() { 
	estimator.free();
    norm.free();
}

/**************************************************************************//**
 *  Determine the basic sampling condition.
 * 
 *  We only sample when:
 *  1) frequency > 0
 *  2) numSampled is a even factor of frequency
 *  3) The configuration (diagonal/off diagonal) matches the estimator type.
 *  4) If we are canonical, then only when we match the desired number of
 *  particles
******************************************************************************/
bool EstimatorBase::baseSample() {

    /* Increment the number of attempted samplings */
    numSampled++;

    if (!frequency)
        return false;
    if ((numSampled % frequency) != 0)
        return false;
    if (!(path.worm.isConfigDiagonal == diagonal))
        return false;
    if (!canonical)
        return true;
	if (path.worm.getNumBeadsOn() == numBeads0)
        return true;

    return false;
}

/**************************************************************************//**
 *  Sample the estimator.
 * 
 *  Here we simply call accumulate every frequency times sample is called,  
 *  depending on whether we are measuring a diagonal or off-diagonal estimator.
 *  If frequency == 0, we don't bother to measure it.
******************************************************************************/
void EstimatorBase::sample() {

    if (baseSample()) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
}

/**************************************************************************//**
 *  Initialize estimator. 
 *
 *  Initilize the estimator and normalization arrays.
******************************************************************************/
void EstimatorBase::initialize(int _numEst) {

    numEst = _numEst;
	estimator.resize(numEst);
	norm.resize(numEst);
	norm = 1.0;
	reset();
}

/**************************************************************************//**
 *  Prepare estimator for i/o. 
 *
 *  Assign the output file and write a header depending on whether or not
 *  we are restarting the simulation.
******************************************************************************/
void EstimatorBase::prepare() {

    /* Provided that we will perform at least one measurement, open an output
     * file and possibly write a header */
    if (frequency > 0) {
        /* Assign the output file pointer */
        outFilePtr = &(communicate()->file(label)->stream());

        /* Write the header to disk if we are not restarting */
        if (!constants()->restart()) {
            (*outFilePtr) << header;
            if (endLine)
                (*outFilePtr) << endl;
        }
    }

}

/*************************************************************************//**
 *  Reset numAccumulated and the estimator to 0.
******************************************************************************/
void EstimatorBase::reset() {
	numAccumulated = 0;
	estimator = 0.0;
}

/*************************************************************************//**
 *  Restart the measurment process from a previous state
******************************************************************************/
void EstimatorBase::restart(const uint32 _numSampled, const uint32 _totNumAccumulated) {
	numSampled = _numSampled;
	totNumAccumulated = _totNumAccumulated;
	reset();
}

/*************************************************************************//**
 *  Output the estimator value to disk.  
 *
 *  We only attempt an output if we have accumulated at least numBins 
 *  measurements.  This works because the final 'scalar' estimator is the 
 *  diagonal estimator which gets measured for both diaagonal and off diagonal 
 *  configurations.
******************************************************************************/
void EstimatorBase::output() {

	/* Average */
	estimator *= (norm/(1.0*numAccumulated));

	/* Now write the estimator values to disk */
	for (int n = 0; n < numEst; n++) 
		(*outFilePtr) << format("%16.8E") % estimator(n);

	if (endLine)
		(*outFilePtr) << endl;

	/* Reset all values */
	reset();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ENERGY ESTIMATOR CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 * We measure the total Kinetic, Potential and E-mu*N as well as the kinetic,
 * potential and total energy per particle.
******************************************************************************/
EnergyEstimator::EnergyEstimator (const Path &_path, ActionBase *_actionPtr,
		int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label), actionPtr(_actionPtr) {

	/* Set estimator name and header, we will always report the energy
	 * first, hence the comment symbol*/
	name = "Energy";
	header = str(format("#%15s%16s%16s%16s%16s%16s%16s") 
			% "K" % "V" % "E" % "E_mu" % "K/N" % "V/N" % "E/N");
    endLine = false;
    initialize(7);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
EnergyEstimator::~EnergyEstimator() { 
}

/*************************************************************************//**
 *  Accumluate the energy.
 *
 *  We use the thermodynamic estimator for the kinetic energy and the
 *  potential estimator for the potential energy. A possible shift is 
 *  made due to a tail correction.
******************************************************************************/
void EnergyEstimator::accumulate() {

	double totK = 0.0;
	double totV = 0.0;
	double gV2 = 0.0;

	int numParticles  = path.getTrueNumParticles();
	int numTimeSlices = path.numTimeSlices;

	/* The total tail correction */
	double tailV = (1.0*numParticles*numParticles/path.boxPtr->volume)
		* actionPtr->potentialPtr->getTailV();

	/* The kinetic normalization factor */
	double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * numTimeSlices);

	/* The action correction normalization factor */
	double normAct = 2.0 * constants()->lambda() * constants()->tau() * constants()->tau()
		/ (1.0 * numTimeSlices);

	/* The classical contribution to the kinetic energy per particle 
	 * including the chemical potential */
	double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;

	/* We first calculate the kinetic energy.  Even though there
	 * may be multiple mixing and swaps, it doesn't matter as we always
	 * just advance one time step at a time, as taken care of through the
	 * linking arrays.  This has been checked! */
	beadLocator beadIndex;
	dVec vel;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
			beadIndex = slice,ptcl;
			vel = path.getVelocity(beadIndex);
			totK -= dot(vel,vel);
		}
	}

	/* Normalize the accumulated link-action part */
	totK *= kinNorm;

	/* Now we compute the potential and kinetic energy.  We use an operator estimater
	 * for V and the thermodynamic estimator for K */
	int eo;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		eo = (slice % 2);
//		totV += (actionPtr->pFactor[eo]) * (actionPtr->potentialPtr->V(slice));
		if (actionPtr->eFactor[eo] > EPS) 
			gV2  += (actionPtr->eFactor[eo]) * (actionPtr->potentialPtr->gradVSquared(slice));
		else if (eo==0)
			totV  += actionPtr->potentialPtr->V(slice);
	}

//	/* Normalize the action correction and the total potential*/
//	gV2 *= normAct;
//	totV /= (1.0 * numTimeSlices);
//	totV += 2.0*gV2;

	/* Normalize the action correction and the total potential*/
	gV2 *= normAct;
	totV /= (0.5 * numTimeSlices);

	/* Perform all the normalizations and compute the individual energy terms */
	totK  += (classicalKinetic + gV2);

	totV += tailV;

	/* Now we accumulate the average total, kinetic and potential energy, 
	 * as well as their values per particles. */
	estimator(0) += totK;
	estimator(1) += totV;
	estimator(2) += totK + totV;

	estimator(3) += totK + totV - constants()->mu()*numParticles;

	estimator(4) += totK/(1.0*numParticles);
	estimator(5) += totV/(1.0*numParticles);
	estimator(6) += (totK + totV)/(1.0*numParticles);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// NUM PARTICLES ESTIMATOR CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  We measure the number of particles, (number of particles)^2 and the
 *  density of particles.
******************************************************************************/
NumberParticlesEstimator::NumberParticlesEstimator (const Path &_path, 
        int _frequency, string _label) :
	EstimatorBase(_path,_frequency,_label) {

	/* Set estimator name and header */
	name = "Number Particles";
	header = str(format("%16s%16s%16s") % "N" % "N^2" % "density");
    endLine = false;
    initialize(3);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
NumberParticlesEstimator::~NumberParticlesEstimator() { 
}

/*************************************************************************//**
 * Accumulate the number of particles and density.
******************************************************************************/
void NumberParticlesEstimator::accumulate() {
	int numParticles = path.getTrueNumParticles();
	estimator(0) += 1.0*numParticles;
	estimator(1) += 1.0*numParticles*numParticles;
	estimator(2) += 1.0*numParticles/path.boxPtr->volume;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// NUMBER DISTRIBUTION ESTIMATOR CLASS ---------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  The number of particles probability distribution is setup.  We allow 
 *  for particle fluctuations up to 50 particles away from the initial value
 *  specified.
******************************************************************************/
NumberDistributionEstimator::NumberDistributionEstimator (const Path &_path, 
        int _frequency, string _label) :
	EstimatorBase(_path,_frequency,_label) {

	/* For now, we assume 50 particles on each side of the mean */
	particleShift = 50;
	startParticleNumber = max(constants()->initialNumParticles()-particleShift,0);
	if (startParticleNumber == 0)
		endParticleNumber = 2*particleShift;
	else
		endParticleNumber = constants()->initialNumParticles() + particleShift;
	maxNumParticles = endParticleNumber - startParticleNumber + 1;

	/* If our number of paticles is too small, we compensate */
	if ((constants()->initialNumParticles() - particleShift) < 0)
		particleShift = constants()->initialNumParticles();

    initialize(maxNumParticles);

	/* Set estimator name and header */
	name = "Number Distribution";
	header = str(format("#%15d") % startParticleNumber);
	for (int n = startParticleNumber+1; n <= endParticleNumber; n++) 
		header.append(str(format("%16d") % n));
}

/*************************************************************************//**
 * Destructor.
******************************************************************************/
NumberDistributionEstimator::~NumberDistributionEstimator() { 
}

/*************************************************************************//**
 * Accumulate the number probability distribution
******************************************************************************/
void NumberDistributionEstimator::accumulate() {

	/* Get the correct particle Number index, and increment the corresponding bin */
	int index = path.getTrueNumParticles()-constants()->initialNumParticles() 
        + particleShift;
	if (index >= 0 && index < maxNumParticles)
		estimator(index) += 1.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PARTICLE DENSITY ESTIMATOR CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  A full NDIM-dimensional particle density hisgram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
ParticlePositionEstimator::ParticlePositionEstimator (const Path &_path,
        int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

    initialize(path.boxPtr->numGrid);

    /* Set estimator name and header. */
    name = "Particle Position";
    header = str(format("#%15d\n#%15s") % NGRIDSEP % "density");

    /* The normalization: 1/(dV*M) */
    for (int n = 0; n < numEst; n++)
        norm(n) = 1.0/(1.0*path.numTimeSlices*path.boxPtr->gridBoxVolume(n));
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
ParticlePositionEstimator::~ParticlePositionEstimator() { 
}

/*************************************************************************//**
 *  Overload the output of the base class so that a running average
 *  is kept rather than keeping all data.
******************************************************************************/
void ParticlePositionEstimator::output() {

    /* Prepare the position file for writing over old data */
    communicate()->file(label)->reset();

    (*outFilePtr) << header;
    if (endLine)
        (*outFilePtr) << endl;

	/* Now write the running average of the estimator to disk */
	for (int n = 0; n < numEst; n++) 
        (*outFilePtr) << format("%16.8E\n") % 
            (norm(n)*estimator(n)/totNumAccumulated);

    communicate()->file(label)->rename();
}

/*************************************************************************//**
 *  Accumulate a histogram of all particle positions, with output 
 *  being the running average of the density per grid space.
******************************************************************************/
void ParticlePositionEstimator::accumulate() {

    beadLocator beadIndex;

    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;

            /* update our particle position histogram */
            int n = path.boxPtr->gridIndex(path(beadIndex));
            estimator(n) += 1.0;
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SUPERFLUID FRACTION ESTIMATOR CLASS ---------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Measure the superfluid fraction from the value of the winding number in
 *  all directions with periodic boundary conditions.  We also measure
 *  the individual values of the winding numbers and their probability
 *  distribution.
******************************************************************************/
SuperfluidFractionEstimator::SuperfluidFractionEstimator (const Path &_path, 
		int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

	windMax = 10;
	/* We compute a bunch of estimators here, the superfluid fraction, the winding
	 * number in all possible dimensions, and the winding number histograms up to 
	 * windMax windings. These are all diagonal estimators and we have our own
	 * output file.*/
    initialize(4 + 2*windMax + 1 + 1);

	/* Set estimator name */
	name = "Superfluid Fraction";
	header = str(format("#%15s%16s%16s%16s") % "rho_s/rho" % "W^2(x)" % "W^2(y)" % "W^2(z)");
	for (int w = -windMax; w <= windMax; w++)
		header += str(format("%11sP(%+1d)") % " " % w);
    header += str(format("%16s") % "Area_rho_s");

	/* The pre-factor for the superfluid density is always the same */
	norm(0) = constants()->T() / (2.0 * sum(path.boxPtr->periodic) * constants()->lambda());

    /* The pre-factor for the area esimator */
    norm(5+2*windMax) = 0.5*constants()->T()*constants()->numTimeSlices()/constants()->lambda();
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
SuperfluidFractionEstimator::~SuperfluidFractionEstimator() { 
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  We measure the winding number W, as well as the W^2/N which is used to 
 *  compute the superfluid fraction.  The probability distribution of winding
 *  number fluctuations is also computed.
******************************************************************************/
void SuperfluidFractionEstimator::accumulate() {

	int numTimeSlices= path.numTimeSlices;
	double locW2oN = 0.0;

	/* Sum up the winding number over all particles */
	beadLocator beadIndex;
    double Az, I;
    dVec pos1,pos2;

    Az = I = 0.0;
	dVec W,vel;
	W = 0.0;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            /* The winding number estimator */
			beadIndex = slice,ptcl;
			vel = path.getVelocity(beadIndex);
			W += vel;

            /* The area estimator */
			pos1 = path(beadIndex);
			pos2 = path(path.next(beadIndex));
			Az += pos1[0]*pos2[1]-pos2[0]*pos1[1];
            I +=  pos1[0]*pos2[0] + pos1[1]*pos2[1];

		}
	}

	/* Scale by the periodicity of the boundary conditions */
	W *= path.boxPtr->periodic;

	/* Compute the locally scaled W^2/N */
	locW2oN = dot(W,W)/(1.0*path.getTrueNumParticles());

	/* The average winding number squared */
	estimator(0) += locW2oN;

	/* Scale by the length of the system in each dimension*/
	W *= path.boxPtr->sideInv;

	/* The individual winding numbers, we always store 3 values regardless
	 * of the dimensionality to ensure output file consistency */
	int i;
	for (i = 0; i < NDIM; i++)
		estimator(1+i) += W[i]*W[i];
	for (int j = i; j < 3; j++)
		estimator(1+j) += 0.0;

	/* Calcluate the winding number probablity in the NDIM^th dimension */
	int n = 0;
	for (int p = -windMax; p <= windMax; p++) {
		if (abs(W[NDIM-1]-1.0*p) < 0.2)
			estimator(4+n) += 1.0;
		++n;
	}

    /* The Area Estimator */
    estimator(5+2*windMax) += Az*Az/I;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LOCAL SUPERFLUID DENSITY ESTIMATOR CLASS ----------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Initialize the variables needed to measure the local superfluid density 
 *  as described in: 
 *  @see E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).
 *
******************************************************************************/
LocalSuperfluidDensityEstimator::LocalSuperfluidDensityEstimator
(const Path &_path, int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

    /* This is a 'local' histogram estimator so we use the defined grid */
    initialize(2*path.boxPtr->numGrid);

    /* Set estimator name and header. */
    name = "Local Superfluid";
    header = str(format("#%15d\n") % NGRIDSEP);
    header += str(format("#%15s%16s") % "W:rho_s" % "A:rho_s");

    /* The normalization: 1/(dV) */
//    for (int n = 0; n < numEst/2; n++) {
//        norm(n) = 1.0/(path.boxPtr->gridBoxVolume(n));
//        norm(n+numEst/2) = 1.0/(path.boxPtr->gridBoxVolume(n));
//    }

    /* The normalization constant for the winding and area estimator. */
    for (int n = 0; n < numEst/2; n++) {
        double dV = path.boxPtr->gridBoxVolume(n);
        norm(n) = 0.5*constants()->T()/(dV*constants()->lambda());
        norm(n+numEst/2) = 0.5*constants()->T()/(dV*constants()->lambda());
    }

    /* Initialize the local arrays */
    locW.resize(numEst/2);
    locW = 0.0;

    locAz.resize(numEst/2);
    locAz = 0.0;

}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
LocalSuperfluidDensityEstimator::~LocalSuperfluidDensityEstimator() { 
    locW.free();
    locAz.free();
}

/*************************************************************************//**
 *  Overload the output of the base class so that a running average
 *  is kept rather than keeping all data.
******************************************************************************/
void LocalSuperfluidDensityEstimator::output() {

    /* Prepare the position file for writing over old data */
    communicate()->file(label)->reset();

    (*outFilePtr) << header;
    if (endLine)
        (*outFilePtr) << endl;

	/* Now write the running average of the estimator to disk */
	for (int n = 0; n < numEst/2; n++) 
        (*outFilePtr) << format("%16.8E%16.8E\n") % 
            (norm(n)*estimator(n)/totNumAccumulated) %
            (norm(n+numEst/2)*estimator(n+numEst/2)/totNumAccumulated);

    communicate()->file(label)->rename();
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  INSERT REFERENCES AND DESCRIPTION
 *  
 *  
******************************************************************************/
void LocalSuperfluidDensityEstimator::accumulate() {

	int numTimeSlices = path.numTimeSlices;
    locAz = 0.0;
    locW = 0.0;

	beadLocator beadIndex;
    double Az, I, rp2;
    dVec pos1,pos2;
	dVec W,vel;
	W = 0.0;

    Az = I = 0.0;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
			pos1 = path(beadIndex);
			pos2 = path(path.next(beadIndex));
            int n = path.boxPtr->gridIndex(pos1);

            /*  The distance from the z-axis squared */
            rp2 = pos1[0]*pos1[0] + pos1[1]*pos1[1];

            /* The winding number estimator */
			vel = path.getVelocity(beadIndex)*path.boxPtr->periodic;
			//W += vel;

            /* The local part of the winding number */
            locW(n) += dot(vel,vel);

            /* The z-component of the area estimator */
			Az += pos1[0]*pos2[1] - pos2[0]*pos1[1];
            //I +=  pos1[0]*pos2[0] + pos1[1]*pos2[1];

            /* Get the local part of the path area */
            locAz(n) += (pos1[0]*pos2[1] - pos2[0]*pos1[1])/rp2;
		}
	}

    locW /= 1.0*path.getTrueNumParticles();
    //locAz *= path.getTrueNumParticles()*Az/I;
    locAz *= Az;
    for (int n = 0; n < numEst/2; n++) {
        //estimator(n) += dot(locW(n),W);
        estimator(n) += locW(n);
        estimator(n+numEst/2) += locAz(n);
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// DIAGONAL FRACTION ESTIMATOR CLASS -----------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 *
 *  This esimator determines the fraction of the simulation spent in the 
 *  diagonal ensemble (no worms).
******************************************************************************/
DiagonalFractionEstimator::DiagonalFractionEstimator (const Path &_path, 
		int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

    initialize(1);

	/* Set estimator name */
	name = "Diagonal Fraction";
	header = str(format("%16s") % "diagonal");

}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
DiagonalFractionEstimator::~DiagonalFractionEstimator() { 
}

/*************************************************************************//**
 *  Here we just test and see if we are diagonal,  if so, accumulate the
 *  counter.
******************************************************************************/
void DiagonalFractionEstimator::accumulate() {
	if (path.worm.isConfigDiagonal)
		estimator(0) += 1.0;
}

/*************************************************************************//**
 *  Overload sampling to make sure it is always done, regardless of ensemble.
******************************************************************************/
void DiagonalFractionEstimator::sample() {

	numSampled++;

    if (frequency && ((numSampled % frequency) == 0)) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// WORM ESTIMATOR CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  We measure various worm properties including the relative worm length, 
 *  relative worm gap, worm cost, the separation between head and tail, and
 *  the number of beads / number of time slices.
******************************************************************************/
WormPropertiesEstimator::WormPropertiesEstimator (const Path &_path, 
        int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label) {

	/* We measure the average worm length, gap and cost.  It is an off-diagonal
	 * estimator that is output to its own file */
    initialize(5);
    diagonal = false;

	/* Set estimator name */
	name = "Worm Properties";
	header = str(format("#%15s%16s%16s%16s%16s") % "rel-worm-len" % 
			"rel-worm-gap" % "worm-cost" % "head-tail-sep" % "particles");
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
WormPropertiesEstimator::~WormPropertiesEstimator() { 
}

/*************************************************************************//**
 *  Accumulate the length, gap and cost of the current worm.
******************************************************************************/
void WormPropertiesEstimator::accumulate() {
	estimator(0) += 1.0*path.worm.length / (1.0*path.numTimeSlices);
	estimator(1) += 1.0*path.worm.gap / (1.0*path.numTimeSlices);
	double r = dot(path.worm.sep,path.worm.sep);
	if (path.worm.gap == 0)
		estimator(2) += 0.0;
	else
		estimator(2) += 1.0*r*constants()->fourLambdaTauInv()/(1.0*path.worm.gap);
	estimator(3) += 1.0*sqrt(r);
	estimator(4) += 1.0*path.worm.getNumBeadsOn()/(1.0*constants()->numTimeSlices());
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PERMUTATION CYCLE ESTIMATOR CLASS -----------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 * Constructor.
 * 
 * We setup the permuatation cycle estimator to measure cycles which can 
 * contain up to 40 particles.
******************************************************************************/
PermutationCycleEstimator::PermutationCycleEstimator(const Path &_path, 
		int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

	/* We just choose arbitrarily to only count cycles including up to 40 particles */
	maxNumCycles = 40;
	
	/* The permutatoin cycle estimator has its own file, and consists of 
	 * maxNumCycles permutation cycles */
    initialize(maxNumCycles);

	/* Set estimator name and header, which contains the permutation cycle
	 * numbers */
	name = "Permutation Cycle";
	header = str(format("#%15d") % 1);
	for (int n = 2; n <= maxNumCycles; n++) 
		header.append(str(format("%16d") % n));
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
PermutationCycleEstimator::~PermutationCycleEstimator() { 
	doBead.free();
}

/*************************************************************************//**
 * Accumulate permuation cycle.
 * 
 * Go through each worldline and count how many beads it includes. We must be
 * careful to avoid any double counting, which is acomplished by marking
 * beads we have touched.
******************************************************************************/
void PermutationCycleEstimator::accumulate() {

	/* The start bead for each world line, and the moving index */
	beadLocator startBead;
	beadLocator beadIndex;

	int numParticles = path.getTrueNumParticles();
	int numWorldlines = path.numBeadsAtSlice(0);
	double cycleNorm = 1.0 / (1.0*numParticles);

	/* We create a local vector, which determines whether or not we have
	 * already included a bead at slice 0*/
	doBead.resize(numWorldlines);
	doBead = true;

	/* We go through each particle/worldline */
	for (int n = 0; n < numWorldlines; n++) {

		/* The initial bead to be moved */
		startBead = 0,n;

		/* We make sure we don't try to touch the same worldline twice */
		if (doBead(n)) {

			/* Mark the beads as touched and increment the number of worldlines */
			beadIndex = startBead;

			/* The world line length, we simply advance until we have looped back on 
			 * ourselves. */
			int wlLength = 0;
			do {
				wlLength++;

				/* We turn off any zero-slice beads we have touched */
				if (beadIndex[0]==0)
					doBead(beadIndex[1]) = false;

				beadIndex = path.next(beadIndex);
			} while (!all(beadIndex==startBead));

			/* Accumulte the cycle length counter */
			int cycleNum = int(wlLength / path.numTimeSlices);
			if ((cycleNum > 0) && (cycleNum <= maxNumCycles)) 
				estimator(cycleNum-1) += 1.0*cycleNum*cycleNorm;
		} // doBead

	} // n
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ONE BODY DENSITY MATRIX ESTIMATOR CLASS -----------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  The one body density matrix estimator is initialized.  We measure NOBDMSEP
 *  positions, out to the maximum separation in the sample (which may depend
 *  on the type of simulation cell).
******************************************************************************/
OneBodyDensityMatrixEstimator::OneBodyDensityMatrixEstimator (Path &_path,
		ActionBase *_actionPtr, MTRand &_random, int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label), 
	lpath(_path),
	actionPtr(_actionPtr),
	random(_random)
{

	sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());

	/* We chooose the maximum separation to be sqrt(NDIM)*min(L)/2 */
	dR = 0.5*sqrt(sum(path.boxPtr->periodic))*(blitz::min(path.boxPtr->side)) / (1.0*NOBDMSEP);

	/* This is an off-diagonal estimator*/
    initialize(NOBDMSEP);
    diagonal = false;

	/* Set estimator name */
	name = "One Body Density Matrix";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NOBDMSEP; n++) 
		header.append(str(format("%16.3E") % (n*dR)));

	numReps = 5;
	norm = 1.0 / (1.0*numReps);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
OneBodyDensityMatrixEstimator::~OneBodyDensityMatrixEstimator() { 
}

/*************************************************************************//**
 *  Sample the OBDM.
 * 
 *  We overload the sample method for the one body density matrix as we
 *  only want to measure when the gap is not too large, otherwise we will be
 *  dominated by tiny close probabilities.
******************************************************************************/
void OneBodyDensityMatrixEstimator::sample() {
	numSampled++;
	if ( frequency && 
         ((numSampled % frequency) == 0) &&
         (path.worm.isConfigDiagonal == diagonal) &&
         (path.worm.gap > 0) && (path.worm.gap <= constants()->Mbar())  &&
         (actionPtr->eFactor[(lpath.worm.tail[0] % 2)] < EPS) ) {

		/* If we are canonical, we want the closed configuration to be close
		 * to our ideal one */
		if ( (!canonical) || 
			 (abs(path.worm.getNumBeadsOn()+path.worm.gap-numBeads0) <= 2) ) {
			totNumAccumulated++;
			numAccumulated++;
			accumulate();
		}
	}
}

/*************************************************************************//**
 *  Return a dimensionally dependent random vector of length r.  
 *
 *  If we are in 3D in a cylinder geometry, we only shift in the z-direction.
 *  @param r The length of the random vector
 *  @return a random NDIM-vector of length r
******************************************************************************/
inline dVec OneBodyDensityMatrixEstimator::getRandomVector(const double r) {
	dVec rVec;
	rVec = 0.0;
#if NDIM==1
	if (random.rand() < 0.5)
		rVec = r;
	else
		rVec = -r;
#elif NDIM==2
	double theta = 2.0*M_PI*random.rand();
	rVec[0] = r*cos(theta);
	rVec[1] = r*sin(theta);
#elif NDIM==3
	if (lpath.boxPtr->name == "Prism") {
		double theta = 2.0*M_PI*random.rand();
		double phi   = M_PI*random.rand();
		rVec[0] = r*cos(theta)*sin(phi);
		rVec[1] = r*sin(theta)*sin(phi);
		rVec[2] = r*cos(phi);
	} 
	else {
		if (random.rand() < 0.5)
			rVec[NDIM-1] = r;
		else
			rVec[NDIM-1] = -r;
	}
#endif
	return rVec;
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
dVec OneBodyDensityMatrixEstimator::newStagingPosition(const beadLocator &neighborIndex, 
		const beadLocator &endIndex, const int stageLength, const int k) {

    /* The rescaled value of lambda used for staging */
    double f1 = 1.0 * (stageLength - k - 1);
    double f2 = 1.0 / (1.0*(stageLength - k));
    double sqrtLambdaKTau = sqrt2LambdaTau * sqrt(f1 * f2);

	/* We find the new 'midpoint' position which exactly samples the kinetic 
	 * density matrix */
	neighborPos = lpath(neighborIndex);
	newRanPos = lpath(endIndex) - neighborPos;
	lpath.boxPtr->putInBC(newRanPos);
	newRanPos *= f2;
	newRanPos += neighborPos;

	/* This is the random kick around that midpoint */
	for (int i = 0; i < NDIM; i++)
		newRanPos[i] = random.randNorm(newRanPos[i],sqrtLambdaKTau);

    lpath.boxPtr->putInside(newRanPos);

    return newRanPos;
}

/*************************************************************************//**
 *  Accumulate the OBDM.
 * 
 *  We perform a fake close move, where the head of the worm is advanced to
 *  a position a distance 'r' away from the tail but at the same time slice.
 *  The probability of excepting such a move is equal (up to normalization)
 *  to the one body density matrix.
******************************************************************************/
void OneBodyDensityMatrixEstimator::accumulate() {

	oldTailPos = lpath(lpath.worm.tail);
	oldAction = actionPtr->potentialAction(lpath.worm.tail);

	/* We make a list of all the beads involved in the move, adding them
	 * as we go. */
	beadLocator beadIndex;
	beadIndex = lpath.worm.head;
	dVec pos;
	pos = 0.0;
	for (int k = 0; k < (lpath.worm.gap-1); k++) 
		beadIndex = lpath.addNextBead(beadIndex,pos);

	/* Perform the final connection to the tail*/
	lpath.next(beadIndex) = lpath.worm.tail;
	lpath.prev(lpath.worm.tail) = beadIndex;

	for (int p = 0; p < numReps; p++) {

		/* Now we loop through all possible separations, evaluating the potential
		 * action */
		for (int n = 0; n < NOBDMSEP; n++) {

			newAction = 0.0;
			++numAttempted;

			/* Assign the new displaced tail position */
			newTailPos = oldTailPos + getRandomVector(n*dR);
			lpath.boxPtr->putInside(newTailPos);
			lpath.updateBead(lpath.worm.tail,newTailPos);

			/* Compute the free particle density matrix */
			rho0Norm = actionPtr->rho0(lpath.worm.head,lpath.worm.tail,lpath.worm.gap);

			/* action shift coming from a finite chemical potential */
			double muShift = lpath.worm.gap*constants()->mu()*constants()->tau();

			/* Starting from the head, we generate new positions for the beads via
			 * staging, and accumulate the potential action */
			beadIndex = lpath.worm.head;
			int k = 0;
			do {
				if (!all(beadIndex==lpath.worm.head) && !all(beadIndex==lpath.worm.tail)) {
					lpath.updateBead(beadIndex,
							newStagingPosition(path.prev(beadIndex),lpath.worm.tail,lpath.worm.gap,k));
					++k;
				}

				newAction += actionPtr->potentialAction(beadIndex);

				beadIndex = lpath.next(beadIndex);
			} while (!all(beadIndex==lpath.next(lpath.worm.tail)));

			double expAction = exp(-newAction + oldAction + muShift);

			estimator(n) += rho0Norm*expAction;

			/* Record the probability of accepting the move */
			if (random.randExc() < rho0Norm*expAction)
				++numAccepted;

		} // end for n

	} // end for k

	/* Now we must undo any damge we have caused by reverting the tail to its previous position,
	 * and turning off all intermediate beads */
	lpath.updateBead(lpath.worm.tail,oldTailPos);

	/* Delete all the beads that were added. */
	beadIndex = lpath.next(lpath.worm.head);
	while (!all(beadIndex==lpath.worm.tail)) {
		beadIndex = lpath.delBeadGetNext(beadIndex);
	}
	lpath.next(lpath.worm.head) = XXX;
	lpath.prev(lpath.worm.tail) = XXX;
}

/*************************************************************************//**
 *  For the one body density matrix estimator, we would like to output
 *  the acceptance information for the accumulate trial move.
******************************************************************************/
void OneBodyDensityMatrixEstimator::outputFooter() {

	(*outFilePtr) << format("# accepted: %16.8E attempted: %16.8E ratio: %16.4E\n") 
		% (1.0*numAccepted) % (1.0*numAttempted) % (1.0*numAccepted/(1.0*numAttempted));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PAIR CORRELATION ESTIMATOR CLASS ------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  For the pair correlation function, we measure NPCFSEP positions to get
 *  high enough data density to observe possible osscilations.  The normalization
 *  depends on dimension.
******************************************************************************/
PairCorrelationEstimator::PairCorrelationEstimator (const Path &_path, 
		ActionBase *_actionPtr, int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label),
	actionPtr(_actionPtr)
{
	/* The spatial discretization */
	dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

	/* This is a diagonal estimator that gets its own file */
	initialize(NPCFSEP);

	/* Set estimator name */
	name = "Pair Correlation Function";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NPCFSEP; n++) 
		header.append(str(format("%16.3E") % ((n)*dR)));

	/* The normalization factor for the pair correlation function depends 
	 * on the dimensionality */
//	TinyVector<double,3> gNorm;
//	gNorm[0] = 0.5;
//	gNorm[1] = 1.0/(4.0*M_PI);
//	gNorm[2] = 1.0/(8.0*M_PI);
//	norm(0) = 1.0;
//	for (int n = 1; n < NPCFSEP; n++)
//		norm(n) = (gNorm[NDIM-1]*path.boxPtr->volume) / (dR*pow(n*dR,NDIM-1));

	/* The normalization factor for the pair correlation function depends 
	 * on the dimensionality, and container type */
	if (path.boxPtr->name == "Cylinder") {
		for (int n = 0; n < NPCFSEP; n++)
			norm(n) = 0.5*path.boxPtr->side[NDIM-1] / dR;
	}
	else {
		TinyVector<double,3> gNorm;
		gNorm[0] = 1.0;
		gNorm[1] = 1.0/(M_PI);
		gNorm[2] = 3.0/(2.0*M_PI);
		double dV;
		for (int n = 0; n < NPCFSEP; n++) {
			dV = pow((n+1)*dR,NDIM)-pow(n*dR,NDIM);
			norm(n) = 0.5*(gNorm[NDIM-1]*path.boxPtr->volume) / dV;
		}
	}
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
PairCorrelationEstimator::~PairCorrelationEstimator() { 
}

/*************************************************************************//**
 *  We simply update the local value of the esitimator with the value from
 *  our potential reference.  
 *
 *  We only compute this for N > 1.
******************************************************************************/
void PairCorrelationEstimator::accumulate() {
	int numParticles = path.getTrueNumParticles();
	double lnorm = 1.0*(numParticles-1)/(1.0*numParticles);
	if (numParticles > 1) {
		estimator += lnorm*(1.0*actionPtr->potentialPtr->sepHist / 
			(1.0*sum(actionPtr->potentialPtr->sepHist)));
	}
	else
		estimator += 1.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RADIAL DENSITY ESTIMATOR CLASS --------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  The radial density is binned into NRADSEP boxes.
******************************************************************************/
RadialDensityEstimator::RadialDensityEstimator (const Path &_path, 
        int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label)
{
	/* The spatial discretization */
	dR  = 0.5*path.boxPtr->side[0] / (1.0*NRADSEP);

	/* This is a diagonal estimator that gets its own file */
	initialize(NRADSEP);

	/* Set estimator name */
	name = "Radial Density";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NRADSEP; n++) 
		header.append(str(format("%16.3E") % ((n)*dR)));

	norm = 1.0 / (path.boxPtr->side[NDIM-1]*path.numTimeSlices);
	for (int n = 0; n < NRADSEP; n++) 
		norm(n) /= (M_PI*(2*n+1)*dR*dR);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RadialDensityEstimator::~RadialDensityEstimator() { 
}

/*************************************************************************//**
 *  Accumulate a histogram of all particle distances from the axis.
******************************************************************************/
void RadialDensityEstimator::accumulate() {

	dVec pos;
	beadLocator beadIndex;
	for (int slice = 0; slice < path.numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
			beadIndex = slice,ptcl;
			pos = path(beadIndex);
			int k = int(sqrt(pos[0]*pos[0]+pos[1]*pos[1])/dR);
			if (k < NRADSEP)
				estimator(k) += 1.0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CYLINDER ESTIMATORS ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

/*************************************************************************//**
 * Determine if a position is inside the cutoff radius
******************************************************************************/
inline bool include(const dVec &r, double maxR) {
	return (r[0]*r[0] + r[1]*r[1] < maxR*maxR);
}

/*************************************************************************//**
 * Count the number of particles inside a given radius.
 *
 * Here we arbitrarily only count slice 0 particles.
******************************************************************************/
int num1DParticles(const Path &path, double maxR) {
	int tot = 0;
	dVec r;
	for (int ptcl = 0; ptcl < path.numBeadsAtSlice(0); ptcl++) {
		r = path(0,ptcl);
		if (include(path(0,ptcl),maxR))
			tot++;
	}
	return tot;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER ENERGY ESTIMATOR CLASS -------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 * We measure the total Kinetic, Potential and E-mu*N as well as the kinetic,
 * potential and total energy per particle.
******************************************************************************/
CylinderEnergyEstimator::CylinderEnergyEstimator (const Path &_path, ActionBase *_actionPtr,
		double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label), actionPtr(_actionPtr) {

	/* Assign the cut-off radius */
	maxR = _maxR;

	/* We compute three diagonal estimators, kinetic, potential and total energy
	 * per particle */
	initialize(7);
    endLine = false;

	/* Set estimator name and header, we will always report the energy
	 * first, hence the comment symbol*/
	name = "Cyl Energy";
	header = str(format("#%15s%16s%16s%16s%16s%16s%16s") 
			% "K" % "V" % "E" % "E_mu" % "K/N" % "V/N" % "E/N");
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderEnergyEstimator::~CylinderEnergyEstimator() { 
}

/*************************************************************************//**
 *  Accumluate the energy.
 *
 *  We use the thermodynamic estimator for the kinetic energy and the
 *  potential estimator for the potential energy. A possible shift is 
 *  made due to a tail correction.
******************************************************************************/
void CylinderEnergyEstimator::accumulate() {

	double totK = 0.0;
	double totV = 0.0;
	double gV2 = 0.0;

	int numParticles  = num1DParticles(path,maxR);
	int numTimeSlices = path.numTimeSlices;

	/* The kinetic normalization factor */
	double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * numTimeSlices);

	/* The action correction normalization factor */
	double normAct = 2.0 * constants()->lambda() * constants()->tau() * constants()->tau()
		/ (1.0 * numTimeSlices);

	/* The classical contribution to the kinetic energy per particle 
	 * including the chemical potential */
	double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;

	/* We first calculate the kinetic energy.  Even though there
	 * may be multiple mixing and swaps, it doesn't matter as we always
	 * just advance one time step at a time, as taken care of through the
	 * linking arrays.  This has been checked! */
	beadLocator beadIndex;
	dVec vel;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
			beadIndex = slice,ptcl;
			if (include(path(beadIndex),maxR)) {
				vel = path.getVelocity(beadIndex);
				totK -= dot(vel,vel);
			}
		}
	}

	/* Normalize the accumulated link-action part */
	totK *= kinNorm;

	/* Now we compute the potential and kinetic energy.  We use an operator estimater
	 * for V and the thermodynamic estimator for K */
	int eo;
	for (int slice = 0; slice < numTimeSlices; slice++) {
		eo = (slice % 2);
		if (actionPtr->eFactor[eo] > EPS) 
			gV2  += (actionPtr->eFactor[eo]) * (actionPtr->potentialPtr->gradVSquared(slice,maxR));
		else if (eo==0)
			totV  += actionPtr->potentialPtr->V(slice,maxR);
	}

	/* Normalize the action correction and the total potential*/
	gV2 *= normAct;
	totV /= (0.5 * numTimeSlices);

	/* Perform all the normalizations and compute the individual energy terms */
	totK  += (classicalKinetic + gV2);

	/* Now we accumulate the average total, kinetic and potential energy, 
	 * as well as their values per particles, provided we have at least one
     * particle in the central region. */
    if (numParticles > 0) {
        estimator(0) += totK;
        estimator(1) += totV;
        estimator(2) += totK + totV;

        estimator(3) += totK + totV - constants()->mu()*numParticles;

        estimator(4) += totK/(1.0*numParticles);
        estimator(5) += totV/(1.0*numParticles);
        estimator(6) += (totK + totV)/(1.0*numParticles);
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER NUM PARTICLES ESTIMATOR CLASS ------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  We measure the number of particles, (number of particles)^2 and the
 *  density of particles.
******************************************************************************/
CylinderNumberParticlesEstimator::CylinderNumberParticlesEstimator (const Path &_path, 
		double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

	/* Assign the cut-off radius */
	maxR = _maxR;

	/* We compute three diagonal estimators, the total number of particles,
	 * total number of particles squared and density. */
	initialize(3);

	/* Set estimator name and header */
	name = "Cyl Number Particles";
	header = str(format("%16s%16s%16s") % "N" % "N^2" % "density");
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderNumberParticlesEstimator::~CylinderNumberParticlesEstimator() { 
}

/*************************************************************************//**
 * Accumulate the number of particles and density.
******************************************************************************/
void CylinderNumberParticlesEstimator::accumulate() {
	int numParticles = num1DParticles(path,maxR);
	estimator(0) += 1.0*numParticles;
	estimator(1) += 1.0*numParticles*numParticles;
	estimator(2) += 1.0*numParticles/path.boxPtr->side[NDIM-1];
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER NUMBER DISTRIBUTION ESTIMATOR CLASS ------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  The number of particles probability distribution is setup.  We allow 
 *  for particle fluctuations up to 50 particles away from the initial value
 *  specified.
******************************************************************************/
CylinderNumberDistributionEstimator::CylinderNumberDistributionEstimator 
	(const Path &_path, double _maxR, int _frequency, string _label) : 
        EstimatorBase(_path,_frequency,_label) {

	/* Assign the cut-off radius */
	maxR = _maxR;

	/* For now, we assume a maximum of 200 total particles. */
	maxNumParticles = 200;
	initialize(maxNumParticles);

	/* Set estimator name and header */
	name = "Cyl Number Distribution";
	header = str(format("#%15d") % 0);
	for (int n = 1; n < maxNumParticles; n++) 
		header.append(str(format("%16d") % n));
}

/*************************************************************************//**
 * Destructor.
******************************************************************************/
CylinderNumberDistributionEstimator::~CylinderNumberDistributionEstimator() { 
}

/*************************************************************************//**
 * Accumulate the number probability distribution
******************************************************************************/
void CylinderNumberDistributionEstimator::accumulate() {

	/* Get the correct particle Number index, and increment the corresponding bin */
	int index = num1DParticles(path,maxR);
	if (index >= 0 && index < maxNumParticles)
		estimator(index) += 1.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER SUPERFLUID FRACTION ESTIMATOR CLASS ------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Measure the superfluid fraction from the value of the winding number in
 *  all directions with periodic boundary conditions.  We also measure
 *  the individual values of the winding numbers and their probability
 *  distribution.
******************************************************************************/
CylinderSuperfluidFractionEstimator::CylinderSuperfluidFractionEstimator (const Path &_path, 
		double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_frequency,_label) {

	/* Assign the cut-off radius */
	maxR = _maxR;

	windMax = 10;
	/* We compute a bunch of estimators here, the superfluid fraction, the winding
	 * number in all possible dimensions, and the winding number histograms up to 
	 * windMax windings. These are all diagonal estimators and we have our own
	 * output file.*/
	initialize(4+2*windMax+1);

	/* Set estimator name */
	name = "Cyl Superfluid Fraction";
	header = str(format("#%15s%16s%16s%16s") % "rho_s/rho" % "W^2(x)" % "W^2(y)" % "W^2(z)");
	for (int w = -windMax; w <= windMax; w++)
		header += str(format("%11sP(%+1d)") % " " % w);

	/* The pre-factor for the superfluid density is always the same */
	norm(0) = constants()->T() / (2.0 * sum(path.boxPtr->periodic) * constants()->lambda());
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderSuperfluidFractionEstimator::~CylinderSuperfluidFractionEstimator() { 
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  We measure the winding number W, as well as the W^2/N which is used to 
 *  compute the superfluid fraction.  The probability distribution of winding
 *  number fluctuations is also computed.
******************************************************************************/
void CylinderSuperfluidFractionEstimator::accumulate() {

	double locW2oN = 0.0;

	/* Sum up the winding number over all particles */
	dVec W,locW,vel;
	W = 0.0;
	locW = 0.0;

	/* The start bead for each world line, and the moving index */
	beadLocator startBead;
	beadLocator beadIndex;

	int numWorldlines = path.numBeadsAtSlice(0);

	/* We create a local vector, which determines whether or not we have
	 * already included a bead at slice 0*/
	doBead.resize(numWorldlines);
	doBead = true;

	/* Needed to ensure an included world line */
	bool includeWorldline = true;

	/* We go through each particle/worldline */
	for (int n = 0; n < numWorldlines; n++) {

		/* The initial bead to be moved */
		startBead = 0,n;

		/* We make sure we don't try to touch the same worldline twice */
		if (doBead(n)) {

			/* Mark the beads as touched and increment the number of worldlines */
			beadIndex = startBead;

			/* Go through all worldlines, summing up the winding number */
			locW = 0.0;
			includeWorldline = true;
			do {

				/* We turn off any zero-slice beads we have touched */
				if (beadIndex[0]==0)
					doBead(beadIndex[1]) = false;

				/* If the bead is inside our cutoff radius, include its winding */
				if (include(path(beadIndex),maxR)) {
					vel = path.getVelocity(beadIndex);
					locW += vel;
				}
				else
					includeWorldline = false;

				beadIndex = path.next(beadIndex);
			} while (!all(beadIndex==startBead));

			if (includeWorldline)
				W += locW;

		} // doBead
	} // worldLine
				
	/* Scale by the periodicity of the boundary conditions */
	W *= path.boxPtr->periodic;

	/* Compute the locally scaled W^2/N */
    int numParticles = num1DParticles(path,maxR);
    if (numParticles > 0)
        locW2oN = dot(W,W)/(1.0*numParticles);
    else
        locW2oN = 0.0;

	/* The average winding number squared */
	estimator(0) += locW2oN;

	/* Scale by the length of the system in each dimension*/
	W *= path.boxPtr->sideInv;

	/* The individual winding numbers, we always store 3 values regardless
	 * of the dimensionality to ensure output file consistency */
	int i;
	for (i = 0; i < NDIM; i++)
		estimator(1+i) += W[i]*W[i];
	for (int j = i; j < 3; j++)
		estimator(1+j) += 0.0;

	/* Calcluate the winding number probablity in the NDIM^th dimension */
	int n = 0;
	for (int p = -windMax; p <= windMax; p++) {
		if (abs(W[NDIM-1]-1.0*p) < 0.2)
			estimator(4+n) += 1.0;
		++n;
	}
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER ONE BODY DENSITY MATRIX ESTIMATOR CLASS --------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  The one body density matrix estimator is initialized.  We measure NOBDMSEP
 *  positions, out to the maximum separation in the sample (which may depend
 *  on the type of simulation cell).
******************************************************************************/
CylinderOneBodyDensityMatrixEstimator::CylinderOneBodyDensityMatrixEstimator (Path &_path,
		ActionBase *_actionPtr, MTRand &_random, double _maxR, int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label), 
	lpath(_path),
	actionPtr(_actionPtr),
	random(_random)
{
	/* Assign the cut-off radius */
	maxR = _maxR;

	sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());

	/* We chooose the maximum separation to be sqrt(NDIM)*L/2 */
	dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NOBDMSEP);

	/* This is an off-diagonal estimator that gets its own file */
	initialize(NOBDMSEP);
    diagonal = false;

	/* Set estimator name */
	name = "Cyl One Body Density Matrix";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NOBDMSEP; n++) 
		header.append(str(format("%16.3E") % (n*dR)));

	numReps = 10;
	norm = 1.0 / (1.0*numReps);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderOneBodyDensityMatrixEstimator::~CylinderOneBodyDensityMatrixEstimator() { 
}

/*************************************************************************//**
 *  Sample the OBDM.
 * 
 *  We overload the sample method for the one body density matrix as we
 *  only want to measure when the gap is not too large, otherwise we will be
 *  dominated by tiny close probabilities.
******************************************************************************/
void CylinderOneBodyDensityMatrixEstimator::sample() {
	numSampled++;
	if ( frequency && 
         ((numSampled % frequency) == 0) && 
         (path.worm.isConfigDiagonal == diagonal) &&
         (path.worm.gap > 0) && (path.worm.gap <= constants()->Mbar())  &&
         (actionPtr->eFactor[(lpath.worm.tail[0] % 2)] < EPS) ) {
		
		/* We only attempt the partial-close if both the head and tail
		 * are within the include region. */
		if ( include(path(path.worm.head),maxR) && include(path(path.worm.tail),maxR) ) {
			totNumAccumulated++;
			numAccumulated++;
			accumulate();
		}
	}
}

/*************************************************************************//**
 *  Return a dimensionally dependent random vector of length r.  
 *
 *  If we are in 3D in a cylinder geometry, we only shift in the z-direction.
 *  @param r The length of the random vector
 *  @return a random NDIM-vector of length r
******************************************************************************/
inline dVec CylinderOneBodyDensityMatrixEstimator::getRandomVector(const double r) {
	dVec rVec;
	rVec = 0.0;
	if (random.rand() < 0.5)
		rVec[NDIM-1] = r;
	else
		rVec[NDIM-1] = -r;

	return rVec;
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
dVec CylinderOneBodyDensityMatrixEstimator::newStagingPosition(const beadLocator &neighborIndex, 
		const beadLocator &endIndex, const int stageLength, const int k) {

    /* The rescaled value of lambda used for staging */
    double f1 = 1.0 * (stageLength - k - 1);
    double f2 = 1.0 / (1.0*(stageLength - k));
    double sqrtLambdaKTau = sqrt2LambdaTau * sqrt(f1 * f2);

	/* We find the new 'midpoint' position which exactly samples the kinetic 
	 * density matrix */
	neighborPos = lpath(neighborIndex);
	newRanPos = lpath(endIndex) - neighborPos;
	lpath.boxPtr->putInBC(newRanPos);
	newRanPos *= f2;
	newRanPos += neighborPos;

	/* This is the random kick around that midpoint */
	for (int i = 0; i < NDIM; i++)
		newRanPos[i] = random.randNorm(newRanPos[i],sqrtLambdaKTau);

    lpath.boxPtr->putInside(newRanPos);

    return newRanPos;
}

/*************************************************************************//**
 *  Accumulate the OBDM.
 * 
 *  We perform a fake close move, where the head of the worm is advanced to
 *  a position a distance 'r' away from the tail but at the same time slice.
 *  The probability of excepting such a move is equal (up to normalization)
 *  to the one body density matrix.
******************************************************************************/
void CylinderOneBodyDensityMatrixEstimator::accumulate() {

	oldTailPos = lpath(lpath.worm.tail);
	oldAction = actionPtr->potentialAction(lpath.worm.tail);

	/* We make a list of all the beads involved in the move, adding them
	 * as we go. */
	beadLocator beadIndex;
	beadIndex = lpath.worm.head;
	dVec pos;
	pos = 0.0;
	for (int k = 0; k < (lpath.worm.gap-1); k++) 
		beadIndex = lpath.addNextBead(beadIndex,pos);

	/* Perform the final connection to the tail*/
	lpath.next(beadIndex) = lpath.worm.tail;
	lpath.prev(lpath.worm.tail) = beadIndex;

	for (int p = 0; p < numReps; p++) {

		/* Now we loop through all possible separations, evaluating the potential
		 * action */
		for (int n = 0; n < NOBDMSEP; n++) {

			newAction = 0.0;
			++numAttempted;

			/* Assign the new displaced tail position */
			newTailPos = oldTailPos + getRandomVector(n*dR);
			lpath.boxPtr->putInside(newTailPos);
			lpath.updateBead(lpath.worm.tail,newTailPos);

			/* Compute the free particle density matrix */
			rho0Norm = actionPtr->rho0(lpath.worm.head,lpath.worm.tail,lpath.worm.gap);

			/* action shift coming from a finite chemical potential */
			double muShift = lpath.worm.gap*constants()->mu()*constants()->tau();

			/* Starting from the head, we generate new positions for the beads via
			 * staging, and accumulate the potential action */
			beadIndex = lpath.worm.head;
			int k = 0;
			do {
				if (!all(beadIndex==lpath.worm.head) && !all(beadIndex==lpath.worm.tail)) {
					lpath.updateBead(beadIndex,
							newStagingPosition(path.prev(beadIndex),lpath.worm.tail,lpath.worm.gap,k));
					++k;
				}

				newAction += actionPtr->potentialAction(beadIndex);

				beadIndex = lpath.next(beadIndex);
			} while (!all(beadIndex==lpath.next(lpath.worm.tail)));

			double expAction = exp(-newAction + oldAction + muShift);
			estimator(n) += rho0Norm*expAction;

			/* Record the probability of accepting the move */
			if (random.randExc() < expAction)
				++numAccepted;

		} // end for n

	} // end for k

	/* Now we must undo any damge we have caused by reverting the tail to its previous position,
	 * and turning off all intermediate beads */
	lpath.updateBead(lpath.worm.tail,oldTailPos);

	/* Delete all the beads that were added. */
	beadIndex = lpath.next(lpath.worm.head);
	while (!all(beadIndex==lpath.worm.tail)) {
		beadIndex = lpath.delBeadGetNext(beadIndex);
	}
	lpath.next(lpath.worm.head) = XXX;
	lpath.prev(lpath.worm.tail) = XXX;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER PAIR CORRELATION ESTIMATOR CLASS ---------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  For the pair correlation function, we measure NPCFSEP positions to get
 *  high enough data density to observe possible osscilations.  The normalization
 *  depends on dimension.
******************************************************************************/
CylinderPairCorrelationEstimator::CylinderPairCorrelationEstimator (const Path &_path, 
		ActionBase *_actionPtr, double _maxR, int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label),
	actionPtr(_actionPtr)
{
	/* Assign the cut-off radius */
	maxR = _maxR;

	/* The spatial discretization */
	dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

	/* This is a diagonal estimator that gets its own file */
	initialize(NPCFSEP);

	/* Set estimator name */
	name = "Cyl Pair Correlation Function";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NPCFSEP; n++) 
		header.append(str(format("%16.3E") % ((n)*dR)));

	/* The normalization factor for the pair correlation function */
	norm = 0.5*path.boxPtr->side[NDIM-1] / dR;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderPairCorrelationEstimator::~CylinderPairCorrelationEstimator() { 
}

/*************************************************************************//**
 *  We simply update the local value of the esitimator with the value from
 *  our potential reference.  
 *
 *  We only compute this for N1D > 1.
******************************************************************************/
void CylinderPairCorrelationEstimator::accumulate() {
	int N1D = num1DParticles(path,maxR);
    if (N1D > 1) {
        double lnorm = 1.0*sum(actionPtr->potentialPtr->cylSepHist);
        lnorm /= 1.0*(N1D-1)/(1.0*N1D);
        estimator += 1.0*actionPtr->potentialPtr->cylSepHist / (1.0*lnorm);
    }
}

/**************************************************************************//**
 *  Sample the estimator.
 * 
 *  Here we overload the cylinder pair correlation function estimator, as
 *  we only measure when we have some relevant particle separations.
******************************************************************************/
void CylinderPairCorrelationEstimator::sample() {
    numSampled++;

    if (baseSample() && (sum(actionPtr->potentialPtr->cylSepHist) > 0)) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER RADIAL POTENTIAL ESTIMATOR CLASS ---------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  We compute the effective radial potential by averaging the interaction + 
 *  external potential felt by the central chain of particles.
******************************************************************************/
CylinderRadialPotentialEstimator::CylinderRadialPotentialEstimator (const Path &_path, 
		ActionBase *_actionPtr, MTRand &_random, double _maxR, int _frequency, string _label) : 
	EstimatorBase(_path,_frequency,_label),
	actionPtr(_actionPtr),
	random(_random)
{
	/* Assign the cut-off radius */
	maxR = _maxR;

	/* The spatial discretization */
	dR  = 0.5*path.boxPtr->side[0] / (1.0*NRADSEP);

	/* This is a diagonal estimator that gets its own file */
	initialize(NRADSEP);
	radPot.resize(NRADSEP);

	/* Set estimator name */
	name = "Cyl Radial Potential";

	/* The header is the first line which contains the spatial separations */
	header = str(format("#%15.3E") % 0.0);
	for (int n = 1; n < NRADSEP; n++) 
		header.append(str(format("%16.3E") % ((n)*dR)));
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderRadialPotentialEstimator::~CylinderRadialPotentialEstimator() { 
}

/*************************************************************************//**
 *  We determine what the effective radial potential is inside the cutoff 
 *  distance.
******************************************************************************/
void CylinderRadialPotentialEstimator::accumulate() {

	double totV = 0.0;
	dVec r1,r2; 		// The two bead positions
	dVec sep;			// The bead separation

	beadLocator bead2;	// The bead locator

	/* Choose a random position */
	for (int n = 0; n < NRADSEP; n++) {

		totV = 0.0;

		double theta = 2.0*M_PI*random.rand();
		r1[0] = n*dR*cos(theta);
		r1[1] = n*dR*sin(theta);
		r1[2] = path.boxPtr->side[2]*(-0.5 + random.randExc());

		/* We sum up the external and interaction energy over all slices*/
		for (int slice = 0; slice < path.numTimeSlices; slice++) {
			bead2[0] = slice;

			int numParticles = path.numBeadsAtSlice(slice);

			/* Sum over particles */
			for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

				r2 = path(bead2);
				if (!include(r2,maxR)) {
					sep = r2 - r1;
					path.boxPtr->putInBC(sep);
					totV += actionPtr->potentialPtr->interactionPtr->V(sep);
				} // bead2 is outside maxR
			} // bead2
		} // slice

		totV /= 1.0*path.numTimeSlices;
		totV += actionPtr->potentialPtr->externalPtr->V(r1);

		estimator(n) += totV;
	} // n
}

/*************************************************************************//**
 *  We compute the total potential for all particles inside the cutoff as a
 *  function of their position.
 *
 *  We only compute this for N > 1.
******************************************************************************/
void CylinderRadialPotentialEstimator::accumulate1() {

	double totV = 0.0;
	dVec r1,r2; 		// The two bead positions
	dVec sep;			// The bead separation
	double rad1,rad2;	// The two bead radii
	int nR;				// The bin number

	beadLocator bead1,bead2;	// The bead locators
	bool found1,found2;			// Are the beads in the central chain
	found1 = found2 = false;
	radPot = 0.0;
	int numFound1 = 0;

	/* We sum up the external and interaction energy over all slices*/
	for (int slice = 0; slice < path.numTimeSlices; slice++) {
		bead1[0] = bead2[0] = slice;

		int numParticles = path.numBeadsAtSlice(slice);

		/* Sum over particles */
		for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {

			r1 = path(bead1);
			rad1   = r1[0]*r1[0] + r1[1]*r1[1];
			found1 = (rad1 < maxR*maxR);

			/* If the first particle is in the central chain, looks for its
			 * interacting partners */
			if (found1) {
				totV = actionPtr->potentialPtr->externalPtr->V(r1);
				numFound1 ++;

				/* We don't have to worry about double counting here, as
				 * we never allow two particles in the central chain to 
				 * interact */
				for (bead2[1] = 0; bead2[1] < numParticles; bead2[1]++) {

					r2 = path(bead2);
					rad2   = r2[0]*r2[0] + r2[1]*r2[1];
					found2 = (rad2 < maxR*maxR);

					/* Only accumulate the energy if bead1 is in the central
					 * chain while bead2 is not */
					if (!found2) {
						sep = path.getSeparation(bead2,bead1);
						totV += actionPtr->potentialPtr->interactionPtr->V(sep);
					} // !found2
				} // bead2
				
				nR = int(sqrt(rad1)/dR);
				if (nR < NRADSEP)
					radPot(nR) += totV;
			} // found1
		} // bead1
	} // slice

	radPot /= (1.0*numFound1);
	estimator += radPot;
}
