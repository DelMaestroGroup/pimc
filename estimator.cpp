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
#include "factory.h"

bool init(vector<string>& names, string name) {
    names.push_back(name);
    return true;
}

/**************************************************************************//**
 * Setup the estimator factory.
******************************************************************************/
EstimatorFactory estimatorFactory;
#define REGISTER_ESTIMATOR(NAME,TYPE) \
    const string TYPE::name = NAME;\
    bool reg ## TYPE = estimatorFactory()->Register<TYPE>(TYPE::name);\
    bool init ## TYPE = init(EstimatorFactory::names,TYPE::name);

/**************************************************************************//**
 * Estimator naming conventions:
 *
 * 1) be as descriptive as possible
 * 2) prepend the keyword "cylinder" for cylinder estimators
 * 3) prepend the keyword "pigs" for PIGS estimators
 * 4) include the word "multi" for multi-path estimators 
******************************************************************************/
REGISTER_ESTIMATOR("energy",EnergyEstimator);
REGISTER_ESTIMATOR("virial",VirialEnergyEstimator);
REGISTER_ESTIMATOR("number particles",NumberParticlesEstimator);
REGISTER_ESTIMATOR("number distribution",NumberDistributionEstimator);
REGISTER_ESTIMATOR("time",TimeEstimator);
REGISTER_ESTIMATOR("particle position",ParticlePositionEstimator);
REGISTER_ESTIMATOR("bipartition density",BipartitionDensityEstimator);
REGISTER_ESTIMATOR("planar density rho",PlaneParticlePositionEstimator);
REGISTER_ESTIMATOR("superfluid fraction",SuperfluidFractionEstimator);
REGISTER_ESTIMATOR("planar winding rhos/rho",PlaneWindingSuperfluidDensityEstimator);
REGISTER_ESTIMATOR("planar area rhos/rho",PlaneAreaSuperfluidDensityEstimator);
REGISTER_ESTIMATOR("radial winding rhos/rho",RadialWindingSuperfluidDensityEstimator);
REGISTER_ESTIMATOR("radial area rhos/rho",RadialAreaSuperfluidDensityEstimator);
REGISTER_ESTIMATOR("local superfluid",LocalSuperfluidDensityEstimator);
REGISTER_ESTIMATOR("diagonal fraction",DiagonalFractionEstimator);
REGISTER_ESTIMATOR("worm properties",WormPropertiesEstimator);
REGISTER_ESTIMATOR("permutation cycle",PermutationCycleEstimator);
REGISTER_ESTIMATOR("local permutation",LocalPermutationEstimator);
REGISTER_ESTIMATOR("one body density matrix",OneBodyDensityMatrixEstimator);
REGISTER_ESTIMATOR("pair correlation function",PairCorrelationEstimator);
REGISTER_ESTIMATOR("static structure factor",StaticStructureFactorEstimator);
REGISTER_ESTIMATOR("intermediate scattering function",IntermediateScatteringFunctionEstimator);
REGISTER_ESTIMATOR("radial density",RadialDensityEstimator);
REGISTER_ESTIMATOR("cylinder energy",CylinderEnergyEstimator);
REGISTER_ESTIMATOR("cylinder number particles",CylinderNumberParticlesEstimator);
REGISTER_ESTIMATOR("cylinder number distribution",CylinderNumberDistributionEstimator);
REGISTER_ESTIMATOR("cylinder linear density",CylinderLinearDensityEstimator);
REGISTER_ESTIMATOR("cylinder superfluid fraction",CylinderSuperfluidFractionEstimator);
REGISTER_ESTIMATOR("cylinder one body density matrix",CylinderOneBodyDensityMatrixEstimator);
REGISTER_ESTIMATOR("cylinder pair correlation function",CylinderPairCorrelationEstimator);
REGISTER_ESTIMATOR("cylinder radial potential",CylinderRadialPotentialEstimator);
REGISTER_ESTIMATOR("cylinder linear potential",CylinderLinearPotentialEstimator);
REGISTER_ESTIMATOR("cylinder potential energy",PotentialEnergyEstimator);
REGISTER_ESTIMATOR("pigs kinetic energy",KineticEnergyEstimator);
REGISTER_ESTIMATOR("pigs energy",PigsEnergyEstimator);
REGISTER_ESTIMATOR("pigs thermodynamic energy",PigsThermoEnergyEstimator);
REGISTER_ESTIMATOR("pigs total energy",TotalEnergyEstimator);
REGISTER_ESTIMATOR("pigs thermodynamic potential energy",ThermoPotentialEnergyEstimator);
REGISTER_ESTIMATOR("pigs positions",PositionEstimator);
REGISTER_ESTIMATOR("pigs particle resolved positions",ParticleResolvedPositionEstimator);
REGISTER_ESTIMATOR("pigs particle correlations",ParticleCorrelationEstimator);
REGISTER_ESTIMATOR("pigs velocity",VelocityEstimator);
REGISTER_ESTIMATOR("pigs subregion occupation",SubregionOccupationEstimator);
REGISTER_ESTIMATOR("pigs one body density matrix",PIGSOneBodyDensityMatrixEstimator);

/**************************************************************************//**
 * Setup the estimator factory for multi path estimators.
******************************************************************************/
/* Factory<EstimatorBase* (Path &, Path&, ActionBase *, ActionBase*, MTRand &, double)> MultiEstimatorFactory; */
MultiEstimatorFactory multiEstimatorFactory;
const string SwapEstimator::name = "pigs multi swap";
bool regSwap = multiEstimatorFactory()->Register<SwapEstimator>(SwapEstimator::name);
bool initSwap = init(MultiEstimatorFactory::names,SwapEstimator::name);

const string EntPartEstimator::name = "pigs multi entanglement of particles";
bool regEntPart = multiEstimatorFactory()->Register<EntPartEstimator>(EntPartEstimator::name);
bool initEntPart = init(MultiEstimatorFactory::names,EntPartEstimator::name);

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ESTIMATOR BASE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 * Initialize and specify the output file.
 *
 * We take in some options that might only be needed by derived classes
 * (_actionPtr, _random, _maxR) in order to facilitate instantiation via a
 * factory.
******************************************************************************/
EstimatorBase::EstimatorBase(const Path &_path, ActionBase *_actionPtr, 
        const MTRand &_random, double _maxR, int _frequency, string _label) :
    path(_path),
    actionPtr(_actionPtr),
    random(_random),
    maxR(_maxR),
    frequency(_frequency),
    label(_label),
    numSampled(0),
    numAccumulated(0),
    totNumAccumulated(0),
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
    /* cout << getName() << "\t" << frequency << endl; */
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

/*************************************************************************//**
*  AppendLabel
******************************************************************************/
void EstimatorBase::appendLabel(string append) {
    label = label + append;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// TIME ESTIMATOR CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 * Used to measure the time between bins
******************************************************************************/
TimeEstimator::TimeEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Set estimator name and header */
    header = str(format("%16s%16s") % "us" % "mcsteps");
    endLine = false;
    initialize(2);
}

/*************************************************************************//**
 *  Overload sampling to make sure it is always done, regardless of ensemble.
******************************************************************************/
void TimeEstimator::sample() {

    numSampled++;

    if (frequency && ((numSampled % frequency) == 0)) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
}

/*************************************************************************//**
 *  Grab the initial time.
******************************************************************************/
void TimeEstimator::accumulate() {
    if (totNumAccumulated == 1)
        time_begin = std::chrono::high_resolution_clock::now();
}

/*************************************************************************//**
 *  Grab the final time and write to disk.  
******************************************************************************/
void TimeEstimator::output() {
    time_end = std::chrono::high_resolution_clock::now();
    estimator(0) = 0.001*std::chrono::duration_cast<std::chrono::nanoseconds>(
            time_end - time_begin).count()/numAccumulated;
    estimator(1) = 1.0*numAccumulated;
        
    // close the loop
    time_begin = time_end;

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
 * potential and total energy per particle using the thermodynamic and operator
 * estimators.
******************************************************************************/
EnergyEstimator::EnergyEstimator (const Path &_path, ActionBase *_actionPtr,
        const MTRand &_random, double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Set estimator name and header, we will always report the energy
     * first, hence the comment symbol*/
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
    double totVop = 0.0;
    double totV = 0.0;

    int numParticles  = path.getTrueNumParticles();
    int numTimeSlices = path.numTimeSlices;

    /* The total tail correction */
    double tailV = (1.0*numParticles*numParticles/path.boxPtr->volume)
        * actionPtr->interactionPtr->tailV;

    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * numTimeSlices);

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
    double t1 = 0.0;
    double t2 = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        eo = (slice % 2);
        t1 += actionPtr->derivPotentialActionLambda(slice);
        t2 += actionPtr->derivPotentialActionTau(slice);
        if (eo==0) 
            totVop  += actionPtr->potential(slice);
    }

    t1 *= constants()->lambda()/(constants()->tau()*numTimeSlices);
    t2 /= 1.0*numTimeSlices;

    /* Normalize the action correction and the total potential*/
    totVop /= (0.5 * numTimeSlices);

    /* Perform all the normalizations and compute the individual energy terms */
    totK += (classicalKinetic + t1);
    totV = t2 - t1 + tailV;

    totVop += tailV;

    totV = totVop;

    /* Now we accumulate the average total, kinetic and potential energy, 
     * as well as their values per particles. */
    estimator(0) += totK;
    estimator(1) += totV;
    estimator(2) += totK + totV;

    estimator(3) += totK + totV - constants()->mu()*numParticles;

    if (numParticles > 0) {
        estimator(4) += totK/(1.0*numParticles);
        estimator(5) += totV/(1.0*numParticles);
        estimator(6) += (totK + totV)/(1.0*numParticles);
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// VIRIAL ENERGY ESTIMATOR CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 * We measure the total Kinetic, Potential and E-mu*N as well as the kinetic,
 * potential and total energy per particle.
 *
 *
******************************************************************************/
VirialEnergyEstimator::VirialEnergyEstimator (const Path &_path, ActionBase *_actionPtr,
        const MTRand &_random, double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Set estimator name and header, we will always report the energy
     * first, hence the comment symbol*/
    header = str(format("#%15s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s") 
            % "K_op" % "K_cv" % "V_op" % "V_cv" % "E" % "E_mu" % "K_op/N" % "K_cv/N" % "V_op/N"
            % " V_cv/N" % "E/N" % "EEcv*Beta^2"% "Ecv*Beta" % "dEdB" % "CvCov1"
            % "CvCov2" % "CvCov3" % "E_th" % "P");
    endLine = false;
    initialize(19);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
VirialEnergyEstimator::~VirialEnergyEstimator() { 
}

/*************************************************************************//**
 *  Accumluate the energy.
 *
 * We use the thermodynamic estimator for the kinetic energy and the
 * potential estimator for the potential energy. A possible shift is 
 * made due to a tail correction.
 *
 * We measure the Centroid Virial energy estimator and the potential energy
 * estimator via the operator method.  The potential is subtracted off
 * from the total energy to get the kinetic energy.  A tail correction term
 * is included for the potential energy due to periodic BC.
 *
 * We also measure the centroid virial specific heat terms so that the 
 * specific heat may be computed post-process.
 *
******************************************************************************/
void VirialEnergyEstimator::accumulate() {

    /* Set up potential operator energy, thermodynamic energy,
     * centroid virial energy, centroid virial kinetic energy (see Jang Jang Voth),
     * and terms that go into both energy estimators. */
    double totVop = 0.0;
    double thermE = 0.0; // = thermTerm1 +  + T5 + tailV
    double totEcv = 0.0; // = T1+T2+T3+T4+T5+tailV
    double Kcv = 0.0; // = T1+T2+T3+T4+virKinTerm
    double T2 = 0.0;
    double T3 = 0.0;
    double T4 = 0.0;
    double virKinTerm = 0.0;

    int numParticles  = path.getTrueNumParticles();
    int numTimeSlices = path.numTimeSlices;
    double beta = 1.0*numTimeSlices*constants()->tau();
    int virialWindow = constants()->virialWindow();

    /* The total tail correction */
    double tailV = (1.0*numParticles*numParticles/path.boxPtr->volume)
        * actionPtr->interactionPtr->tailV;

    /* The constant term from the thermodynamic energy. */
    double thermTerm1 = (0.5 * NDIM / constants()->tau()) * numParticles;

    /* The constant term from the centroid virial energy. */ 
    double T1 = 0.5 * NDIM * numParticles / (1.0*virialWindow*constants()->tau());

    /* Calculate the exchange energy for centroid virial energy.
     * Also computes kinetic piece of thermodynamic energy which 
     * fluctuates wildly. */
    double exchangeNorm = 1.0/(4.0*virialWindow*pow(constants()->tau(),2)
            *constants()->lambda()*numTimeSlices);

    /* Compute the thermodynamic pressure */
    double Pressure = NDIM*numParticles;
    double P2,P3;
    P2 = P3 = 0.0;

    beadLocator bead1, beadNext, beadNextOld;
    dVec vel1, vel2;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            bead1 = slice,ptcl; // current bead
            vel2 = path.getVelocity(bead1);
            beadNextOld = bead1;
            vel1 = 0.0;
            /* get r_{current + window} - r_{current} */
            for (int gamma = 1; gamma <= virialWindow; gamma++) {
                beadNext = path.next(bead1, gamma);
                vel1 += path.getSeparation(beadNext, beadNextOld);
                beadNextOld = beadNext;
            }
            T2 -= dot(vel1,vel2);

            /* compute second term of thermodynamic estimator */
            thermE -= dot(vel2,vel2);
        }
    }   
    P2 = thermE;
    T2 *= exchangeNorm;

    /* add second pressure term. */
    P2 /= (2.0*constants()->lambda()*constants()->tau()*numTimeSlices);
    Pressure += P2;

    /* Compute T3, T4, T5, operator potential energy, 
     * and the kinetic energy correction. */
    int eo;
    double T5 = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        eo = (slice % 2);
        T3 += actionPtr->deltaDOTgradUterm1(slice);
        T4 += actionPtr->deltaDOTgradUterm2(slice);
        T5 += actionPtr->derivPotentialActionTau(slice);
        virKinTerm += actionPtr->virKinCorr(slice);
        if (eo==0) 
            totVop  += actionPtr->potential(slice);

        P3 += actionPtr->rDOTgradUterm1(slice)
            + actionPtr->rDOTgradUterm2(slice);
    }

    P3 *= (1.0/(2.0*numTimeSlices));
    Pressure -= P3;

    /* end pressure calculation */
    Pressure /= (NDIM*constants()->tau()*constants()->V());

    totVop /= (0.5 * numTimeSlices);
    totVop += tailV;

    T3 /= (2.0*beta);
    T4 /= (1.0*beta);
    T5 /= (1.0*numTimeSlices);
    virKinTerm /= (0.5*beta);

    /* Normalize second term of the thermodynamic energy then
     * add the rest of the terms to it. */
    thermE *= constants()->fourLambdaTauInv() / (constants()->tau() * numTimeSlices);
    thermE += thermTerm1;
    thermE += T5;
    thermE += tailV;
 
    /* Compute total centroid virial energy */
    totEcv = T1 + T2 + T3 + T4 + T5 + tailV;

    /* Compute centroid virial kinetic energy */
    Kcv = T1 + T2 + T3 + T4 + virKinTerm;

    /* Compute dE/d(\beta) for centroid virial specific heat:
     * C_V^{CV}/(k_B \beta^2) = <E^2> - <E>^2 - <dEdB> */
    double dEdB = (-1.0*T1 - 2.0*T2 + 2.0*T4)/constants()->tau(); // checked

    for (int slice = 0; slice<numTimeSlices; slice++){
        dEdB += actionPtr->secondderivPotentialActionTau(slice)/(1.0*numTimeSlices);
    }
    dEdB *= beta*beta/(1.0*numTimeSlices);

    /* accumulate all energy estimators. */
    estimator(0) += totEcv - totVop; // operator kinetic energy
    estimator(1) += Kcv; // centroid virial kinetic energy
    estimator(2) += totVop; // operator potential energy
    estimator(3) += totEcv - Kcv; //centroid virial potential energy
    estimator(4) += totEcv; // total energy
    estimator(5) += totEcv - constants()->mu()*numParticles;
    estimator(6) += (totEcv - totVop)/(1.0*numParticles);
    estimator(7) += Kcv/(1.0*numParticles);
    estimator(8) += totVop/(1.0*numParticles);
    estimator(9) += (totEcv - Kcv)/(1.0*numParticles);
    estimator(10) += totEcv/(1.0*numParticles);

    /* accumulate specific heat estimators. */
    estimator(11) += totEcv*thermE*beta*beta;
    estimator(12) += totEcv*beta;
    estimator(13) += dEdB;
    estimator(14) += totEcv*thermE*beta*beta*totEcv*beta;
    estimator(15) += totEcv*beta*dEdB;
    estimator(16) += totEcv*thermE*beta*beta*dEdB;

    /* thermodynamic energy */
    estimator(17) += thermE;

    /* Pressure */
    estimator(18) += Pressure;

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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Set estimator name and header */
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label)  : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    initialize(path.boxPtr->numGrid);

    /* Set estimator name and header. */
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
// BIPARTITION DENSITY ESTIMATOR CLASS ---------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Measure the density inside of the constrained region (film)
 *  as well as the unconstrained region (3D bulk) for the excluded volume
 *  (Gasparini) geometry.
******************************************************************************/
BipartitionDensityEstimator::BipartitionDensityEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Set estimator name and header*/
    header = str(format("#%15s%16s") % "film dens" % "bulk dens");
    endLine = true;
    initialize(2);
    norm(0) = 1.0/(path.numTimeSlices);
    norm(1) = 1.0/(path.numTimeSlices);
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
BipartitionDensityEstimator::~BipartitionDensityEstimator() { 
    // empty destructor
}

/*************************************************************************//**
 *  Accumulate density properties from both regions.
 *
 *  We measure the density in both regions of interest.
******************************************************************************/
void BipartitionDensityEstimator::accumulate() {

    /* label the lengths of the sides of the simulation cell */
    dVec lside;
    for (int i = 0; i < NDIM; i++)
        lside[i] = path.boxPtr->side[i];

    /* read in the exclusion lengths */
    Array<double,1> excLens (actionPtr->externalPtr->getExcLen());
    double excZ = excLens(1);
    
    /* determine volume of film region and bulk region */
    double bulkVol ((lside[2]-2.0*excZ)*lside[0]*lside[1]);
    double filmArea (lside[0]*2.0*excZ);

    /* keep track of number of particles in bulk and film */
    int filmNum (0);
    int bulkNum (0);
    
    dVec pos;
    beadLocator beadIndex;
    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;

            pos = path(beadIndex);
            if (pos[2] > excZ)
                bulkNum += 1;
            else if (pos[2] < -excZ)
                bulkNum += 1;
            else
                filmNum += 1;
        }
    }
    /* update estimator with density in both regions*/
    estimator(0) += 1.0*filmNum/(1.0*filmArea);
    estimator(1) += 1.0*bulkNum/(1.0*bulkVol);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PLANE PARTICLE DENSITY ESTIMATOR CLASS ------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  A full NDIM-dimensional particle density hisgram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
PlaneParticlePositionEstimator::PlaneParticlePositionEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    numGrid = (2*NGRIDSEP)*(2*NGRIDSEP);

    /* The spatial discretization */
    for (int i = 0; i < NDIM; i++)
        dl[i]  = path.boxPtr->side[i] / (2.0*NGRIDSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % (1.0*n)));

    /* Compute the area of a grid box */
    double A = 1.0;
    for (int i = 0; i < NDIM-1; i++)
        A *= dl[i];

    norm = 1.0/(1.0*path.numTimeSlices*A*path.boxPtr->side[NDIM-1]);
    side = path.boxPtr->side;
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
PlaneParticlePositionEstimator::~PlaneParticlePositionEstimator() { 
}


/*************************************************************************//**
 *  Accumulate a histogram of all particle positions, with output 
 *  being the running average of the density per grid space.
******************************************************************************/
void PlaneParticlePositionEstimator::accumulate() {

    beadLocator beadIndex;
    dVec pos;

    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            int index = 0;
            for (int i = 0; i < NDIM-1; i++) {  
                int scale = 1;
                for (int j = i+1; j < NDIM-1; j++) 
                    scale *= 2*NGRIDSEP;
                index += scale*static_cast<int>(abs(pos[i] + 0.5*side[i] - EPS ) / (dl[i] + EPS));
            }

            /* update our particle position histogram */
            if (index < numGrid)
                estimator(index) += 1.0;
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    windMax = 10;
    /* We compute a bunch of estimators here, the superfluid fraction, the winding
     * number in all possible dimensions, and the winding number histograms up to 
     * windMax windings. These are all diagonal estimators and we have our own
     * output file.*/
    initialize(4 + 2*windMax + 1 + 1);

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
// PLANE WINDING SUPERFLUID DENSITY ESTIMATOR CLASS --------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Initialize the variables needed to measure the radially averaged winding 
 *  number superfluid density.
 *  @see E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).
 *
******************************************************************************/
PlaneWindingSuperfluidDensityEstimator::PlaneWindingSuperfluidDensityEstimator
(const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    numGrid = (2*NGRIDSEP)*(2*NGRIDSEP);

    /* The spatial discretization */
    dx  = path.boxPtr->side[0] / (2.0*NGRIDSEP);
    dy  = path.boxPtr->side[1] / (2.0*NGRIDSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % (1.0*n)));

    norm = 0.5 * constants()->T()/(dx*dy*path.boxPtr->side[NDIM-1]*constants()->lambda());

    /* Initialize the local arrays */
    locWz.resize(numGrid);
    locWz = 0.0;

    side = path.boxPtr->side;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
PlaneWindingSuperfluidDensityEstimator::~PlaneWindingSuperfluidDensityEstimator() { 
    locWz.free();
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  
******************************************************************************/
void PlaneWindingSuperfluidDensityEstimator::accumulate() {

    int numTimeSlices = path.numTimeSlices;

    beadLocator beadIndex;
    double Wz;
    dVec pos1,pos2;
    dVec vel;

    Wz = 0.0;
    locWz = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
            pos1 = path(beadIndex);
            pos2 = path(path.next(beadIndex));

            int i = static_cast<int>(abs(pos1[0] + 0.5*side[0] - EPS ) / (dx + EPS));
            int j = static_cast<int>(abs(pos1[1] + 0.5*side[1] - EPS ) / (dy + EPS));
            int k = 2*NGRIDSEP*j + i;

            /* The winding number estimator */
            vel = path.getVelocity(beadIndex)*path.boxPtr->periodic;
            Wz += vel[NDIM-1];

            /* The local part of the winding number */
            if (k < numGrid)
                locWz(k) += vel[NDIM-1];
        }
    }

    estimator += locWz*Wz;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PLANE AREA SUPERFLUID DENSITY ESTIMATOR CLASS -----------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Initialize the variables needed to measure the radially averaged winding 
 *  number superfluid density.
 *  @see E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).
 *
******************************************************************************/
PlaneAreaSuperfluidDensityEstimator::PlaneAreaSuperfluidDensityEstimator
(const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    numGrid = (2*NGRIDSEP)*(2*NGRIDSEP);

    /* The spatial discretization */
    dx  = path.boxPtr->side[0] / (2.0*NGRIDSEP);
    dy  = path.boxPtr->side[1] / (2.0*NGRIDSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % (1.0*n)));

    norm = 0.5 * constants()->T()/(dx*dy*path.boxPtr->side[NDIM-1]*constants()->lambda());

    /* Initialize the local arrays */
    locAz.resize(numGrid);
    locAz = 0.0;

    side = path.boxPtr->side;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
PlaneAreaSuperfluidDensityEstimator::~PlaneAreaSuperfluidDensityEstimator() { 
    locAz.free();
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  
******************************************************************************/
void PlaneAreaSuperfluidDensityEstimator::accumulate() {

    int numTimeSlices = path.numTimeSlices;

    beadLocator beadIndex;
    double tAz,Az,rp2;
    dVec pos1,pos2;

    Az = 0.0;
    locAz = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
            pos1 = path(beadIndex);
            pos2 = path(path.next(beadIndex));

            int i = static_cast<int>(abs(pos1[0] + 0.5*side[0] - EPS ) / (dx + EPS));
            int j = static_cast<int>(abs(pos1[1] + 0.5*side[1] - EPS ) / (dy + EPS));
            int k = 2*NGRIDSEP*j + i;

            /*  The distance from the z-axis squared */
            rp2 = pos1[0]*pos1[0] + pos1[1]*pos1[1];
            if (rp2 < dx*dx)
                rp2 = 0.25*dx*dx;

            /* The z-component of the area estimator */
            tAz = pos1[0]*pos2[1] - pos2[0]*pos1[1];

            /* The total area */
            Az += tAz;

            /* The local part of the area */
            if (k < numGrid)
                locAz(k) += tAz/rp2;
        }
    }

    estimator += locAz*Az;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RADIAL WINDING SUPERFLUID DENSITY ESTIMATOR CLASS -------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Initialize the variables needed to measure the radially averaged winding 
 *  number superfluid density.
 *  @see E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).
 *
******************************************************************************/
RadialWindingSuperfluidDensityEstimator::RadialWindingSuperfluidDensityEstimator
(const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    numGrid = NGRIDSEP;

    /* The spatial discretization */
    dR  = 0.5*path.boxPtr->side[0] / (1.0*numGrid);

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % ((n)*dR)));

    norm = 0.5 * constants()->T()/constants()->lambda();
    for (int n = 0; n < numGrid; n++) 
        norm(n) /= (M_PI*(2*n+1)*dR*dR*path.boxPtr->side[NDIM-1]);

    /* Initialize the local arrays */
    locWz.resize(numGrid);
    locWz = 0.0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RadialWindingSuperfluidDensityEstimator::~RadialWindingSuperfluidDensityEstimator() { 
    locWz.free();
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  
******************************************************************************/
void RadialWindingSuperfluidDensityEstimator::accumulate() {

    int numTimeSlices = path.numTimeSlices;

    beadLocator beadIndex;
    double Wz;
    dVec pos1,pos2;
    dVec vel;

    Wz = 0.0;
    locWz = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
            pos1 = path(beadIndex);
            pos2 = path(path.next(beadIndex));

            int k = int(sqrt(pos1[0]*pos1[0]+pos1[1]*pos1[1])/dR);

            /* The winding number estimator */
            vel = path.getVelocity(beadIndex)*path.boxPtr->periodic;
            Wz += vel[NDIM-1];

            /* The local part of the winding number */
            if (k < numGrid)
                locWz(k) += vel[NDIM-1];
        }
    }

    estimator += locWz*Wz;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RADIAL AREA SUPERFLUID DENSITY ESTIMATOR CLASS ----------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor
 * 
 *  Initialize the variables needed to measure the radially averaged winding 
 *  number superfluid density.
 *  @see E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).
 *
******************************************************************************/
RadialAreaSuperfluidDensityEstimator::RadialAreaSuperfluidDensityEstimator
(const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    numGrid = NGRIDSEP;

    /* The spatial discretization */
    dR  = 0.5*path.boxPtr->side[0] / (1.0*numGrid);

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % ((n)*dR)));

    norm = 0.5 * constants()->T()/constants()->lambda();
    for (int n = 0; n < numGrid; n++) 
        norm(n) /= (M_PI*(2*n+1)*dR*dR*path.boxPtr->side[NDIM-1]);

    /* Initialize the local arrays */
    locAz.resize(numGrid);
    locAz = 0.0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
RadialAreaSuperfluidDensityEstimator::~RadialAreaSuperfluidDensityEstimator() { 
    locAz.free();
}

/*************************************************************************//**
 *  Accumulate superfluid properties.
 *
 *  
******************************************************************************/
void RadialAreaSuperfluidDensityEstimator::accumulate() {

    int numTimeSlices = path.numTimeSlices;

    beadLocator beadIndex;
    double Az,rp2,tAz;
    dVec pos1,pos2;

    Az = 0.0;
    locAz = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
            pos1 = path(beadIndex);
            pos2 = path(path.next(beadIndex));
            rp2 = pos1[0]*pos1[0] + pos1[1]*pos1[1];

            int k = int(sqrt(rp2)/dR);

            /* The z-component of the area estimator */
            tAz = pos1[0]*pos2[1] - pos2[0]*pos1[1];

            /* The total area */
            Az += tAz;

            /* The local part of the winding number */
            if (k < numGrid)
                locAz(k) += tAz/rp2;
        }
    }

    estimator += locAz*Az;
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
(const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label):
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* This is a 'local' histogram estimator so we use the defined grid */
    numGrid = path.boxPtr->numGrid;
    initialize(3*numGrid);

    /* The smallest allowed radius */
    dR = 0.5*path.boxPtr->side[0]/numGrid;

    header = str(format("#%15d\n") % NGRIDSEP);
    header += str(format("#%15s%16s%16s") % "W:rho_s" % "A:rho_s" % "A^2");

    double C = 0.5*constants()->T()/constants()->lambda();

    /* The normalization constant for the local winding and area estimators. */
    for (int n = 0; n < numGrid; n++) {
        double dV = path.boxPtr->gridBoxVolume(n);
        norm(n) = C/dV;
        norm(n+numGrid) = C/dV;
        norm(n+2*numGrid) = C/dV;
    }

    /* Initialize the local arrays */
    locWz.resize(numGrid);
    locWz = 0.0;

    locAz.resize(numGrid);
    locAz = 0.0;

    locA2.resize(numGrid);
    locA2 = 0.0;

}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
LocalSuperfluidDensityEstimator::~LocalSuperfluidDensityEstimator() { 
    locWz.free();
    locAz.free();
    locA2.free();
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
    for (int n = 0; n < numGrid; n++) {
        for (int i = 0; i < int(numEst/numGrid); i++)
            (*outFilePtr) << format("%16.8E") % 
                (norm(n+i*numGrid)*estimator(n+i*numGrid)/totNumAccumulated);
        (*outFilePtr) << endl;
    }

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
    locA2 = 0.0;
    locWz = 0.0;

    beadLocator beadIndex;
    double Az,rp2,Wz;
    double tAz;
    dVec pos1,pos2;
    dVec vel;

    Az = Wz = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {

            beadIndex = slice,ptcl;
            pos1 = path(beadIndex);
            pos2 = path(path.next(beadIndex));
            int n = path.boxPtr->gridIndex(pos1);

            /*  The distance from the z-axis squared */
            //rp2 = pos1[0]*pos1[0] + pos1[1]*pos1[1];
            rp2 = pos1[0]*pos2[0] + pos1[1]*pos2[1];
            if (abs(rp2) < dR*dR)
                rp2 = dR*dR;

            /* The winding number estimator */
            vel = path.getVelocity(beadIndex)*path.boxPtr->periodic;
            Wz += vel[NDIM-1];

            /* The local part of the winding number */
            locWz(n) += vel[NDIM-1];

            /* The z-component of the area estimator */
            tAz = pos1[0]*pos2[1] - pos2[0]*pos1[1];

            /* The total area */
            Az += tAz;

            /* The two local components */
            locA2(n) += tAz; 
            locAz(n) += tAz/rp2;
        }
    }

    locWz *= Wz;
    locAz *= Az;
    locA2 *= Az;

    for (int n = 0; n < path.boxPtr->numGrid; n++) {
        estimator(n) += locWz(n);
        estimator(n+numGrid) += locAz(n);
        estimator(n+2*numGrid) += locA2(n);
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    initialize(1);

    header = str(format("%16s") % "diagonal");
    endLine = false;
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* We measure the average worm length, gap and cost.  It is an off-diagonal
     * estimator that is output to its own file */
    initialize(5);
    diagonal = false;

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
 * We setup the permutation cycle estimator to measure cycles which can 
 * contain up to 40 particles.
******************************************************************************/
PermutationCycleEstimator::PermutationCycleEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label)  : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* We just choose arbitrarily to only count cycles including up to 40 particles */
    maxNumCycles = 40;
    
    /* The permutatoin cycle estimator has its own file, and consists of 
     * maxNumCycles permutation cycles */
    initialize(maxNumCycles);

    /* Set estimator name and header, which contains the permutation cycle
     * numbers */
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

    double cycleNorm;

    if (numParticles > 0)
        cycleNorm = 1.0 / (1.0*numParticles);
    else
        cycleNorm = 0.0;

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
// LOCAL PERMUTATION CYCLE ESTIMATOR CLASS -----------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 * Constructor.
 * 
 * Local particle permutation number density histogram.
 *
 * We setup the permuatation cycle estimator to measure cycles which can 
 * contain up to 40 particles.
******************************************************************************/
LocalPermutationEstimator::LocalPermutationEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label)  : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* We just choose arbitrarily to only count cycles including up to 40 particles */
    maxNumCycles = 40;
    
    /* The local permutation cycle estimator has its own file, and consists of 
     * maxNumCycles permutation cycles */
    initialize(path.boxPtr->numGrid);

    /* vector to hold number of worldlines put into a grid space */
    numBeadInGrid.resize(estimator.size());

    /* Set estimator header */
    header = str(format("#%15d") % NGRIDSEP);

}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
LocalPermutationEstimator::~LocalPermutationEstimator() { 
    doBead.free();
}

/*************************************************************************//**
 *  Overload the output of the base class so that a running average
 *  is kept rather than keeping all data.
******************************************************************************/
void LocalPermutationEstimator::output() {

    /* Prepare the position file for writing over old data */
    communicate()->file(label)->reset();

    (*outFilePtr) << header;
    if (endLine)
        (*outFilePtr) << endl;

    /* Now write the running average of the estimator to disk */
    for (int n = 0; n < numEst; n++)
        (*outFilePtr) << format("%16.6E\n") %
            (estimator(n)/totNumAccumulated);

    communicate()->file(label)->rename();
}

/*************************************************************************//**
 * Accumulate permuation cycle.
 * 
 * Go through each worldline and count how many beads it includes. We must be
 * careful to avoid any double counting, which is acomplished by marking
 * beads we have touched.
 *
 * Once beads in a worldline are counted, loop through and label the 
 * positions by the corresponding permutation.
******************************************************************************/
void LocalPermutationEstimator::accumulate() {

    /* The start bead for each world line, and the moving index */
    beadLocator startBead;
    beadLocator beadIndex;

    int numWorldlines = path.numBeadsAtSlice(0);

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
            } while (!all(beadIndex==startBead)); // up to here, we have computed WL length only.

            /* Accumulate the cycle length counter */
            int cycleNum = int(wlLength / path.numTimeSlices);
            
            /* Loop through worldline again, this time binning the appropriate
             * permutation number (- 1) corresponding to its spatial coordinates.*/
            do {
                int nn = path.boxPtr->gridIndex(path(beadIndex));
                if ((cycleNum > 0) && (cycleNum <= maxNumCycles)){
                    estimator(nn) += (1.0*cycleNum - 1.0);
                    numBeadInGrid(nn) += 1;
                }

                beadIndex = path.next(beadIndex);
            } while (!all(beadIndex==startBead));
        } // doBead
    } // n

    /* Correct for multiple worldlines being in the same gridpoint
     * and compute normalization factor. */
    for (int i(0); i<int(estimator.size()); i++){
        if ((estimator(i) > 0) && (numBeadInGrid(i) > 0)){
            estimator(i) /= numBeadInGrid(i);
        }
    }

    /* resize array to store number of beads per grid */
    numBeadInGrid.resize(0);
    numBeadInGrid.resize(estimator.size());
        
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label),
    lpath(_path)
{

    sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());

    /* We chooose the maximum separation to be sqrt(NDIM)*min(L)/2 */
    dR = 0.5*sqrt(sum(path.boxPtr->periodic))*(blitz::min(path.boxPtr->side)) / (1.0*NOBDMSEP);

    /* This is an off-diagonal estimator*/
    initialize(NOBDMSEP);
    diagonal = false;

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
 *
 *  !!NB!! We only measure the OBDM when the tail is on an even time slice
******************************************************************************/
void OneBodyDensityMatrixEstimator::sample() {
    numSampled++;
    if ( frequency && 
         ((numSampled % frequency) == 0) &&
         (path.worm.isConfigDiagonal == diagonal) &&
         (path.worm.gap > 0) && (path.worm.gap <= constants()->Mbar())  &&
         ( (lpath.worm.tail[0] % 2) == 0) ) {

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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The spatial discretization */
    dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NPCFSEP);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < NPCFSEP; n++) 
        header.append(str(format("%16.3E") % ((n)*dR)));

    /* The normalization factor for the pair correlation function depends 
     * on the dimensionality */
//  TinyVector<double,3> gNorm;
//  gNorm[0] = 0.5;
//  gNorm[1] = 1.0/(4.0*M_PI);
//  gNorm[2] = 1.0/(8.0*M_PI);
//  norm(0) = 1.0;
//  for (int n = 1; n < NPCFSEP; n++)
//      norm(n) = (gNorm[NDIM-1]*path.boxPtr->volume) / (dR*pow(n*dR,NDIM-1));

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
        estimator += lnorm*(1.0*actionPtr->sepHist / 
            (1.0*sum(actionPtr->sepHist)));
    }
    else
        estimator += 0.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// STATIC STRUCTURE FACTOR ESTIMATOR CLASS -----------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Compute the static structure factor for a fixed set of q vectors.
 *  N.B. BROKEN: Doesn't work in general dimensions. 
******************************************************************************/
StaticStructureFactorEstimator::StaticStructureFactorEstimator(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    numq = 10;
    q.resize(numq);
    for (int n = 0; n < numq; n++)
        q(n) = (2.0*M_PI/path.boxPtr->side[0])*n;

    /* Initialize the accumulator intermediate scattering function*/
    sf.resize(numq);
    sf = 0.0;

    /* This is a diagonal estimator that gets its own file */
    initialize(numq);

    /* The magnitude of q */
    header = str(format("#%15.6E") % 0.0);
    for (int n = 1; n < numq; n++) 
        header.append(str(format("%16.6E") % q(n)));

    /* Utilize imaginary time translational symmetry */
    norm = 1.0/constants()->numTimeSlices();
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
StaticStructureFactorEstimator::~StaticStructureFactorEstimator() { 
}

/*************************************************************************//**
 *  measure the intermediate scattering function for each value of the 
 *  imaginary time separation tau.
 *
 *  We only compute this for N > 1 due to the normalization.
******************************************************************************/
void StaticStructureFactorEstimator::accumulate() {

    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();

    beadLocator bead1,bead2;  // The bead locator
    sf = 0.0; // initialize

    dVec pos1,pos2;      // The two bead positions
    dVec lq;             // The value of the wave-vector
    lq = 0.0;
    
    for (int nq = 0; nq < numq; nq++) {
        lq[0] = q(nq);

        /* Average over all initial time slices */
        for (int slice = 0; slice < numTimeSlices; slice++) {

            bead1[0] = bead2[0] = slice;

            for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(slice); bead1[1]++) {

                pos1 = path(bead1);
                double lq1 = dot(lq,pos1);

                for (bead2[1] = 0; bead2[1] < path.numBeadsAtSlice(slice); bead2[1]++) {

                    pos2 = path(bead2);
                    double lq2 = dot(lq,pos2);

                    sf(nq) += cos(lq1)*cos(lq2) + sin(lq1)*sin(lq2);

                } // bead2[1]
            } // bead1[1]
        } //slice
    } // nq

    estimator += sf/numParticles; 
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// INTERMEDIATE SCATTERING FUNCTION ESTIMATOR CLASS --------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Measure the intermediate scattering function at a hardcoded wavevector
 *  N.B. BROKEN: Doesn't work in general dimensions. 
******************************************************************************/
IntermediateScatteringFunctionEstimator::IntermediateScatteringFunctionEstimator(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    int numTimeSlices = constants()->numTimeSlices();

    /* these are the hard-coded q-vectors for now */
    numq = 3;
    Array <double, 1> qMag(numq);         // the q-vector magnitudes
    /* qMag.resize(numq); */
    qMag = 0.761,1.75,1.81;

    /* initialize the number of vectors with each magnitude */
    numqVecs.resize(numq);               
    numqVecs = 0;

    /* The allowable error in the q-vector magnitude */
    double eps = 2.0*M_PI/min(path.boxPtr->side)/sqrt(NDIM);
    eps *= eps;

    /* Determine the set of q-vectors that have these magintudes */
    for (int nq = 0; nq < numq; nq++)
    {
        double cq = qMag(nq);
        vector <dVec> qvecs;

        int maxComp = ceil(cq*blitz::min(path.boxPtr->side)/(2.0*M_PI))+1;
        int maxNumQ = ipow(2*maxComp + 1,NDIM);
        
        iVec qi;
        for (int n = 0; n < maxNumQ; n++) {
            for (int i = 0; i < NDIM; i++) {
                int scale = 1;
                for (int j = i+1; j < NDIM; j++) 
                    scale *= (2*maxComp + 1);
                qi[i] = (n/scale) % (2*maxComp + 1);

                /* Wrap into the appropriate winding sector */
                qi[i] -= (qi[i] > maxComp)*(2*maxComp + 1);
            }
            dVec qd = 2.0*M_PI*qi/path.boxPtr->side;

            /* Store the winding number */
            if (abs(dot(qd,qd)-cq*cq) < eps) {
                qvecs.push_back(qd);
                numqVecs(nq)++;
            }
        }
        q.push_back(qvecs);
    }

    /* get more accurate q-magnitudes */
    for (int nq = 0; nq < numq; nq++) 
        qMag(nq) = sqrt(dot(q[nq][0],q[nq][0]));

    /* Initialize the accumulator for the intermediate scattering function*/
    /* N.B. for now we hard-code three wave-vectors */
    isf.resize(numq*numTimeSlices);
    isf = 0.0;

    /* This is a diagonal estimator that gets its own file */
    initialize(numq*numTimeSlices);

    /* the q-values */
    header = str(format("#%15.6E") % qMag(0));
    for (int n = 1; n < numq; n++)
        header.append(str(format("%16.6E") % qMag(n)));
    header.append("\n");

    /* The imaginary time values */
    header.append(str(format("#%15.6E") % 0.0));
    for (int n = 1; n < numTimeSlices; n++) 
        header.append(str(format("%16.6E") % (constants()->tau()*n)));

    for (int nq = 1; nq < numq; nq++) {
        for (int n = 0; n < numTimeSlices; n++) 
            header.append(str(format("%16.6E") % (constants()->tau()*n)));
    }

    /* utilize imaginary time translational symmetry */
    norm = 1.0/numTimeSlices;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
IntermediateScatteringFunctionEstimator::~IntermediateScatteringFunctionEstimator() { 
}

/*************************************************************************//**
 *  measure the intermediate scattering function for each value of the 
 *  imaginary time separation tau.
 *
 *  We only compute this for N > 1 due to the normalization.
******************************************************************************/
void IntermediateScatteringFunctionEstimator::accumulate() {

    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();

    beadLocator bead1,bead2;  // The bead locator
    isf = 0.0; // initialize
    dVec pos1,pos2;         // The two bead positions

    /* q-magnitudes */
    for (int nq = 0; nq < numq; nq++) {

        /* q-vectors */
        for (const auto &cqvec : q[nq]) {

            /* Average over all initial time slices */
            for (bead1[0] = 0; bead1[0] < numTimeSlices; bead1[0]++) {

                /* compute for each tau separation */
                for (int tausep = 0;  tausep < numTimeSlices; tausep++){

                    bead2[0] = (bead1[0] + tausep) % numTimeSlices;

                    for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(bead1[0]); bead1[1]++) {

                        pos1 = path(bead1);
                        double lq1 = dot(cqvec,pos1);

                        for (bead2[1] = 0; bead2[1] < path.numBeadsAtSlice(bead2[0]); bead2[1]++) {

                            pos2 = path(bead2);
                            double lq2 = dot(cqvec,pos2);

                            isf(nq*numTimeSlices + tausep) += (cos(lq1)*cos(lq2) + sin(lq1)*sin(lq2))/numqVecs(nq);
                        } // bead2[1]
                    } // bead1[1]
                } //tausep   
            } //bead1[0]
        } //q-vectors
    } //q-magnitudes

    estimator += isf/numParticles;
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The spatial discretization */
    dR  = 0.5*path.boxPtr->side[0] / (1.0*NRADSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NRADSEP);

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
    double rsq;
    beadLocator beadIndex;
    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);
            rsq = 0.0;
            for (int i = 0; i < NDIM-1; i++)
                rsq += pos[i]*pos[i]; 
            int k = int(sqrt(rsq)/dR);
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
CylinderEnergyEstimator::CylinderEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    /* We compute three diagonal estimators, kinetic, potential and total energy
     * per particle */
    initialize(7);
    endLine = false;

    /* Set estimator header, we will always report the energy
     * first, hence the comment symbol*/
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

    int numParticles  = num1DParticles(path,maxR);
    int numTimeSlices = path.numTimeSlices;

    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * numTimeSlices);

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
    double t1 = 0.0;
    double t2 = 0.0;
    for (int slice = 0; slice < numTimeSlices; slice++) {
        eo = (slice % 2);
        t1 += actionPtr->derivPotentialActionLambda(slice,maxR);
        t2 += actionPtr->derivPotentialActionTau(slice,maxR);
        if (eo==0)
            totV  += actionPtr->potential(slice,maxR);
    }

    /* Normalize the action correction and the total potential*/
    t1 *= constants()->lambda()/(constants()->tau()*numTimeSlices);
    t2 /= 1.0*numTimeSlices;
    totV /= (0.5 * numTimeSlices);

    /* Perform all the normalizations and compute the individual energy terms */
    totK  += (classicalKinetic + t1);

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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* We compute three diagonal estimators, the total number of particles,
     * total number of particles squared and density. */
    initialize(3);

    /* Set estimator header */
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
    (const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    /* For now, we assume a maximum of 200 total particles. */
    maxNumParticles = 200;
    initialize(maxNumParticles);

    /* Set estimator header */
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
// CYLINDER LINEAR DENSITY ESTIMATOR CLASS -----------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Setup a vector estimator which measures the linear density along the 
 *  pore axis.
******************************************************************************/
CylinderLinearDensityEstimator::CylinderLinearDensityEstimator
    (const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The length of the cylinder */
    Lz = path.boxPtr->side[NDIM-1];

    /* The spatial discretization */
    dz = Lz / (1.0*NRADSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NRADSEP);

    /* The header is the first line which contains the spatial positions */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < NRADSEP; n++) 
        header.append(str(format("%16.3E") % (n*dz)));

    /* The normalization factor for the linear density*/
    norm = 1.0/(dz * constants()->numTimeSlices());
}

/*************************************************************************//**
 * Destructor.
******************************************************************************/
CylinderLinearDensityEstimator::~CylinderLinearDensityEstimator() { 
}

/*************************************************************************//**
 * Accumulate the linear density.
******************************************************************************/
void CylinderLinearDensityEstimator::accumulate() {

    dVec pos;
    beadLocator beadIndex;
    /* visit each bead */
    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            /* If we are inside the cutoff cylinder, accumulate the density 
             * histogram */
            if (include(pos,maxR)) {
                int k = int((0.5*Lz + pos[NDIM-1])/dz);
                if (k < NRADSEP)
                    estimator(k) += 1.0;
            }
        }
    }
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    windMax = 10;
    /* We compute a bunch of estimators here, the superfluid fraction, the winding
     * number in all possible dimensions, and the winding number histograms up to 
     * windMax windings. These are all diagonal estimators and we have our own
     * output file.*/
    initialize(4+2*windMax+1);

    /* Set estimator header */
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
CylinderOneBodyDensityMatrixEstimator::CylinderOneBodyDensityMatrixEstimator 
  (Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
   int _frequency, string _label) : 
      EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label),
      lpath(_path)
{
    sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());

    /* We chooose the maximum separation to be sqrt(NDIM)*L/2 */
    dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NOBDMSEP);

    /* This is an off-diagonal estimator that gets its own file */
    initialize(NOBDMSEP);
    diagonal = false;

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
         ( (lpath.worm.tail[0] % 2) == 0) ) {
        
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The spatial discretization */
    dR = 0.5*sqrt(sum(path.boxPtr->periodic))*path.boxPtr->side[NDIM-1] / (1.0*NPCFSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NPCFSEP);

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
        double lnorm = 1.0*sum(actionPtr->cylSepHist);
        lnorm /= 1.0*(N1D-1)/(1.0*N1D);
        estimator += 1.0*actionPtr->cylSepHist / (1.0*lnorm);
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

    if (baseSample() && (sum(actionPtr->cylSepHist) > 0)) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER LINEAR POTENTIAL ESTIMATOR CLASS ---------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  We compute the effective potential V(z) felt by a fictitious particle
 *  located along the axis of the cylinder (r=0).
******************************************************************************/
CylinderLinearPotentialEstimator::CylinderLinearPotentialEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The length of the cylinder */
    Lz = path.boxPtr->side[NDIM-1];

    /* The spatial discretization */
    dz = Lz / (1.0*NRADSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NRADSEP);

    /* The header is the first line which contains the spatial separations */
    header = str(format("#%15.3E") % 0.0);
    for (int n = 1; n < NRADSEP; n++) 
        header.append(str(format("%16.3E") % ((n)*dz)));
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderLinearPotentialEstimator::~CylinderLinearPotentialEstimator() { 
}

/*************************************************************************//**
 *  We determine what the effective potential along the axis of the pore.
******************************************************************************/
void CylinderLinearPotentialEstimator::accumulate() {

    double totV = 0.0;
    dVec r1,r2;         // The two bead positions
    r1 = 0.0;

    dVec sep;           // The bead separation
    beadLocator bead2;  // The bead locator

    /* sample all positions along the pore */
    for (int n = 0; n < NRADSEP; n++) {

        r1[NDIM-1] = -0.5*Lz + n*dz;

        /* Get the external potential */
        totV = 0.0;

        /* We sum up the interaction energy over all slices*/
        int numBeads = 0;
        for (int slice = 0; slice < path.numTimeSlices; slice++) {
            bead2[0] = slice;

            /* Sum over particles */
            for (bead2[1] = 0; bead2[1] < path.numBeadsAtSlice(slice); bead2[1]++) {

                r2 = path(bead2);
                if (!include(r2,maxR)) {
                    sep = r2 - r1;
                    path.boxPtr->putInBC(sep);
                    totV += actionPtr->interactionPtr->V(sep);
                    numBeads++;

                } // bead2 is outside maxR
            } // bead2
        } // slice

        totV /= 1.0*numBeads;

        /* Add the constant piece from the external potential */
        totV += actionPtr->externalPtr->V(r1);

        estimator(n) += totV;
    } // n
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
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The spatial discretization */
    dR  = 0.5*path.boxPtr->side[0] / (1.0*NRADSEP);

    /* This is a diagonal estimator that gets its own file */
    initialize(NRADSEP);
    radPot.resize(NRADSEP);

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
    dVec r1,r2;         // The two bead positions
    dVec sep;           // The bead separation

    beadLocator bead2;  // The bead locator

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
                    totV += actionPtr->interactionPtr->V(sep);
                } // bead2 is outside maxR
            } // bead2
        } // slice

        totV /= 1.0*path.numTimeSlices;
        totV += actionPtr->externalPtr->V(r1);

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
    dVec r1,r2;         // The two bead positions
    dVec sep;           // The bead separation
    double rad1,rad2;   // The two bead radii
    int nR;             // The bin number

    beadLocator bead1,bead2;    // The bead locators
    bool found1,found2;         // Are the beads in the central chain
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
                totV = actionPtr->externalPtr->V(r1);
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
                        totV += actionPtr->interactionPtr->V(sep);
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// BEGIN PIGS ESTIMATORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// POTENTIAL ENERGY ESTIMATOR CLASS ------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Setup the potential energy estimator.  We measure it on every other time slice.
******************************************************************************/
PotentialEnergyEstimator::PotentialEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* We measure on every other time slice */
    initialize( (constants()->numTimeSlices()-1)/2 +1);

    /* Set estimator header */
    header = str(format("#%15f") % 0.0 );
    for (int n = 2; n < constants()->numTimeSlices(); n+=2)
        header.append(str(format("%16f") % (n*constants()->tau()) ));
}

/*************************************************************************//**
 * Destructor.
******************************************************************************/
PotentialEnergyEstimator::~PotentialEnergyEstimator() { 
    // empty destructor
}

/*************************************************************************//**
 * Accumulate the potential energy
******************************************************************************/
void PotentialEnergyEstimator::accumulate() {

    /* The total tail correction */
    double tailV = (1.0*path.getTrueNumParticles()*path.getTrueNumParticles()
                    /path.boxPtr->volume)*actionPtr->interactionPtr->tailV;
    
    /* We use a simple operator estimator for V. */
    for (int slice = 0; slice <= path.numTimeSlices; slice+=2)
        estimator(slice/2) += actionPtr->potential(slice) + tailV;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// KINETIC ENERGY ESTIMATOR CLASS ----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 *
 *  Setup the kinetic energy estimator.  We measure it on every other time slice.
 ******************************************************************************/
KineticEnergyEstimator::KineticEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on every other time slice */
    initialize( (constants()->numTimeSlices()-1)/2 );
    
    /* Set estimator header */
    header = str(format("#%15f") % constants()->tau());
    for (int n = 2; n < (constants()->numTimeSlices()-1); n+=2)
        header.append(str(format("%16f") % ((n+1)*constants()->tau()) ));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
KineticEnergyEstimator::~KineticEnergyEstimator() {
    // empty destructor
}


/*************************************************************************//**
* Accumulate the potential energy
******************************************************************************/
void KineticEnergyEstimator::accumulate() {
    
    int numTimeSlices = constants()->numTimeSlices();
    int numParticles  = path.getTrueNumParticles();
    
    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau()*2.0);
    
    /* The classical contribution to the kinetic energy per particle
     * including the chemical potential */
    double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;
    
    beadLocator beadIndex;
    dVec vel,pos;
    for (int slice = 0; slice < (numTimeSlices-1); slice+=2) {
        double K = 0.0;
        for (int eo = 0; eo < 2; eo++){
            for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice+eo); ptcl++) {
                beadIndex = slice+eo,ptcl;
                vel = path.getVelocity(beadIndex);
                K -= dot(vel,vel);
            }
        }
        /* Normalize the accumulated link-action part */
        K *= kinNorm;
        
        /* Copmute the correction to the accumulated link-action part */
        double t1 = 0.0;
        for (int eo = 0; eo < 2; eo++){
            t1 += actionPtr->derivPotentialActionLambda(slice+eo);
        }

        t1 *= constants()->lambda()/(2.0*constants()->tau());
        
        /* Perform all the normalizations and compute the individual energy terms */
        K  += (classicalKinetic + t1);
        
        estimator(slice/2) += K;
    }
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PIGS ENERGY ESTIMATOR CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
* We measure the total Kinetic, Potential and E-mu*N as well as the kinetic,
* potential and total energy per particle.
******************************************************************************/
PigsEnergyEstimator::PigsEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* Set estimator header, we will always report the energy
     * first, hence the comment symbol*/
    header = str(format("#%15s%16s%16s%16s%16s%16s%16s")
                 % "K" % "V" % "E" % "E_mu" % "K/N" % "V/N" % "E/N");
    endLine = false;
    initialize(7);
}

/*************************************************************************//**
*  Destructor.
******************************************************************************/
PigsEnergyEstimator::~PigsEnergyEstimator() {
}

/*************************************************************************//**
*  Accumluate the energy.
*
*  We use the thermodynamic estimator for the kinetic energy and the
*  potential estimator for the potential energy. A possible shift is
*  made due to a tail correction.
******************************************************************************/
void PigsEnergyEstimator::accumulate() {
    
    double totK = 0.0;
    double totVop = 0.0;
    double totV = 0.0;
    
    int numParticles  = path.getTrueNumParticles();
    int numTimeSlices = path.numTimeSlices;
    
    /* The total tail correction */
    double tailV = (1.0*numParticles*numParticles/path.boxPtr->volume)
    * actionPtr->interactionPtr->tailV;
    
    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * 4.0);
    
    /* The classical contribution to the kinetic energy per particle
     * including the chemical potential */
    double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;
    
    /* We first calculate the kinetic energy.  Even though there
     * may be multiple mixing and swaps, it doesn't matter as we always
     * just advance one time step at a time, as taken care of through the
     * linking arrays.  This has been checked! */
    beadLocator beadIndex;
    dVec vel;
    int midSlice = (numTimeSlices-1)/2;
    for (int slice = midSlice-2; slice < midSlice+2; slice++) {
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
    double t1 = 0.0;
    double t2 = 0.0;
    for (int slice = midSlice-2; slice < midSlice+2; slice++) {
        t1 += actionPtr->derivPotentialActionLambda(slice);
        t2 += actionPtr->derivPotentialActionTau(slice);
    }
    totVop  += actionPtr->potential(midSlice);
    
    t1 *= constants()->lambda()/(constants()->tau()*4.0);
    t2 /= 4.0;
    
    /* Perform all the normalizations and compute the individual energy terms */
    totK += (classicalKinetic + t1);
    totV = t2 - t1 + tailV;
    
    totVop += tailV;
    
    totV = totVop;
    
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
// PIGS NON-LOCAL ENERGY ESTIMATOR CLASS -------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
* We measure the total, Kinetic and Potential energy.
******************************************************************************/
PigsThermoEnergyEstimator::PigsThermoEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* Set estimator name and header, we will always report the energy
     * first, hence the comment symbol*/
    header = str(format("#%15s%16s%16s") % "K" % "V" % "E");
    endLine = true;
    initialize(3);
}

/*************************************************************************//**
*  Destructor.
******************************************************************************/
PigsThermoEnergyEstimator::~PigsThermoEnergyEstimator() {
}

/*************************************************************************//**
*  Accumluate the energy.
*
*  We use the thermodynamic estimator for both the kinetic energy and the
*  potential energy. A possible shift is
*  made due to a tail correction.
******************************************************************************/

void PigsThermoEnergyEstimator::accumulate() {
    
    /* Set the period of the action */
    int actionPeriod = 1;
    if( (actionPtr->getActionName() == "gsf")||(actionPtr->getActionName() == "li_broughton") )
        actionPeriod = 2;
    
    double totK = 0.0;
    double totV = 0.0;
    
    int numParticles  = path.getTrueNumParticles();
    int numTimeSlices = path.numTimeSlices;
    
    /* The total tail correction */
    double tailV = (1.0*numParticles*numParticles/path.boxPtr->volume)
                        * actionPtr->interactionPtr->tailV;
    
    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau() * (2.0*actionPeriod));
    
    /* The classical contribution to the kinetic energy per particle
     * including the chemical potential */
    double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;
    
    /* We first calculate the kinetic energy.  Even though there
     * may be multiple mixing and swaps, it doesn't matter as we always
     * just advance one time step at a time, as taken care of through the
     * linking arrays.  This has been checked! */
    beadLocator beadIndex;
    dVec vel;
    int midSlice = (numTimeSlices-1)/2;
    for (int slice = midSlice-actionPeriod; slice < midSlice+actionPeriod; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            vel = path.getVelocity(beadIndex);
            totK -= dot(vel,vel);
        }
    }
    
    /* Normalize the accumulated link-action part */
    totK *= kinNorm;
    
    /* Now we compute the potential and kinetic energy.  We use an operator estimator
     * for V and the thermodynamic estimator for K */
    double t1 = 0.0;
    double t2 = 0.0;
    if(actionPtr->local){
        t1 += (0.5)*actionPtr->derivPotentialActionLambda(midSlice-actionPeriod);
        t2 += (0.5)*actionPtr->derivPotentialActionTau(midSlice-actionPeriod);
        for (int slice = midSlice-actionPeriod+1; slice < midSlice+actionPeriod; slice++) {
            t1 += actionPtr->derivPotentialActionLambda(slice);
            t2 += actionPtr->derivPotentialActionTau(slice);
        }
        t1 += (0.5)*actionPtr->derivPotentialActionLambda(midSlice+actionPeriod);
        t2 += (0.5)*actionPtr->derivPotentialActionTau(midSlice+actionPeriod);
    }else{
        for (int slice = midSlice-actionPeriod; slice < midSlice+actionPeriod; slice++) {
            t1 += actionPtr->derivPotentialActionLambda(slice);
            t2 += actionPtr->derivPotentialActionTau(slice);
        }
    }
    
    t1 *= constants()->lambda()/(constants()->tau()*(2.0*actionPeriod));
    t2 /= (2.0*actionPeriod);
    
    /* Perform all the normalizations and compute the individual energy terms */
    totK += (classicalKinetic + t1);
    totV = t2 - t1 + tailV;
    
    /* Now we accumulate the average total, kinetic and potential energy, 
     * as well as their values per particles. */
    estimator(0) += totK;
    estimator(1) += totV;
    estimator(2) += totK + totV;
    
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// TOTAL ENERGY ESTIMATOR CLASS ----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the total energy estimator.  We measure it on every time slice.
******************************************************************************/
TotalEnergyEstimator::TotalEnergyEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on every other time slice */
    initialize( (constants()->numTimeSlices()-1) );
    
    /* Set estimator header */
    header = str(format("#%15f") % constants()->tau());
    for (int n = 1; n < (constants()->numTimeSlices()-1); n++)
        header.append(str(format("%16f") % (n*constants()->tau()) ));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
TotalEnergyEstimator::~TotalEnergyEstimator() {
    // empty destructor
}


/*************************************************************************//**
* Accumulate the total energy
******************************************************************************/
void TotalEnergyEstimator::accumulate() {
    
    int numTimeSlices = constants()->numTimeSlices();
    int numParticles  = path.getTrueNumParticles();
    
    /* The kinetic normalization factor */
    double kinNorm = constants()->fourLambdaTauInv() / (constants()->tau());
    
    /* The classical contribution to the kinetic energy per particle
     * including the chemical potential */
    double classicalKinetic = (0.5 * NDIM / constants()->tau()) * numParticles;
    
    beadLocator beadIndex;
    dVec vel,pos;
    for (int slice = 0; slice < (numTimeSlices-1); slice++) {
        double K = 0.0;
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            vel = path.getVelocity(beadIndex);
            K -= dot(vel,vel);
        }
        /* Normalize the accumulated link-action part */
        K *= kinNorm;
        
        /* Perform all the normalizations and compute the individual energy terms */
        K  += (classicalKinetic);
        
        double dUdtau = actionPtr->derivPotentialActionTau(slice);
        estimator(slice) += K+dUdtau;
    }
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// THERMODYNAMIC POTENTIAL ENERGY ESTIMATOR CLASS ----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the total energy estimator.  We measure it on every time slice.
******************************************************************************/
ThermoPotentialEnergyEstimator::ThermoPotentialEnergyEstimator  (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label):
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on every other time slice */
    initialize( (constants()->numTimeSlices()-1) );
    
    /* Set estimator header */
    header = str(format("#%15f") % constants()->tau());
    for (int n = 1; n < (constants()->numTimeSlices()-1); n++)
        header.append(str(format("%16f") % (n*constants()->tau()) ));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
ThermoPotentialEnergyEstimator::~ThermoPotentialEnergyEstimator() {
    // empty destructor
}

/*************************************************************************//**
* Accumulate the total energy
******************************************************************************/
void ThermoPotentialEnergyEstimator::accumulate() {
    
    int numTimeSlices = constants()->numTimeSlices();
    
    for (int slice = 0; slice < (numTimeSlices-1); slice++) {        
        double dUdtau = actionPtr->derivPotentialActionTau(slice);
        double dUdlam = actionPtr->derivPotentialActionLambda(slice);
        estimator(slice) += dUdtau - (constants()->lambda()/constants()->tau())*dUdlam;
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// POSITION ESTIMATOR CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Setup the potential energy estimator.  We measure it on each time slice.
******************************************************************************/
PositionEstimator::PositionEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    /* We measure on each time slice */
    initialize(constants()->numTimeSlices());

    /* Set estimator header */
    header = str(format("#%15d") % 0);
    for (int n = 1; n < constants()->numTimeSlices(); n++) 
        header.append(str(format("%16d") % n));
}

/*************************************************************************//**
 * Destructor.
******************************************************************************/
PositionEstimator::~PositionEstimator() { 
    // empty destructor
}

/*************************************************************************//**
 * Accumulate the potential energy
******************************************************************************/
void PositionEstimator::accumulate() {

    /* We use a simple operator estimator for V. */
    beadLocator beadIndex;
    double x;
    for (beadIndex[0] = 0; beadIndex[0] < path.numTimeSlices; ++beadIndex[0]) {
        x = 0.0;
        for (beadIndex[1] = 0; beadIndex[1] <
                path.numBeadsAtSlice(beadIndex[0]); ++beadIndex[1]) 
            x += path(beadIndex)[0];
        estimator(beadIndex[0]) += x;
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PARTICLE RESOLVED POSITION ESTIMATOR CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the potential energy estimator.  We measure it on each time slice.
******************************************************************************/
ParticleResolvedPositionEstimator::ParticleResolvedPositionEstimator(const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label)  :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on each time slice */
    initialize(constants()->initialNumParticles()*2);
    
    /* Set estimator header */
    header = str(format("#%15d") % 0);
    header.append(str(format("%16d") % 0));
    for (int n = 1; n < constants()->initialNumParticles(); n++){
        header.append(str(format("%16d") % n));
        header.append(str(format("%16d") % n));
    }
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
ParticleResolvedPositionEstimator::~ParticleResolvedPositionEstimator() {
    // empty destructor
}

/*************************************************************************//**
* Accumulate the particle resolved positions
******************************************************************************/
void ParticleResolvedPositionEstimator::accumulate() {
    
    /* We use a simple operator estimator for V. */
    beadLocator beadIndex;
    double x;
    
    beadIndex[0] = (path.numTimeSlices-1)/2;
    for (beadIndex[1] = 0; beadIndex[1] <path.numBeadsAtSlice(beadIndex[0]);
         ++beadIndex[1]){
        x = path(beadIndex)[0];
        if ( beadIndex[1] < constants()->initialNumParticles()){
            estimator(beadIndex[1]*2) += x;
            estimator(beadIndex[1]*2+1) += x*x;
        }
    }
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PARTICLE CORRELATION ESTIMATOR CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the potential energy estimator.  We measure it on each time slice.
******************************************************************************/
ParticleCorrelationEstimator::ParticleCorrelationEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on each time slice */
    initialize(constants()->initialNumParticles()-1);
    
    /* Set estimator header */
    header = str(format("#%15d") % 1);
    for (int n = 2; n < constants()->initialNumParticles(); n++)
        header.append(str(format("%16d") % n));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
ParticleCorrelationEstimator::~ParticleCorrelationEstimator() {
    // empty destructor
}

/*************************************************************************//**
* Accumulate the particle resolved positions
******************************************************************************/
void ParticleCorrelationEstimator::accumulate() {
    
    /* We use a simple operator estimator for V. */
    beadLocator beadIndex0,beadIndex;
    
    beadIndex0[0] = (path.numTimeSlices-1)/2;
    beadIndex[0] = beadIndex0[0];
    dVec r;
    r = 0.0;
    
    beadIndex0[1]=0;
    for (beadIndex[1] = 1; beadIndex[1] <path.numBeadsAtSlice(beadIndex[0]);
         ++beadIndex[1]){
        if ( beadIndex[1] < constants()->initialNumParticles())
            r = path.getSeparation(beadIndex0,beadIndex);
            estimator(beadIndex[1]-1) += dot(r,r);
    }
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Velocity ESTIMATOR CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the velocity estimator.  We measure it on each time slice.
******************************************************************************/
VelocityEstimator::VelocityEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on each time slice */
    initialize(constants()->numTimeSlices()-1);
    
    /* Set estimator header */
    header = str(format("#%15d") % 0);
    for (int n = 1; n < constants()->numTimeSlices()-1; n++)
        header.append(str(format("%16d") % n));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
VelocityEstimator::~VelocityEstimator() {
    // empty destructor
}

/*************************************************************************//**
* Accumulate the velocity                                                                            
* ******************************************************************************/
void VelocityEstimator::accumulate() {
    
    beadLocator beadIndex;
    dVec vel;
    
    beadIndex[1] = 0;
    for (beadIndex[0] = 0; beadIndex[0] < (path.numTimeSlices-1); ++beadIndex[0]) {
        if ( (path.breakSlice > 0) && (beadIndex[0] == path.breakSlice)
                                    &&( all(path.next(beadIndex)==XXX) ) ){
            beadLocator nextBead = beadIndex;
            nextBead[0]++;
            vel = path.getSeparation(beadIndex,nextBead);
        }else
            vel = path.getVelocity(beadIndex);
        estimator(beadIndex[0]) += sqrt(dot(vel,vel));
    }
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SubregionOccupation ESTIMATOR CLASS ---------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
*
*  Setup the velocity estimator.  We measure it on each time slice.
******************************************************************************/
SubregionOccupationEstimator::SubregionOccupationEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    
    /* We measure on each time slice */
    initialize(3);
    
    /* Set estimator header */
    header = str(format("#%15s") % "Z");
    header.append(str(format("%16s") % "pA"));
    header.append(str(format("%16s") % "pB"));
}

/*************************************************************************//**
* Destructor.
******************************************************************************/
SubregionOccupationEstimator::~SubregionOccupationEstimator() {
    // empty destructor
}

/*************************************************************************//**
* Accumulate the velocity
* ******************************************************************************/
void SubregionOccupationEstimator::accumulate() {
    
    beadLocator midLIndex, midRIndex;
    double rho0Mid,Z,pA,pB;
    
    midLIndex[0] = path.breakSlice;
    midRIndex[0] = path.breakSlice+1;
    midLIndex[1] = 0;
    midRIndex[1] = 0;
    
    Z = 0.0;
    if( path.inSubregionA(midRIndex) ){
        rho0Mid = actionPtr->rho0(midLIndex,midRIndex,1);
        Z = rho0Mid;
        pA = rho0Mid;
        pB = 0.0;
    }else if( path.inSubregionB(midRIndex) ){
        Z = 1.0;
        pA = 0.0;
        pB = 1.0;
    }else{
        Z = 0.0;
        pA = 0.0;
        pB = 0.0;
        cout << "Error!! Invalid configuration for SubregionOccupatoin!!" << endl;
    }
    
    estimator(0) += Z;
    estimator(1) += pA;
    estimator(2) += pB;
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
PIGSOneBodyDensityMatrixEstimator::PIGSOneBodyDensityMatrixEstimator (Path &_path,
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label),
    lpath(_path)
{

    sqrt2LambdaTau = sqrt(2.0 * constants()->lambda() * constants()->tau());

    /* We chooose the maximum separation to be sqrt(NDIM)*min(L)/2 */
    dR = 0.5*sqrt(sum(path.boxPtr->periodic))*(blitz::min(path.boxPtr->side)) / (1.0*NOBDMSEP);

    /* This is an off-diagonal estimator*/
    initialize(NOBDMSEP);
    diagonal = false;

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
PIGSOneBodyDensityMatrixEstimator::~PIGSOneBodyDensityMatrixEstimator() {
}

/*************************************************************************//**
 *  Sample the OBDM.
 *
 *  We overload the sample method for the one body density matrix as we
 *  only want to measure when the gap is not too large, otherwise we will be
 *  dominated by tiny close probabilities.
******************************************************************************/
void PIGSOneBodyDensityMatrixEstimator::sample() {
    numSampled++;
    if ( frequency &&
         ((numSampled % frequency) == 0)
         //  && (path.worm.isConfigDiagonal == diagonal) &&
         //(path.worm.gap > 0) && (path.worm.gap <= constants()->Mbar())  &&
         //(actionPtr->eFactor[(lpath.worm.tail[0] % 2)] < EPS)
       ){

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
inline dVec PIGSOneBodyDensityMatrixEstimator::getRandomVector(const double r) {
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
 * Accumulate the OBDM.
 *
 * We perform a fake close move, where the head of the worm is advanced to
 * a position a distance 'r' away from the tail but at the same time slice.
 * The probability of excepting such a move is equal (up to normalization)
 * to the one body density matrix.
*****************************************************************************/
void PIGSOneBodyDensityMatrixEstimator::accumulate() {
   
   beadLocator beadIndexL,beadIndexR;
   beadIndexL[0] = lpath.breakSlice;
   beadIndexL[1] = path.brokenWorldlinesL.front();
   beadIndexR[0] = lpath.breakSlice+1;
   beadIndexR[1] = path.brokenWorldlinesR.front();

   oldTailPos = lpath(beadIndexR);
   //oldAction = actionPtr->potentialAction(beadIndexR);

   dVec pos;
   pos = 0.0;

   /* Connection the broken beads*/
   //lpath.next(beadIndexL) = beadIndexR;
   //lpath.prev(beadIndexR) = beadIndexL;

    //for (int p = 0; p < numReps; p++) {

        /* Now we loop through all possible separations, evaluating the potential
         * action */
        for (int n = 0; n < NOBDMSEP; n++) {

            newAction = 0.0;
            ++numAttempted;

            /* Assign the new displaced tail position */
            newTailPos = oldTailPos + getRandomVector(n*dR);
            lpath.boxPtr->putInside(newTailPos);
            lpath.updateBead(beadIndexR,newTailPos);

            /* Compute the free particle density matrix */
            
            /* CMH: Need to sample winding sector for general dimension here */
            //double rho0 = actionPtr->rho01D(beadIndexL,beadIndexR,1);
            
            double rho0 = 1.0;

            /* action shift coming from a finite chemical potential */
            //double muShift = lpath.worm.gap*constants()->mu()*constants()->tau();

            /* Copmute the potential the potential action */
            newAction += actionPtr->potentialAction(beadIndexR);

            //double expAction = exp(-newAction + oldAction + muShift);
            double expAction = exp(-0.5*newAction + 0.5*oldAction);
            //double expAction = exp(-newAction);

            estimator(n) += rho0*expAction;

            /* Record the probability of accepting the move */
            if (random.randExc() < rho0*expAction)
                ++numAccepted;

        } // end for n

    //} // end for k

    /* Now we must undo any damge we have caused by reverting the tail to its previous position*/
    lpath.updateBead(beadIndexR,oldTailPos);
    //lpath.next(beadIndexL) = XXX;
    //lpath.prev(beadIndexR) = XXX;
}

/*************************************************************************//**
 *  For the one body density matrix estimator, we would like to output
 *  the acceptance information for the accumulate trial move.
******************************************************************************/
void PIGSOneBodyDensityMatrixEstimator::outputFooter() {

    (*outFilePtr) << format("# accepted: %16.8E attempted: %16.8E ratio: %16.4E\n")
        % (1.0*numAccepted) % (1.0*numAttempted) % (1.0*numAccepted/(1.0*numAttempted));
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// DOUBLED ESTIMATOR BASE CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
******************************************************************************/
DoubledEstimator::DoubledEstimator (const Path &_path, const Path &_path2,
        ActionBase *_actionPtr, ActionBase* _actionPtr2, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label), path2(_path2) 
{
}

/*************************************************************************//**
*  Destructor.
******************************************************************************/
DoubledEstimator::~DoubledEstimator() {
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Swap Estimator Class ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
******************************************************************************/
SwapEstimator::SwapEstimator (Path &_path, Path &_path2, ActionBase *_actionPtr, 
        ActionBase *_actionPtr2, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
DoubledEstimator(_path,_path2,_actionPtr,_actionPtr2,_random,_maxR,_frequency,_label), 
    lpath(_path),lpath2(_path2),
    actionPtr(_actionPtr), 
    actionPtr2(_actionPtr2) 
{
    
    /* Set estimator header */
    header = str(format("#%15s%16s%16s") % "Z" % "S" % "ZS");
    initialize(3);
}

/*************************************************************************//**
*  Destructor.
******************************************************************************/
SwapEstimator::~SwapEstimator() {
}

/*************************************************************************//**
* Accumulate the swap estimator for open paths.
******************************************************************************/
void SwapEstimator::accumulate() {
    if (lpath.breakSlice > 0)
        accumulateOpen();
    else
        accumulateClosed();
}

/*************************************************************************//**
* Accumulate the swap estimator for open paths.
******************************************************************************/
void SwapEstimator::accumulateOpen() {
    
    double S,Z; // Swap and normalization
    
    beadLocator beadIndexL,beadIndexR,beadIndexL2,beadIndexR2;
    beadIndexL[0] = lpath.breakSlice;
    beadIndexR[0] = lpath.breakSlice+1;
    beadIndexL2[0] = lpath2.breakSlice;
    beadIndexR2[0] = lpath2.breakSlice+1;
    
    /* Old tail positions & potential actions */
    vector<dVec> oldTailPos(lpath.brokenWorldlinesR.size());
    vector<dVec> oldTailPos2(lpath2.brokenWorldlinesR.size());

    for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
        beadIndexR[1] = lpath.brokenWorldlinesR[i];
        oldTailPos[i] = lpath(beadIndexR);
    }
    for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
        beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
        oldTailPos2[i] = lpath2(beadIndexR2);
    }
    
    /* Compute direct propagator for Z */
    double oldPotAction,oldPotAction2;
    double rho0Dir, rho0Dir2;
    
    /* Compute direct potential actions for both paths */
    /* loop over broken beads */
    oldPotAction = 0.0;
    for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
        beadIndexR[1] = lpath.brokenWorldlinesR[i];
        oldPotAction += actionPtr->potentialAction(beadIndexR);
    }
    oldPotAction2 = 0.0;
    for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
        beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
        oldPotAction2 += actionPtr2->potentialAction(beadIndexR2);
    }
    
    
    /* Initialize permuation vector */
    vector<int> permutation(lpath.brokenWorldlinesR.size());
    int Nperm = 0;
    double rho0perm;
    for( uint32 i=0; i<permutation.size(); i++)
        permutation[i] = i;
    
    /* Compute direct free particle density matricies for both paths */
    if ( lpath.brokenWorldlinesL.size() == 0){
        rho0Dir = 1.0;
    }else{
        rho0Dir = 0.0;
        /* sum over permutations */
        do {
            /* sum over broken worldlines */
            rho0perm = 1.0;
            for( uint32 i=0; i<permutation.size(); i++){
                beadIndexL[1] = lpath.brokenWorldlinesL[i];
                beadIndexR[1] = lpath.brokenWorldlinesR[permutation[i]];
                rho0perm *= actionPtr->rho0(beadIndexL,beadIndexR,1);
            }
            rho0Dir += rho0perm;
            Nperm++;
        } while (next_permutation(permutation.begin(),permutation.end()));
        rho0Dir *= 1.0/((double)Nperm);
    }
    
    /* Re-initialize permuation vector */
    permutation.resize(lpath2.brokenWorldlinesR.size());
    for( uint32 i=0; i<permutation.size(); i++)
        permutation[i] = i;
    Nperm = 0;
    if ( lpath2.brokenWorldlinesL.size() == 0){
        rho0Dir2 = 1.0;
    }else{
        rho0Dir2 = 0.0;
        /* sum over permutations */
        do {
            /* sum over broken worldlines */
            rho0perm = 1.0;
            for( uint32 i=0; i<permutation.size(); i++){
                beadIndexL2[1] = lpath2.brokenWorldlinesL[i];
                beadIndexR2[1] = lpath2.brokenWorldlinesR[permutation[i]];
                rho0perm *= actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
            }
            rho0Dir2 += rho0perm;
            Nperm++;
        } while (next_permutation(permutation.begin(),permutation.end()));
        rho0Dir2 *= 1.0/((double)Nperm);
    }
    
    /* Compute the swapped propagator for S */
   
    if( (lpath.brokenWorldlinesR.size() == lpath2.brokenWorldlinesR.size()) ){
        
        if ( lpath.brokenWorldlinesR.size() == 0){
            S = 1.0;
            Z = 1.0;
        }else{
            /* Swap the tail positions */
            for( uint32 i=0; i<oldTailPos2.size(); i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                lpath.updateBead(beadIndexR,oldTailPos2[i]);
            }
            for( uint32 i=0; i<oldTailPos.size(); i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                lpath2.updateBead(beadIndexR2,oldTailPos[i]);
            }
            
            /* Compute the swapped free particle density matrix */
            /* Re-initialize permuation vector */
            for( uint32 i=0; i<permutation.size(); i++)
                permutation[i] = i;
            Nperm = 0;
            double rho0Swap = 0.0;
            /* sum over permutations */
            do {
                /* sum over broken worldlines */
                rho0perm = 1.0;
                for( uint32 i=0; i<permutation.size(); i++){
                    beadIndexL[1] = lpath.brokenWorldlinesL[i];
                    beadIndexR[1] = lpath.brokenWorldlinesR[permutation[i]];
                    rho0perm *= actionPtr->rho0(beadIndexL,beadIndexR,1);
                }
                rho0Swap += rho0perm;
                Nperm++;
            } while (next_permutation(permutation.begin(),permutation.end()));
            rho0Swap*= 1.0/((double)Nperm);
            
            /* Re-initialize permuation vector */
            for( uint32 i=0; i<permutation.size(); i++)
                permutation[i] = i;
            Nperm = 0;
            double rho0Swap2 = 0.0;
            /* sum over permutations */
            do {
                /* sum over broken worldlines */
                rho0perm = 1.0;
                for( uint32 i=0; i<permutation.size(); i++){
                    beadIndexL2[1] = lpath2.brokenWorldlinesL[i];
                    beadIndexR2[1] = lpath2.brokenWorldlinesR[permutation[i]];
                    rho0perm *= actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
                }
                rho0Swap2 += rho0perm;
                Nperm++;
            } while (next_permutation(permutation.begin(),permutation.end()));
            rho0Swap2*= 1.0/((double)Nperm);

            /* Compute the potential the potential action */
            /* Compute direct potential actions for both paths */
            /* loop over broken beads */
            double newPotAction = 0.0;
            for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                newPotAction += actionPtr->potentialAction(beadIndexR);
            }
            double newPotAction2 = 0.0;
            for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                newPotAction2 += actionPtr2->potentialAction(beadIndexR2);
            }

            /* Now we must undo any damge we have caused by reverting the tail to its previous position*/
            for( uint32 i=0; i<oldTailPos.size(); i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                lpath.updateBead(beadIndexR,oldTailPos[i]);
            }
            for( uint32 i=0; i<oldTailPos2.size(); i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                lpath2.updateBead(beadIndexR2,oldTailPos2[i]);
            }
            S = rho0Swap*rho0Swap2*exp(-(0.5)*(newPotAction+newPotAction2)+(0.5)*(oldPotAction+oldPotAction2));
            Z = rho0Dir*rho0Dir2;
        }
    }else{
        S = 0.0;
        Z = rho0Dir*rho0Dir2;
    }
    
    /* In the normalization factor we must account for the factor of 1/2 in the 
       open path weight */
    //Z = rho0Dir*rho0Dir2*exp(-(oldPotAction+oldPotAction2));
    //S = rho0Swap*rho0Swap2*exp(-(newPotAction+newPotAction2));
    
    estimator(0) += Z;
    estimator(1) += S;
    estimator(2) += S*Z;
}

/*************************************************************************//**
* Accumulate the swap estimator for closed.
******************************************************************************/
void SwapEstimator::accumulateClosed() {
    
    double S,Z; // Swap and normalization
    
    /* We assume the broken bead is bead 0 */
    int swapSlice = (lpath.numTimeSlices-1)/2-1;
    
    beadLocator beadIndexL,beadIndexR,beadIndexL2,beadIndexR2;
    beadIndexL[0] = swapSlice;
    beadIndexR[0] = swapSlice+1;
    beadIndexL2[0] = swapSlice;
    beadIndexR2[0] = swapSlice+1;
    
    S = 0.0;
    Z = 1.0;
    for (int brokenBeadIndex1=0; brokenBeadIndex1<lpath.getNumParticles(); brokenBeadIndex1++){
        
        beadIndexL[1] = brokenBeadIndex1;
        beadIndexR[1] = brokenBeadIndex1;
        
        /* Old tail positions & potential actions */
        dVec oldTailPos = lpath(beadIndexR);
        double oldPotAction = actionPtr->potentialAction(beadIndexR);
        
        /* Compute direct free particle density matricies */
        double rho0Dir = actionPtr->rho0(beadIndexL,beadIndexR,1);
        
        for (int brokenBeadIndex2=0; brokenBeadIndex2<lpath2.getNumParticles(); brokenBeadIndex2++){
            
            beadIndexL2[1] = brokenBeadIndex2;
            beadIndexR2[1] = brokenBeadIndex2;
            
            /* Old tail positions & potential actions */
            dVec oldTailPos2 = lpath2(beadIndexR2);
            double oldPotAction2 = actionPtr2->potentialAction(beadIndexR2);
            
            /* Compute direct free particle density matricies */
            double rho0Dir2 = actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
            
            /* Swap the tail positions */
            lpath.updateBead(beadIndexR,oldTailPos2);
            lpath2.updateBead(beadIndexR2,oldTailPos);
            
            /* Compute the swapped free particle density matrix */
            double rho0Swap = actionPtr->rho0(beadIndexL,beadIndexR,1);
            double rho0Swap2 = actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
            
            /* Copmute the potential the potential action */
            double newPotAction = actionPtr->potentialAction(beadIndexR);
            double newPotAction2 = actionPtr2->potentialAction(beadIndexR2);
            
            /* Now we must undo any damge we have caused by reverting the tail to its previous position*/
            lpath.updateBead(beadIndexR,oldTailPos);
            lpath2.updateBead(beadIndexR2,oldTailPos2);
           
            S += (rho0Swap*rho0Swap2/(rho0Dir*rho0Dir2))*exp(-(0.5)*(newPotAction+newPotAction2)+(0.5)*(oldPotAction+oldPotAction2));
        }
    }
    
    S /= lpath.getNumParticles()*lpath2.getNumParticles();
    
    estimator(0) += Z;
    estimator(1) += S;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ----------- Entanglement of Particles Estimator Class ---------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
*  Constructor.
******************************************************************************/
EntPartEstimator::EntPartEstimator (Path &_path, Path &_path2,
                              ActionBase *_actionPtr,ActionBase *_actionPtr2,
                              const MTRand &_random, double _maxR,
                              int _frequency, string _label) :
DoubledEstimator(_path,_path2,_actionPtr, _actionPtr2,_random,_maxR,_frequency,_label), 
    lpath(_path),lpath2(_path2),
    actionPtr(_actionPtr), 
    actionPtr2(_actionPtr2) 
{
    stringstream headerSS;
    
    /* Set estimator name and header */
    header = str(format("#%15s") % "Z" );
    for (int n = 1; n < constants()->initialNumParticles(); n++){
        headerSS.str("");
        headerSS << 'Z' << n;
        header.append(str(format("%16s") % headerSS.str().c_str()));
        headerSS.str("");
        headerSS << 'S' << n;
        header.append(str(format("%16s") % headerSS.str().c_str()));
    }
    initialize(1+2*(constants()->initialNumParticles()-1));
}

/*************************************************************************//**
*  Destructor.
******************************************************************************/
EntPartEstimator::~EntPartEstimator() {
}

/*************************************************************************//**
* Accumulate the swap estimator for open paths.
******************************************************************************/
void EntPartEstimator::accumulate() {
    double S,Z; // Swap and normalization
    
    bool nMatch = false; //True if both paths have the same number of broken paths>0
    int nBin;   // number of broken paths if nMatch=True;
    
    beadLocator beadIndexL,beadIndexR,beadIndexL2,beadIndexR2;
    beadIndexL[0] = lpath.breakSlice;
    beadIndexR[0] = lpath.breakSlice+1;
    beadIndexL2[0] = lpath2.breakSlice;
    beadIndexR2[0] = lpath2.breakSlice+1;
    
    /* Old tail positions & potential actions */
    vector<dVec> oldTailPos(lpath.brokenWorldlinesR.size());
    vector<dVec> oldTailPos2(lpath2.brokenWorldlinesR.size());
    
    for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
        beadIndexR[1] = lpath.brokenWorldlinesR[i];
        oldTailPos[i] = lpath(beadIndexR);
    }
    for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
        beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
        oldTailPos2[i] = lpath2(beadIndexR2);
    }
    
    /* Compute direct propagator for Z */
    double oldPotAction,oldPotAction2;
    double rho0Dir, rho0Dir2;
    
    /* Compute direct potential actions for both paths */
    /* loop over broken beads */
    oldPotAction = 0.0;
    for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
        beadIndexR[1] = lpath.brokenWorldlinesR[i];
        oldPotAction += actionPtr->potentialAction(beadIndexR);
    }
    oldPotAction2 = 0.0;
    for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
        beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
        oldPotAction2 += actionPtr2->potentialAction(beadIndexR2);
    }
    
    
    /* Initialize permuation vector */
    vector<int> permutation(lpath.brokenWorldlinesR.size());
    int Nperm = 0;
    double rho0perm;
    for( uint32 i=0; i<permutation.size(); i++)
        permutation[i] = i;
    
    /* Compute direct free particle density matricies for both paths */
    if ( lpath.brokenWorldlinesL.size() == 0){
        rho0Dir = 1.0;
    }else{
        rho0Dir = 0.0;
        /* sum over permutations */
        do {
            /* sum over broken worldlines */
            rho0perm = 1.0;
            for( uint32 i=0; i<permutation.size(); i++){
                beadIndexL[1] = lpath.brokenWorldlinesL[i];
                beadIndexR[1] = lpath.brokenWorldlinesR[permutation[i]];
                rho0perm *= actionPtr->rho0(beadIndexL,beadIndexR,1);
            }
            rho0Dir += rho0perm;
            Nperm++;
        } while (next_permutation(permutation.begin(),permutation.end()));
        rho0Dir *= 1.0/((double)Nperm);
    }
    
    /* Re-initialize permuation vector */
    permutation.resize(lpath2.brokenWorldlinesR.size());
    for( uint32 i=0; i<permutation.size(); i++)
        permutation[i] = i;
    Nperm = 0;
    if ( lpath2.brokenWorldlinesL.size() == 0){
        rho0Dir2 = 1.0;
    }else{
        rho0Dir2 = 0.0;
        /* sum over permutations */
        do {
            /* sum over broken worldlines */
            rho0perm = 1.0;
            for( uint32 i=0; i<permutation.size(); i++){
                beadIndexL2[1] = lpath2.brokenWorldlinesL[i];
                beadIndexR2[1] = lpath2.brokenWorldlinesR[permutation[i]];
                //cout << beadIndexL2[1] << '\t' << beadIndexR2[1] << '\t' << permutation[i] << endl;
                rho0perm *= actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
            }
            rho0Dir2 += rho0perm;
            Nperm++;
        } while (next_permutation(permutation.begin(),permutation.end()));
        rho0Dir2 *= 1.0/((double)Nperm);
    }
    
    /* Compute the swapped propagator for S */
    
    if( (lpath.brokenWorldlinesR.size() == lpath2.brokenWorldlinesR.size()) ){
        
        if ( lpath.brokenWorldlinesR.size() == 0){
            S = 1.0;
            Z = 1.0;
        }else{
            /* Swap the tail positions */
            for( uint32 i=0; i<oldTailPos2.size(); i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                lpath.updateBead(beadIndexR,oldTailPos2[i]);
            }
            for( uint32 i=0; i<oldTailPos.size(); i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                lpath2.updateBead(beadIndexR2,oldTailPos[i]);
            }
            
            /* Compute the swapped free particle density matrix */
            /* Re-initialize permuation vector */
            for( uint32 i=0; i<permutation.size(); i++)
                permutation[i] = i;
            Nperm = 0;
            double rho0Swap = 0.0;
            /* sum over permutations */
            do {
                /* sum over broken worldlines */
                rho0perm = 1.0;
                for( uint32 i=0; i<permutation.size(); i++){
                    beadIndexL[1] = lpath.brokenWorldlinesL[i];
                    beadIndexR[1] = lpath.brokenWorldlinesR[permutation[i]];
                    rho0perm *= actionPtr->rho0(beadIndexL,beadIndexR,1);
                }
                rho0Swap += rho0perm;
                Nperm++;
            } while (next_permutation(permutation.begin(),permutation.end()));
            rho0Swap*= 1.0/((double)Nperm);
            
            /* Re-initialize permuation vector */
            for( uint32 i=0; i<permutation.size(); i++)
                permutation[i] = i;
            Nperm = 0;
            double rho0Swap2 = 0.0;
            /* sum over permutations */
            do {
                /* sum over broken worldlines */
                rho0perm = 1.0;
                for( uint32 i=0; i<permutation.size(); i++){
                    beadIndexL2[1] = lpath2.brokenWorldlinesL[i];
                    beadIndexR2[1] = lpath2.brokenWorldlinesR[permutation[i]];
                    rho0perm *= actionPtr2->rho0(beadIndexL2,beadIndexR2,1);
                }
                rho0Swap2 += rho0perm;
                Nperm++;
            } while (next_permutation(permutation.begin(),permutation.end()));
            rho0Swap2*= 1.0/((double)Nperm);
            
            /* Compute the potential the potential action */
            /* Compute direct potential actions for both paths */
            /* loop over broken beads */
            double newPotAction = 0.0;
            for( uint32 i=0; i<lpath.brokenWorldlinesR.size();i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                newPotAction += actionPtr->potentialAction(beadIndexR);
            }
            double newPotAction2 = 0.0;
            for( uint32 i=0; i<lpath2.brokenWorldlinesR.size();i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                newPotAction2 += actionPtr2->potentialAction(beadIndexR2);
            }
            
            /* Now we must undo any damge we have caused by reverting the tail to its previous position*/
            for( uint32 i=0; i<oldTailPos.size(); i++){
                beadIndexR[1] = lpath.brokenWorldlinesR[i];
                lpath.updateBead(beadIndexR,oldTailPos[i]);
            }
            for( uint32 i=0; i<oldTailPos2.size(); i++){
                beadIndexR2[1] = lpath2.brokenWorldlinesR[i];
                lpath2.updateBead(beadIndexR2,oldTailPos2[i]);
            }
            S = rho0Swap*rho0Swap2*exp(-(0.5)*(newPotAction+newPotAction2)+(0.5)*(oldPotAction+oldPotAction2));
            Z = rho0Dir*rho0Dir2;
            
            nBin = lpath.brokenWorldlinesR.size();
            if( nBin < constants()->initialNumParticles() )
                nMatch = true;
        }
    }else{
        S = 0.0;
        Z = rho0Dir*rho0Dir2;
    }
    
    /* In the normalization factor we must account for the factor of 1/2 in the
     open path weight */
    //Z = rho0Dir*rho0Dir2*exp(-(oldPotAction+oldPotAction2));
    //S = rho0Swap*rho0Swap2*exp(-(newPotAction+newPotAction2));
    
    //cout << rho0Dir << '\t' << rho0Dir2 << '\t' << Z << '\t' << S << endl;
    
    estimator(0) += Z;
    if( nMatch ){
        estimator(nBin*2-1) += Z;
        estimator(nBin*2) += S;
    }
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// END PIGS ESTIMATORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
