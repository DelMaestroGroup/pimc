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

/* This no longer seems to be needed.  Commenting out for now. */
//#ifndef GPU_BLOCK_SIZE
//    #include "special_functions.h"
//#endif

#ifdef GPU_BLOCK_SIZE
    #ifndef USE_CUDA
        #include "estimator.hip.h"
    #endif
    #ifdef USE_CUDA
        #include "estimator.cuh"
    #endif
#endif

/**************************************************************************//**
 * Setup the estimator factory.
******************************************************************************/
EstimatorFactory estimatorFactory;
#define REGISTER_ESTIMATOR(NAME,TYPE) \
    const string TYPE::name = NAME;\
    bool reg ## TYPE = estimatorFactory()->Register<TYPE>(TYPE::name); 

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
REGISTER_ESTIMATOR("commensurate order parameter",CommensurateOrderParameterEstimator);
REGISTER_ESTIMATOR("bipartition density",BipartitionDensityEstimator);
REGISTER_ESTIMATOR("linear density rho",LinearParticlePositionEstimator);
REGISTER_ESTIMATOR("planar density rho",PlaneParticlePositionEstimator);
REGISTER_ESTIMATOR("planar density average rho",PlaneParticleAveragePositionEstimator);
REGISTER_ESTIMATOR("planar potential average Vext",PlaneAverageExternalPotentialEstimator);
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
REGISTER_ESTIMATOR("cylinder static structure factor",CylinderStaticStructureFactorEstimator);
REGISTER_ESTIMATOR("pigs kinetic energy",KineticEnergyEstimator);
REGISTER_ESTIMATOR("pigs total energy",TotalEnergyEstimator);
REGISTER_ESTIMATOR("pigs thermodynamic potential energy",ThermoPotentialEnergyEstimator);
REGISTER_ESTIMATOR("pigs positions",PositionEstimator);
REGISTER_ESTIMATOR("pigs particle resolved positions",ParticleResolvedPositionEstimator);
REGISTER_ESTIMATOR("pigs particle correlations",ParticleCorrelationEstimator);
REGISTER_ESTIMATOR("pigs velocity",VelocityEstimator);
REGISTER_ESTIMATOR("pigs subregion occupation",SubregionOccupationEstimator);
REGISTER_ESTIMATOR("pigs one body density matrix",PIGSOneBodyDensityMatrixEstimator);

/* GPU accelerated estimators */
#ifdef GPU_BLOCK_SIZE
REGISTER_ESTIMATOR("intermediate scattering function gpu",IntermediateScatteringFunctionEstimatorGpu);
REGISTER_ESTIMATOR("static structure factor gpu",StaticStructureFactorGPUEstimator);
#endif

/**************************************************************************//**
 * Setup the estimator factory for multi path estimators.
******************************************************************************/
MultiEstimatorFactory multiEstimatorFactory;
const string SwapEstimator::name = "pigs multi swap";
bool regSwap = multiEstimatorFactory()->Register<SwapEstimator>(SwapEstimator::name);

const string EntPartEstimator::name = "pigs multi entanglement of particles";
bool regEntPart = multiEstimatorFactory()->Register<EntPartEstimator>(EntPartEstimator::name);


/* These functions recursively generate all vectors {(-max[0], ..., -max[NDIM-1]), ..., (-max[0], ..., -max[NDIM-1])} 
 * and returns the list */
/* void generateVectors(vector<int>& current, const iVec &max, unsigned int pos, vector<vector<int>>&results) { */
/*     if (pos == max.size()) { */
/*         results.push_back(current); */
/*         return; */
/*     } */

/*     for (int i = -max[pos]; i <= max[pos]; ++i) { */
/*         current[pos] = i; */
/*         generateVectors(current, max, pos + 1, results); */
/*     } */
/* } */

/* vector<vector<int>> generateAllVectors(const iVec &max) { */
/*     vector<int> current(max.size(), 0); */

/*     vector<vector<int>> results; */
/*     generateVectors(current, max, 0, results); */
/*     return results; */
/* } */


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

    /* A normalization factor for time slices used in estimators */
    sliceFactor.resize(constants()->numTimeSlices(),1.0);

    /*  The slices over which we average quantities PIGS vs. PIMC */
#if PIGS
    int midSlice = (path.numTimeSlices-1)/2;
    startSlice = midSlice - actionPtr->period;
    endSlice = midSlice + actionPtr->period;
    endDiagSlice = endSlice;

    if(actionPtr->local) {
        endDiagSlice = endSlice + 1;
        sliceFactor[startSlice] = 0.5;
        sliceFactor[endSlice] = 0.5;
    }
#else
    startSlice = 0;
    endSlice = path.numTimeSlices;
    endDiagSlice = endSlice;
#endif
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
 *  Initialize estimator. 
 *
 *  Initilize the estimator and normalization arrays.
******************************************************************************/
void EstimatorBase::initialize(vector<string> estLabel) {

    /* Initialize the map linking names to indices */
    for(size_t i = 0; i < estLabel.size(); ++i) 
        estIndex[estLabel[i]] = i;

    /* write the header string */
    header = "";
    for (const auto& label : estLabel)
        header += str(format("%16s") % label);

    numEst = estLabel.size();
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

        /* Write the header to disk if we are not restarting or if this is
         * a new estimator. */
        if (!constants()->restart() || (!communicate()->file(label)->exists())) {
            /* Check to see if the header has already been written to. 
             * Required for scalar estimators combined into one file. */
            if (!communicate()->file(label)->prepared()) {
                header.replace(header.begin(),header.begin()+1,"#");
                communicate()->file(label)->prepare();
            }
             
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
 *  Output a flat estimator value to disk.
 *
 *  Instead of keeping the individual binned averages, here we reset the 
 *  output file and write the current average to disk.
******************************************************************************/
void EstimatorBase::outputFlat() {

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
 *  Output a histogram-style estimator value to disk.
 *
 *  Normalization works differently here and we need to make sure there is at
 *  least 1 measurement per histogram bin.
******************************************************************************/
void EstimatorBase::outputHist() {

    /* Now write the estimator to disk */
    for (int n = 0; n < numEst; n++) { 
        if (abs(norm(n)) > 0.0)
            (*outFilePtr) << format("%16.8E") % (estimator(n)/norm(n));
        else
            (*outFilePtr) << format("%16.8E") % 0.0;
    }

    if (endLine)
        (*outFilePtr) << endl;

    /* Reset all values */
    norm = 0.0;
    reset();
}

/*************************************************************************//**
*  AppendLabel
******************************************************************************/
void EstimatorBase::appendLabel(string append) {
    label = label + append;
}

/*************************************************************************//**
* create a "(x,y,z)" string from a dVec for outputting 
******************************************************************************/
string EstimatorBase::dVecToString(const dVec& v){
    string strVec = "(";
    for (int i = 0; i < NDIM; i++) {
        strVec += str(format("%+15.8E") % v[i]);
        if (i < NDIM-1)
            strVec += ",";
    }
    return strVec + ")";
}

/*************************************************************************//**
*  Get q-vectors for scattering calculations
*  
*  Based on the geometry of the system and a user-defined choice, 
*  (wavevector and wavevector_type command line arguments) return a
*  list of q-vectors where scattering will be computed.
*  
******************************************************************************/
void EstimatorBase::getQVectors(std::vector<dVec> &qValues) {

    dVec q;

    /* Here we get the text from the command line as tokens */ 
    string input = constants()->wavevector();
    string inputType = constants()->wavevectorType();
    char *cstr = new char [input.length()+1];
    std::strcpy (cstr, input.c_str());
    char *token = std::strtok(cstr, " ");

    /* Below we go through the various possible command line options to 
     * create q-vectors and create the list */

    /* Here we have a NDIM-list of q-vectors indexed by integers, 
     * i.e. "1 0 1" or floats, i.e. "0.3 0 0.3" */
    if ((inputType == "int") || (inputType == "float")) {
        while(token) {
            int j = 0;
            for (int i=0; (i<NDIM) && token; i++) {
                if (inputType == "int")
                    q[i] = (2.0*M_PI/path.boxPtr->side[i])*std::atoi(token);
                else
                    q[i] = std::atof(token);

                token = strtok(NULL, " ");
                j = i;
            }
            if (j==(NDIM-1)){
                qValues.push_back(q);
            }
        }
    }

    if ((inputType == "max_int") || (inputType == "max_float")) {
        iVec q_max_int, _q_int;
        dVec q_max;

        while(token) {
            for (int i=0; (i<NDIM) && token; i++) {
                if (inputType == "max_int") {
                    q_max_int[i] = std::abs(std::atoi(token));
                    q_max[i] = q_max_int[i]*2.0*M_PI/path.boxPtr->side[i];
                }
                else {
                    q_max[i] = std::atof(token);
                    q_max_int[i] = 1 + static_cast<int>(q_max[i]*path.boxPtr->side[i]/2.0/M_PI);
                }
                token = strtok(NULL, " ");
            }
        }

        double q_mag_max = sqrt(dot(q,q));
        double q_mag;

        int n_q = 1;
        /* these are the negative maximal wavevectors */
        for (int i = 0; i < NDIM; i++) {
            n_q *= 2*q_max_int[i] + 1;
            _q_int[i] = -q_max_int[i];
            q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
        }

        q_mag = sqrt(dot(q,q));
        if (q_mag <= q_mag_max+EPS)
            qValues.push_back(q);

        int pos = NDIM - 1;
        int count = 0;
        while (count < n_q - 1) {
            if (_q_int[pos] == q_max_int[pos]) {
                _q_int[pos] = -q_max_int[pos];
                pos -= 1;
            } else {
                _q_int[pos] += 1;
                for (int i = 0; i < NDIM; i++) {
                    q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
                }

                q_mag = sqrt(dot(q,q));
                if (q_mag <= q_mag_max) 
                    qValues.push_back(q);
                count += 1;
                pos = NDIM - 1; //increment the innermost loop
            }
        }
    }

    if (constants()->wavevectorType() == "file_int") {

        std::ifstream file(input);
        std::string line;
        // Read one line at a time into the variable line:
        while(std::getline(file, line)) {
            std::vector<int> line_data;
            std::stringstream line_stream(line);
        
            int value;
            // Read an integer at a time from the line
            while(line_stream >> value) {
                // Add the integers from a line to a 1D array (vector)
                line_data.push_back(value);
            }
            PIMC_ASSERT(line_data.size()==NDIM);

            for (int i=0; i < NDIM; i++) {
                q[i] = (2.0*M_PI/path.boxPtr->side[i])*line_data[i];
            }
            qValues.push_back(q);
        }
    }
    if (constants()->wavevectorType() == "file_float") {

        std::ifstream file(input);
        std::string line;
        // Read one line at a time into the variable line:
        while(std::getline(file, line)) {
            std::vector<int> line_data;
            std::stringstream line_stream(line);
        
            int value;
            // Read an integer at a time from the line
            while(line_stream >> value) {
                // Add the integers from a line to a 1D array (vector)
                line_data.push_back(value);
            }
            PIMC_ASSERT(line_data.size()==NDIM);

            for (int i=0; i < NDIM; i++) {
                q[i] = (2.0*M_PI/path.boxPtr->side[i])*line_data[i];
            }
            qValues.push_back(q);
        }
    }

}

/*************************************************************************//**
*  Get q-vectors for scattering calculations
*  
*  Based on the geometry of the system and a user-defined choice, 
*  (wavevector and wavevector_type command line arguments) return a
*  list of q-vectors where scattering will be computed.
*  
******************************************************************************/
vector <vector<dVec> > EstimatorBase::getQVectors2(double dq, double qMax, 
        int& numq, string qGeometry) {

    /* initilize the total number of q-vectors */
    numq = 0;

    /* The q-vectors will end up in this array to be returned by value. */
    vector <vector<dVec> > q;

    double dtheta;
    if (qGeometry == "line") 
        dtheta = M_PI;
    else if (qGeometry == "sphere") {

        /* Number of θ values per q-magnitude, hard-coded for now */
        int numTheta = 24; 
        dtheta = 0.5*M_PI/numTheta;
    } 
    else {
        cerr << "\nERROR: A valid geometry wasn't chosen for q-space."
            << endl << "Action: choose \"line\" or \"sphere\"" << endl;
        exit(0);
    }

    /* Determine the set of q-vectors that have these magnitudes.  */
    for (double cq = 0.0; cq <= qMax + EPS; cq += dq)
    {
        vector <dVec> qvecs;

        /* cq = 0.0 */
        if (abs(cq) < EPS) {
            dVec qd = 0.0;
            qvecs.push_back(qd);
        }

        /* cq > 0.0 */
        else {

            /* First do θ = 0, i.e. along the z-direction */
            dVec qd = 0.0;
            qd[NDIM-1] = cq;
            qvecs.push_back(qd);

/* Can only do a spherical distribution of q-vectors in 3 spatial dimensions */
#if NDIM==3
            /* Now do the rest of the θ values */
            for (double theta = dtheta; theta <= 0.5*M_PI + EPS; theta += dtheta) {
                double dphi = dtheta/sin(theta);
                for (double phi = 0.0; phi <= 0.5*M_PI+EPS; phi += dphi) {

                    dVec qd;
                    qd[0] = cq*sin(theta)*cos(phi);
                    qd[1] = cq*sin(theta)*sin(phi);
                    qd[2] = cq*cos(theta);

                    qvecs.push_back(qd);
                } // phi
            } // theta
#endif 

        } // non-zero q-mags

        /* Add the list of q-vectors at this magnitude */
        numq += qvecs.size();
        q.push_back(qvecs);
    } //q-mags

    /* output */
    /* int totalNumQVecs = 0; */
    /* for (auto [nq,cq] : enumerate(q)) { */
    /*     double qMag = sqrt(blitz::dot(cq.front(),cq.front())); */
    /*     cout << endl << endl << "qmag = " << qMag << endl; */
    /*     for (const auto &cqvec : cq) */ 
    /*         cout << cqvec << endl; */
    /*     totalNumQVecs += cq.size(); */
    /* } */
    /* cout << "Full: " << totalNumQVecs << endl; */
    /* exit(-1); */
    cout << "numQ = " << numq << endl;

    return q;
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

    endLine = false;
    initialize({"us","mcsteps"});
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

    endLine = false;
    initialize({"K","V","V_ext","V_int","E","E_mu","K/N","V/N","E/N"});

    numPPAccumulated = 0;
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
EnergyEstimator::~EnergyEstimator() { 
}

/*************************************************************************//**
 *  Accumluate the energy.
 *
 *  We use the thermodynamic estimator for the kinetic and potential energy.
 *  The operator estimator is used for the potential energy. 
 *  A possible shift is made due to a tail correction.
******************************************************************************/
void EnergyEstimator::accumulate() {

    double totK = 0.0;
    double totV = 0.0;
    blitz::TinyVector<double,2> totVop(0.0);

    int numParticles  = path.getTrueNumParticles();
    int numTimeSlices = endSlice - startSlice;

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
    for (int slice = startSlice; slice < endSlice; slice++) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            vel = path.getVelocity(beadIndex);
            totK -= dot(vel,vel);
        }
    }
    /* Normalize the accumulated link-action part */
    totK *= kinNorm;

    /* Now we compute the potential and kinetic energy using the thermodynamic estimator.
     * The operator estimator is also computed on physical time slices. */ 
    double t1 = 0.0;
    double t2 = 0.0;

    for (int slice = startSlice; slice < endDiagSlice; slice++) {
        t1 += sliceFactor[slice]*actionPtr->derivPotentialActionLambda(slice);
        t2 += sliceFactor[slice]*actionPtr->derivPotentialActionTau(slice);
        if (!(slice % actionPtr->period))
            totVop  += sliceFactor[slice]*actionPtr->potential(slice);
    }

    t1 *= constants()->lambda()/(constants()->tau()*numTimeSlices);
    t2 /= 1.0*numTimeSlices;

    /* Normalize the action correction and the total potential*/
    totVop /= (numTimeSlices/actionPtr->period);

    /* Perform all the normalizations and compute the individual energy terms */
    totK += (classicalKinetic + t1);
    totV = t2 - t1 + tailV;
    totVop[1] += tailV;

    /* Now we accumulate the average total, kinetic and potential energy, 
     * as well as their values per particles. */
    estimator(estIndex["K"]) += totK;
    estimator(estIndex["V"]) += totV;
    estimator(estIndex["V_ext"]) += totVop[0];
    estimator(estIndex["V_int"]) += totVop[1];

    estimator(estIndex["E"]) += totK + totV;
    estimator(estIndex["E_mu"]) += totK + totV - constants()->mu()*numParticles;

    /* Measure "per particle" estimators in the GCE */
    if (numParticles > 0) {
        numPPAccumulated += 1;

        estimator(estIndex["K/N"]) += totK/(1.0*numParticles);
        estimator(estIndex["V/N"]) += totV/(1.0*numParticles);
        estimator(estIndex["E/N"]) += (totK + totV)/(1.0*numParticles);
    }

    /* We have to modify the normalization in this case as we don't 
     * necessrily measure every time accumulate is called */
    if (numAccumulated == constants()->binSize()) {
        norm(estIndex["K/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["V/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["E/N"]) = 1.0*numAccumulated/numPPAccumulated;

        numPPAccumulated = 0;
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
    /* header = str(format("#%15s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s") */ 
    /*         % "K_op" % "K_cv" % "V_op" % "V_cv" % "E" % "E_mu" % "K_op/N" % "K_cv/N" % "V_op/N" */
    /*         % " V_cv/N" % "E/N" % "EEcv*Beta^2"% "Ecv*Beta" % "dEdB" % "CvCov1" */
    /*         % "CvCov2" % "CvCov3" % "E_th" % "P"); */
    /* initialize(19); */

    initialize({"K_op","K_cv","V_op","V_cv","E","E_mu","K_op/N","K_cv/N","V_op/N",
            "V_cv/N","E/N","EEcv*Beta^2","Ecv*Beta","dEdB","CvCov1", "CvCov2","CvCov3",
            "E_th","P"});

    numPPAccumulated = 0;
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
            totVop  += sum(actionPtr->potential(slice));

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
    estimator(estIndex["K_op"]) += totEcv - totVop; // operator kinetic energy
    estimator(estIndex["K_cv"]) += Kcv; // centroid virial kinetic energy
    estimator(estIndex["V_op"]) += totVop; // operator potential energy
    estimator(estIndex["V_cv"]) += totEcv - Kcv; //centroid virial potential energy
    estimator(estIndex["E"]) += totEcv; // total energy
    estimator(estIndex["E_mu"]) += totEcv - constants()->mu()*numParticles;

    /* Measure "per particle" estimators */
    if (numParticles > 0) {
        numPPAccumulated += 1;

        estimator(estIndex["K_op/N"]) += (totEcv - totVop)/(1.0*numParticles);
        estimator(estIndex["K_cv/N"]) += Kcv/(1.0*numParticles);
        estimator(estIndex["V_op/N"]) += totVop/(1.0*numParticles);
        estimator(estIndex["V_cv/N"]) += (totEcv - Kcv)/(1.0*numParticles);
        estimator(estIndex["E/N"]) += totEcv/(1.0*numParticles);
    }

    /* We have to modify the normalization in this case as we don't 
     * necessrily measure every time accumulate is called */
    if (numAccumulated == constants()->binSize()) {
        norm(estIndex["K_op/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["K_cv/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["V_op/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["V_cv/N"]) = 1.0*numAccumulated/numPPAccumulated;
        norm(estIndex["E/N"]) = 1.0*numAccumulated/numPPAccumulated;

        numPPAccumulated = 0;
    }

    /* accumulate specific heat estimators. */
    estimator(estIndex["EEcv*Beta^2"]) += totEcv*thermE*beta*beta;
    estimator(estIndex["Ecv*Beta"]) += totEcv*beta;
    estimator(estIndex["dEdB"]) += dEdB;
    estimator(estIndex["CvCov1"]) += totEcv*thermE*beta*beta*totEcv*beta;
    estimator(estIndex["cVCov2"]) += totEcv*beta*dEdB;
    estimator(estIndex["CvCov3"]) += totEcv*thermE*beta*beta*dEdB;

    /* thermodynamic energy */
    estimator(estIndex["E_th"]) += thermE;

    /* Pressure */
    estimator(estIndex["P"]) += Pressure;
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
    endLine = false;
    initialize({"N","N^2","density"});
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
// COMMENSURATE ORDER PARAMETER ESTIMATOR CLASS ------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  We measure the average static structure factor at a finite-set of wave-
 *  vectors that correspond to the first shell of reciprocal lattice vectors.
 *
 *  @see https://journals.aps.org/prb/abstract/10.1103/PhysRevB.73.085422
 *  
******************************************************************************/
CommensurateOrderParameterEstimator::CommensurateOrderParameterEstimator (
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Get the current carbon-carbon distance */
    double aCC = constants()->aCC();

    /* The reciprocal lattice vectors*/
    dVec G1,G2;
    G1[0] = 2.0*M_PI/(sqrt(3.0)*aCC);
    G1[1] = 2.0*M_PI/(3.0*aCC);
    G1[2] = 0.0;

    G2[0] = -G1[0];
    G2[1] = G1[1];
    G2[2] = 0.0;

    /* For now we hard-code the g-vectors for graphene */
    g.push_back(G1);
    g.push_back(G2);
    g.push_back(G1+G2);

    /* Set estimator name and header */
    endLine = false;
    initialize({"Scom"});

    norm = 1.0/(g.size()*constants()->numTimeSlices());
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CommensurateOrderParameterEstimator::~CommensurateOrderParameterEstimator() { 
}

/*************************************************************************//**
 * Accumulate the number of Commensurate order parameter
******************************************************************************/
void CommensurateOrderParameterEstimator::accumulate() {

    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();

    double _norm = 1.0;

    if (numParticles > 0)
        _norm /= numParticles;

    beadLocator beadIndex;  // The bead locator
    double Scom = 0.0; // initialize
    dVec pos;

    /* Average over all time slices */
    for (beadIndex[0] = 0; beadIndex[0] < numTimeSlices; beadIndex[0]++) {

        /* Average over all beads */
        for (beadIndex[1] = 0; beadIndex[1] < path.numBeadsAtSlice(beadIndex[0]); beadIndex[1]++) {
            pos = path(beadIndex);

            for (const auto &cg : g) {
                Scom += cos(dot(cg,pos));
            }

        } // beadIndex[1]
    } //beadIndex[0]

    estimator(0) += Scom*_norm;
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
 *  A full NDIM-dimensional particle density histogram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
ParticlePositionEstimator::ParticlePositionEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label)  : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    initialize(path.boxPtr->numGrid);

    /* Setup estimator name and header. We include the PIMCID as we 
     * replace the file on every output.*/
    diffLabels = {"dx","dy","dz"};

    header =  str(format("# PIMCID: %s\n") % constants()->id());
    header += "# ESTINF:";
    for (int i = 0; i < NDIM; i++) 
        header += str(format(" %s = %12.6E") % diffLabels[i] % path.boxPtr->gridSize[i]);
    header += str(format(" NGRIDSEP = %d\n") % NGRIDSEP);
    header += str(format("#%15s") % "density");

    /* The normalization: 1/(dV*M) */
    for (int n = 0; n < numEst; n++)
        norm(n) = 1.0/(1.0*(endSlice-startSlice)*(1.0/actionPtr->period) *
                path.boxPtr->gridBoxVolume(n));
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
ParticlePositionEstimator::~ParticlePositionEstimator() { 
}

/*************************************************************************//**
 *  Accumulate a histogram of all particle positions, with output 
 *  being the running average of the density per grid space.
******************************************************************************/
void ParticlePositionEstimator::accumulate() {

    beadLocator beadIndex;

    for (int slice = startSlice; slice < endDiagSlice; slice += actionPtr->period) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;

            /* update our particle position histogram */
            int n = path.boxPtr->gridIndex(path(beadIndex));
            estimator(n) += sliceFactor[slice];
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
    blitz::Array<double,1> excLens (actionPtr->externalPtr->getExcLen());
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
// LINEAR PARTICLE DENSITY ESTIMATOR CLASS ------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  A 1-dimensional particle density histogram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
LinearParticlePositionEstimator::LinearParticlePositionEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* Make sure we have a bin centered at 0.0 */
    numGrid = 2*NGRIDSEP-1;
    dz = path.boxPtr->side[NDIM-1] / numGrid;

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* A local copy of the box size */
    side = path.boxPtr->side;

    /* The header contains information about the grid  */
    header = str(format("# ESTINF: dz = %12.6E NGRIDSEP = %d\n") % dz % numGrid);
    header += str(format("#%15.6E") % (-0.5*side[NDIM-1]));
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.6E") % (n*dz - 0.5*side[NDIM-1])));

    double A = 1.0;
    for (int i = 0; i < NDIM-1; i++)
        A *= side[i];

    norm = 1.0/(1.0*(endSlice-startSlice)*(1.0/actionPtr->period)*A*dz);
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
LinearParticlePositionEstimator::~LinearParticlePositionEstimator() { 
}

/*************************************************************************//**
 *  Accumulate a histogram of all particle positions, with output 
 *  being the running average of the density per grid space.
******************************************************************************/
void LinearParticlePositionEstimator::accumulate() {

    beadLocator beadIndex;
    dVec pos;
    int index;

    for (int slice = startSlice; slice < endDiagSlice; slice += actionPtr->period) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            /* Get the z-index */
            index = static_cast<int>(abs(pos[NDIM-1] + 0.5*side[NDIM-1] - EPS ) / (dz + EPS));

            /* update our particle position histogram */
            if (index < numGrid)
                estimator(index) += sliceFactor[slice];
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PLANE PARTICLE DENSITY ESTIMATOR CLASS ------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  A 2-dimensional particle density histogram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
PlaneParticlePositionEstimator::PlaneParticlePositionEstimator (const Path &_path, 
        ActionBase *_actionPtr, const MTRand &_random, double _maxR, 
        int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {


    /* We choose an odd number to make sure (0,0) is the central 
     * grid box. */
    numLinearGrid = 4*NGRIDSEP+1;
    numGrid = numLinearGrid*numLinearGrid;

    /* The spatial discretization */
    for (int i = 0; i < NDIM; i++)
        dl[i]  = path.boxPtr->side[i] / numLinearGrid;

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header contains information about the grid  */
    header = str(format("# ESTINF: dx = %12.6E dy = %12.6E NGRIDSEP = %d\n") 
            % dl[0] % dl[1] % numLinearGrid);
    header += str(format("#%15.3E") % 0.0);
    for (int n = 1; n < numGrid; n++) 
        header.append(str(format("%16.3E") % (1.0*n)));

    /* Compute the area of a grid box */
    double A = 1.0;
    for (int i = 0; i < NDIM-1; i++)
        A *= dl[i];

    norm = 1.0/((endSlice-startSlice)*(1.0/actionPtr->period)*A*path.boxPtr->side[NDIM-1]);
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

    for (int slice = startSlice; slice < endDiagSlice; slice += actionPtr->period) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            int index = 0;
            for (int i = 0; i < NDIM-1; i++) {  
                int scale = 1;
                for (int j = i+1; j < NDIM-1; j++) 
                    scale *= numLinearGrid;
                index += scale*static_cast<int>(abs(pos[i] + 0.5*side[i] - EPS ) / (dl[i] + EPS));
            }

            /* update our particle position histogram */
            if (index < numGrid)
                estimator(index) += sliceFactor[slice];
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PLANE PARTICLE AVERAGE DENSITY ESTIMATOR CLASS ----------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  A 2-dimensional particle density averaged histogram.
 *
 *  @note Only tested for cubic boxes
******************************************************************************/
PlaneParticleAveragePositionEstimator::PlaneParticleAveragePositionEstimator (
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* We choose an odd number to make sure (0,0) is the central 
     * grid box. */
    numLinearGrid = 4*NGRIDSEP+1;
    numGrid = numLinearGrid*numLinearGrid;

    /* The spatial discretization */
    for (int i = 0; i < NDIM; i++)
        dl[i]  = path.boxPtr->side[i] / numLinearGrid;

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header contains information about the grid  */
    header = str(format("# PIMCID: %s\n") % constants()->id());
    header += str(format("# ESTINF: dx = %12.6E dy = %12.6E NGRIDSEP = %d\n") 
            % dl[0] % dl[1] % numLinearGrid);
    header += str(format("#%15s") % "plane density");

    /* Compute the area of a grid box */
    double A = 1.0;
    for (int i = 0; i < NDIM-1; i++)
        A *= dl[i];

    norm = 1.0/((endSlice-startSlice)*(1.0/actionPtr->period)*A*path.boxPtr->side[NDIM-1]);
    side = path.boxPtr->side;
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
PlaneParticleAveragePositionEstimator::~PlaneParticleAveragePositionEstimator() { 
}

/*************************************************************************//**
 *  Accumulate a histogram of all particle positions, with output 
 *  being the running average of the density per grid space.
******************************************************************************/
void PlaneParticleAveragePositionEstimator::accumulate() {

    beadLocator beadIndex;
    dVec pos;

    for (int slice = startSlice; slice < endDiagSlice; slice += actionPtr->period) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            int index = 0;
            for (int i = 0; i < NDIM-1; i++) {  
                int scale = 1;
                for (int j = i+1; j < NDIM-1; j++) 
                    scale *= numLinearGrid;
                index += scale*static_cast<int>(abs(pos[i] + 0.5*side[i] - EPS ) / (dl[i] + EPS));
            }

            /* update our particle position histogram */
            if (index < numGrid)
                estimator(index) += sliceFactor[slice];
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// PLANE AVERAGE EXTERNAL POTENTIAL ESTIMATOR
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  2-dimensional average external potential.
 *
******************************************************************************/
PlaneAverageExternalPotentialEstimator::PlaneAverageExternalPotentialEstimator (
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) : 
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    /* We choose an odd number to make sure (0,0) is the central 
     * grid box. */
    numLinearGrid = 4*NGRIDSEP+1;
    numGrid = numLinearGrid*numLinearGrid;

    /* The spatial discretization */
    for (int i = 0; i < NDIM; i++)
        dl[i]  = path.boxPtr->side[i] / numLinearGrid;

    /* This is a diagonal estimator that gets its own file */
    initialize(numGrid);

    /* The header contains information about the grid  */
    header = str(format("# PIMCID: %s\n") % constants()->id());
    header += str(format("# ESTINF: dx = %12.6E dy = %12.6E NGRIDSEP = %d\n") 
            % dl[0] % dl[1] % numLinearGrid);
    header += str(format("#%15s") % "plane external potential");

    side = path.boxPtr->side;
}

/*************************************************************************//**
 *  Destructor
******************************************************************************/
PlaneAverageExternalPotentialEstimator::~PlaneAverageExternalPotentialEstimator() { 
}

/*************************************************************************//**
 *  Accumulate the average external potential energy felt by all particles 
 *  as a function of x and y.
******************************************************************************/
void PlaneAverageExternalPotentialEstimator::accumulate() {

    beadLocator beadIndex;
    dVec pos;

    for (int slice = startSlice; slice < endDiagSlice; slice += actionPtr->period) {
        for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
            beadIndex = slice,ptcl;
            pos = path(beadIndex);

            /* Obtain the index of the particle position */
            int index = 0;
            for (int i = 0; i < NDIM-1; i++) {  
                int scale = 1;
                for (int j = i+1; j < NDIM-1; j++) 
                    scale *= numLinearGrid;
                index += scale*static_cast<int>(abs(pos[i] + 0.5*side[i] - EPS ) / (dl[i] + EPS));
            }

            /* update the external potential in each grid box */
            if (index < numGrid) {
                norm(index) += sliceFactor[slice];
                estimator(index) += sliceFactor[slice]*actionPtr->externalPtr->V(pos);
            }
        } //ptcl
    } //slice
}

/*************************************************************************//**
 *  Output a flat estimator value to disk.
 *
 *  Instead of keeping the individual binned averages, here we reset the 
 *  output file and write the current average to disk.
******************************************************************************/
void PlaneAverageExternalPotentialEstimator::output() {

    /* Prepare the position file for writing over old data */
    communicate()->file(label)->reset();

    (*outFilePtr) << header;
    if (endLine)
        (*outFilePtr) << endl;

    /* Now write the running average of the estimator to disk */
    double Vext = 0.0;
    for (int n = 0; n < numEst; n++) { 
        if (abs(norm(n)) > 0.0)
            Vext = estimator(n)/norm(n);
        else
            Vext = 0.0;
        (*outFilePtr) << format("%16.8E\n") % Vext;
    }

    communicate()->file(label)->rename();
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
        header += str(format("%16s") % str(format("P(%+2d)") % w));
    header += str(format("%16s") % "Area_rho_s");

    /* The pre-factor for the superfluid density is always the same */
    W2Norm = constants()->T() / (2.0 * sum(path.boxPtr->periodic) * constants()->lambda());

    /* The pre-factor for the area esimator */
    ANorm = 0.5*constants()->T()*constants()->numTimeSlices()/constants()->lambda();
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

    if (path.getTrueNumParticles() > 0) {
        numPPAccumulated += 1;

        /* The average winding number squared */
        estimator(0) += locW2oN;

        /* The Area Estimator */
        estimator(5+2*windMax) += Az*Az/I;
    }

    /* We need to re-normalize as it may not be measured everytime
     * for small numbers of particles */
    if (numAccumulated == constants()->binSize()) {
        norm(0) = W2Norm*numAccumulated/numPPAccumulated;
        norm(5+2*windMax) = ANorm*numAccumulated/numPPAccumulated;
        numPPAccumulated = 0;
    }

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
/* void LocalPermutationEstimator::output() { */

/*     /1* Prepare the position file for writing over old data *1/ */
/*     communicate()->file(label)->reset(); */

/*     (*outFilePtr) << header; */
/*     if (endLine) */
/*         (*outFilePtr) << endl; */

/*     /1* Now write the running average of the estimator to disk *1/ */
/*     for (int n = 0; n < numEst; n++) */
/*         (*outFilePtr) << format("%16.6E\n") % */
/*             (estimator(n)/totNumAccumulated); */

/*     communicate()->file(label)->rename(); */
/* } */

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
    if (path.boxPtr->name == "XXX") {
        for (int n = 0; n < NPCFSEP; n++)
            norm(n) = 0.5*path.boxPtr->side[NDIM-1] / dR;
    }
    else {
	blitz::TinyVector<double,3> gNorm;
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
    double tot = 1.0*sum(actionPtr->sepHist);
    if ((numParticles > 1) && (tot > 0.0)) 
        estimator += lnorm*(1.0*actionPtr->sepHist / tot );
    else
        estimator += 0.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// STATIC STRUCTURE FACTOR ESTIMATOR CLASS -----------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
//StaticStructureFactorEstimator::StaticStructureFactorEstimator(
//        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
//        double _maxR, int _frequency, string _label) :
//    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
//{
//
//    /* hard-coded number of q-vector magnitudes */
//    numq = 20; 
//
//    /* initialize the number of vectors with each magnitude */
//    numqVecs.resize(numq);               
//    numqVecs = 0;
//
//    /* initialize the magnitude array */
//    qMag.resize(numq);
//    qMag = 0.0;
//
//    /* The maximum q-vector magnitude to consider (hard-coded for now) */
//    double qmax = 4.0; // 1/A
//    double dq = qmax/(numq-1);
//
//    /* Determine the set of magnitudes */
//    Array <double, 1> qMag(numq);         // the q-vector magnitudes
//    for (int nq = 0; nq < numq; nq++)
//        qMag(nq) = nq*dq;
//
//    /* Determine the set of q-vectors that have these magintudes */
//    for (int nq = 0; nq < numq; nq++)
//    {
//        double cq = qMag(nq);
//        vector <dVec> qvecs;
//
//        int maxComp = ceil(cq*blitz::max(path.boxPtr->side)/(2.0*M_PI))+1;
//        int maxNumQ = ipow(2*maxComp + 1,NDIM);
//        
//        iVec qi;
//        for (int n = 0; n < maxNumQ; n++) {
//            for (int i = 0; i < NDIM; i++) {
//                int scale = 1;
//                for (int j = i+1; j < NDIM; j++) 
//                    scale *= (2*maxComp + 1);
//                qi[i] = (n/scale) % (2*maxComp + 1);
//
//                /* Wrap into the appropriate set of integers*/
//                qi[i] -= (qi[i] > maxComp)*(2*maxComp + 1);
//            }
//
//            dVec qd = 2.0*M_PI*qi/path.boxPtr->side;
//
//            /* Store the winding number */
//            double cqMag = sqrt(dot(qd,qd));
//            if ( ( (abs(cq) < EPS) && (cqMag < EPS) ) ||
//                 ( (abs(cq) >= EPS) && (cqMag > (cq-dq)) && (cqMag <= cq) ) ) 
//            {
//                qvecs.push_back(qd);
//                numqVecs(nq)++;
//            }
//        }
//        
//        q.push_back(qvecs);
//    }
//
//

/*************************************************************************//**
 *  Constructor.
 * 
 *  Compute the static structure factor for a fixed set of q vector magnitude.
******************************************************************************/
StaticStructureFactorEstimator::StaticStructureFactorEstimator(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    /* The maximum q-vector magnitude to consider (hard-coded for now) */
    double qMax = 4.0; // 1/Å

    /* We choose dq from the smallest possible q-vector, set by PBC */
    double dq = 2.0*M_PI/path.boxPtr->side[NDIM-1];

    /* Get the desired q-vectors */
    int numq = 0;
    q = getQVectors2(dq,qMax,numq,"sphere");

    /* Determine how many q-vector magnitudes  we have */
    numq = q.size();

    /* Initialize the accumulator intermediate scattering function*/
    sf.resize(numq);
    sf = 0.0;

    /* This is a diagonal estimator that gets its own file */
    initialize(numq);

    /* The magnitude of q */
    header = str(format("#%15.6E") % 0.0);
    for (int nq = 1; nq < numq; nq++)  {
        double qMag = sqrt(blitz::dot(q[nq][0],q[nq][0]));
        header.append(str(format("%16.6E") % (qMag)));
    }
        /* header.append(str(format("%16.6E") % (qMag(nq)-0.5*dq))); */

    /* Utilize imaginary time translational symmetry */
    norm = 1.0/constants()->numTimeSlices();

    /* Normalize by the number of q-vecs (except when there are none) */
    for (int nq = 0; nq < numq; nq++) {
        if (q[nq].size() > 0)
            norm(nq) /= q[nq].size();
    }
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
StaticStructureFactorEstimator::~StaticStructureFactorEstimator() { 
    sf.free();
}

/*************************************************************************//**
 *  Measure the static structure factor for each value for each q-magnitude
 *
******************************************************************************/
void StaticStructureFactorEstimator::accumulate() {

    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();

    beadLocator bead1,bead2;  // The bead locator
    sf = 0.0; // initialize

    /* q-magnitudes */
    for (auto [nq,cq] : enumerate(q)) {

        /* q-vectors with that q-magnitude*/
        for (const auto &cqvec : cq) {
    
            /* Average over all time slices */
            for (bead1[0] = 0; bead1[0] < numTimeSlices; bead1[0]++) {

                /* This is an equal imaginary time estimator */
                bead2[0] = bead1[0];

                for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(bead1[0]); bead1[1]++) {

                    /* The bead1 = bead2 part */
                    sf(nq) += 1.0;

                    for (bead2[1] = bead1[1]+1; bead2[1] < path.numBeadsAtSlice(bead2[0]); bead2[1]++) {
                        sf(nq) += 2*cos(dot(path.getSeparation(bead1,bead2),cqvec));
                    } // bead2[1]
                } // bead1[1]
            } //bead1[0]

        } // q-vectors
    } // q-magnitudes

    estimator += sf/numParticles; 
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// STATIC STRUCTURE FACTOR GPU ESTIMATOR CLASS -------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/*************************************************************************//**
 *  Constructor.
 * 
 *  A GPU accelerated static structure factor estimator.
 *  
******************************************************************************/
#ifdef GPU_BLOCK_SIZE
StaticStructureFactorGPUEstimator::StaticStructureFactorGPUEstimator(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{

    /* The maximum q-vector magnitude to consider (hard-coded for now) */
    double qMax = 4.0; // 1/Å

    /* We choose dq from the smallest possible q-vector, set by PBC */
    double dq = 2.0*M_PI/path.boxPtr->side[NDIM-1];

    /* Get the desired q-vectors */
    q = getQVectors2(dq,qMax,numq,"sphere");

    /* For this estimator we measure the scattering at each vector q-values
     * so we need to flatten the magnitude ordered list */
    /* NOTE: I can probably do away with this and simply flatten the 2D array */
    qValues_dVec.resize(numq);
    int nq = 0;
    for (const auto &cq : q) {
        for (const auto &cqvec : cq) {
            qValues_dVec(nq) = cqvec;
            nq += 1;
        }
    }

    /* Initialize the accumulator for the static structure factor */
    ssf.resize(numq);
    ssf = 0.0;

    // Create multiple gpu streams
    stream_array.resize(MAX_GPU_STREAMS);
    for (int i = 0; i < MAX_GPU_STREAMS; i++) {
        #ifndef USE_CUDA
        HIP_ASSERT(hipStreamCreate(&stream_array(i)));
        #endif
        #ifdef USE_CUDA
        CUDA_ASSERT(cudaStreamCreate(&stream_array(i)));
        #endif
    }

    /* This is a diagonal estimator that gets its own file */
    initialize(numq);

    /* The header consists of an extra line of the possible q-values */
    header = str(format("# ESTINF: num_q = %d; ") % numq);
    for (int n = 0; n < numq; n++)
        header += dVecToString(qValues_dVec(n)) + " ";
    header += "\n";

    /* We index the q-vectors with an integer */
    header += str(format("#%15d") % 0);
    for (int n = 1; n < numq; n++) 
        header += str(format("%16d") % n);

    /* utilize imaginary time translational symmetry */
    norm = 0.5/constants()->numTimeSlices();

    bytes_beads = NDIM*(1 + constants()->initialNumParticles())*sizeof(double);
    bytes_ssf = ssf.size()*sizeof(double);
    bytes_qvecs = NDIM*numq*sizeof(double);
    #ifndef USE_CUDA
        HIP_ASSERT(hipMalloc(&d_ssf, bytes_ssf)); // Allocate memory for ssf on GPU
        HIP_ASSERT(hipMalloc(&d_qvecs, bytes_qvecs)); // Allocate memory for qvecs on GPU
        HIP_ASSERT(hipMemcpy(d_qvecs, qValues.data(), bytes_qvecs, hipMemcpyHostToDevice )); // Copy qvecs data to gpu
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaMalloc(&d_ssf, bytes_ssf)); // Allocate memory for ssf on GPU
        CUDA_ASSERT(cudaMalloc(&d_qvecs, bytes_qvecs)); // Allocate memory for qvecs on GPU
        CUDA_ASSERT(cudaMemcpy(d_qvecs, qValues_dVec.data(), bytes_qvecs, cudaMemcpyHostToDevice )); // Copy qvecs data to gpu
    #endif
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
StaticStructureFactorGPUEstimator::~StaticStructureFactorGPUEstimator() { 
    for (int i = 0; i < MAX_GPU_STREAMS; i++) {
        #ifndef USE_CUDA
            HIP_ASSERT(hipStreamDestroy(stream_array(i)));
        #endif
        #ifdef USE_CUDA
            CUDA_ASSERT(cudaStreamDestroy(stream_array(i)));
        #endif
    }
    ssf.free();
    stream_array.free();

    // Release device memory
    #ifndef USE_CUDA
        HIP_ASSERT(hipFree(d_beads));
        HIP_ASSERT(hipFree(d_qvecs));
        HIP_ASSERT(hipFree(d_ssf));
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaFree(d_beads));
        CUDA_ASSERT(cudaFree(d_qvecs));
        CUDA_ASSERT(cudaFree(d_ssf));
    #endif
}

/*************************************************************************//**
 *  Measure the static structure factor for each q-vector
 *
 *  We only compute this for N > 1 due to the normalization.
******************************************************************************/
void StaticStructureFactorGPUEstimator::accumulate() {

    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();

    /* Return to these and check if we need them */
    int number_of_beads = numParticles*numTimeSlices;
    int NNM = number_of_beads*numParticles;
    int beta_over_two_idx = numTimeSlices/2;

    double _inorm = 1.0/numParticles;

    /* We need to copy over the current beads array to the device */
    auto beads_extent = path.get_beads_extent();
    int full_number_of_beads = beads_extent[0]*beads_extent[1];
    int full_numTimeSlices = beads_extent[0];
    int full_numParticles = beads_extent[1];

    /* Size, in bytes, of beads array */
    size_t bytes_beads_new = NDIM*full_number_of_beads*sizeof(double);

    #ifndef USE_CUDA
        if (bytes_beads_new > bytes_beads) {
            bytes_beads = bytes_beads_new;
            HIP_ASSERT(hipFree(d_beads));
            HIP_ASSERT(hipMalloc(&d_beads, bytes_beads)); // Allocate memory for beads on GPU
        }
        HIP_ASSERT(hipMemcpy( d_beads, path.get_beads_data_pointer(), bytes_beads, hipMemcpyHostToDevice )); // Copy beads data to gpu
        HIP_ASSERT(hipMemset(d_ssf, 0, bytes_ssf)); // Set initial ssf data to zero
    #endif
    #ifdef USE_CUDA
        if (bytes_beads_new > bytes_beads) {
            bytes_beads = bytes_beads_new;
            CUDA_ASSERT(cudaFree(d_beads));
            CUDA_ASSERT(cudaMalloc(&d_beads, bytes_beads)); // Allocate memory for beads on GPU
        }
        CUDA_ASSERT(cudaMemcpy( d_beads, path.get_beads_data_pointer(), bytes_beads, cudaMemcpyHostToDevice )); // Copy beads data to gpu
        CUDA_ASSERT(cudaMemset(d_ssf, 0, bytes_ssf)); // Set initial ssf data to zero
    #endif

    int grid_size = (NNM + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;

    int stream_idx;
    for (int nq = 0; nq < numq; nq++) {
        int Mi = 0;
        /* stream_idx = (nq*numTimeSlices + Mi) % MAX_GPU_STREAMS; */
        stream_idx = (nq + Mi) % MAX_GPU_STREAMS;

        #ifndef USE_CUDA
        hipLaunchKernelGGL(gpu_isf, dim3(grid_size), dim3(GPU_BLOCK_SIZE), 0, stream_array(stream_idx),
                d_ssf, d_qvecs, d_beads, nq, Mi, _inorm, numq, numTimeSlices, numParticles,
                number_of_beads, full_numTimeSlices, full_numParticles, full_number_of_beads, NNM, beta_over_two_idx);
        #endif
        #ifdef USE_CUDA
        cuda_wrapper::gpu_isf_wrapper(dim3(grid_size), dim3(GPU_BLOCK_SIZE), stream_array(stream_idx),
                d_ssf, d_qvecs, d_beads, nq, Mi, _inorm, numq, numTimeSlices, numParticles,
                number_of_beads, full_numTimeSlices, full_numParticles, full_number_of_beads, NNM, beta_over_two_idx);
        #endif
    }
    #ifndef USE_CUDA
        HIP_ASSERT(hipDeviceSynchronize());
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaDeviceSynchronize());
    #endif

    //// Copy ssf data back to host
    #ifndef USE_CUDA
        HIP_ASSERT(hipMemcpy(ssf.data(), d_ssf, bytes_ssf, hipMemcpyDeviceToHost)); //Only copy up to beta/2 back to host
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaMemcpy(ssf.data(), d_ssf, bytes_ssf, cudaMemcpyDeviceToHost)); //Only copy up to beta/2 back to host
    #endif

    estimator += ssf;

}
#endif

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
    blitz::Array <double, 1> qMag(numq);         // the q-vector magnitudes
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

        /* Make sure we have some q-vectors */
        if (qvecs.size() < 1) {
            cerr << "\nERROR: Intermediate Scattering function: "
                 << "No valid q-vectors were added to the list for measurment." 
                 << endl << "Action: modify q-magintudes." << endl;
            exit(0);
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
// INTERMEDIATE SCATTERING FUNCTION GPU ESTIMATOR CLASS ----------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Measure the intermediate scattering function at wavevectors determined
 *  from the wavevector and wavevector_type command line arguments
******************************************************************************/
#ifdef GPU_BLOCK_SIZE
IntermediateScatteringFunctionEstimatorGpu::IntermediateScatteringFunctionEstimatorGpu(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) {

    int numTimeSlices = constants()->numTimeSlices();

    dVec q;
    string input = constants()->wavevector();
    char * cstr = new char [input.length()+1];
    std::strcpy (cstr, input.c_str());
    char *token = strtok(cstr, " ");
    // parse wavevectorType to determine how to handle isf_input
    if (constants()->wavevectorType() == "int") {
        std::cout << "wavevectorType = int" << std::endl;
        while(token) {
            int j = 0;
            for (int i=0; (i<NDIM) && token; i++) {
                q[i] = (2.0*M_PI/path.boxPtr->side[i])*std::atoi(token);
                token = strtok(NULL, " ");
                j = i;
            }
            if (j==(NDIM-1)){
                qValues.push_back(q);
            }
        }
    }

    if (constants()->wavevectorType() == "float") {
        std::cout << "wavevectorType = float" << std::endl;
        while(token) {
            int j = 0;
            for (int i=0; (i<NDIM) && token; i++){
                q[i] = std::atof(token);
                token = strtok(NULL, " ");
                j = i;
            }
            if (j==(NDIM-1)){
                qValues.push_back(q);
            }
        }
    }

    if (constants()->wavevectorType() == "max_int") {
        std::cout << "wavevectorType = max-int" << std::endl;
        iVec q_int;
        while(token) {
            for (int i=0; (i<NDIM) && token; i++) {
                q_int[i] = std::abs(std::atoi(token));
                token = strtok(NULL, " ");
            }
        }
        
        int _q_int[NDIM];
        int n_q = 1;
        for (int i = 0; i < NDIM; i++) {
            n_q *= 2*q_int[i] + 1;
            _q_int[i] = -q_int[i];
            q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
        }
        qValues.push_back(q);

        int pos = NDIM - 1;
        int count = 0;
        while (count < n_q - 1) {
            if (_q_int[pos] == q_int[pos]) {
                _q_int[pos] = -q_int[pos];
                pos -= 1;
            } else {
                _q_int[pos] += 1;
                for (int i = 0; i < NDIM; i++) {
                    q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
                }
                qValues.push_back(q);
                count += 1;
                pos = NDIM - 1; //increment the innermost loop
            }
        }
    }

    if (constants()->wavevectorType() == "max_float") {
        std::cout << "wavevectorType = max-float" << std::endl;
        iVec q_int;
        double q_mag_max = 0.0;
        double q_mag;
        while(token) {
            for (int i=0; (i<NDIM) && token; i++) {
                q_mag_max = std::atof(token);
                token = strtok(NULL, " ");
            }
        }
        
        int _q_int[NDIM];
        int n_q = 1;
        for (int i = 0; i < NDIM; i++) {
            q_int[i] = 1 + static_cast<int>(q_mag_max*path.boxPtr->side[i]/2.0/M_PI);
            n_q *= 2*q_int[i] + 1;
            _q_int[i] = -q_int[i];
            q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
        }
        q_mag = sqrt(dot(q,q));
        if (q_mag <= q_mag_max) {
            qValues.push_back(q);
        }

        int pos = NDIM - 1;
        int count = 0;
        while (count < n_q - 1) {
            if (_q_int[pos] == q_int[pos]) {
                _q_int[pos] = -q_int[pos];
                pos -= 1;
            } else {
                _q_int[pos] += 1;
                for (int i = 0; i < NDIM; i++) {
                    q[i] = _q_int[i]*2.0*M_PI/path.boxPtr->side[i];
                }
                q_mag = sqrt(dot(q,q));
                if (q_mag <= q_mag_max) {
                    qValues.push_back(q);
                }
                count += 1;
                pos = NDIM - 1; //increment the innermost loop
            }
        }
    }

    if (constants()->wavevectorType() == "file_int") {
        std::cout << "wavevectorType = file-int" << std::endl;

        std::ifstream file(input);
        std::string line;
        // Read one line at a time into the variable line:
        while(std::getline(file, line)) {
            std::vector<int> line_data;
            std::stringstream line_stream(line);
        
            int value;
            // Read an integer at a time from the line
            while(line_stream >> value) {
                // Add the integers from a line to a 1D array (vector)
                line_data.push_back(value);
            }
            PIMC_ASSERT(line_data.size()==NDIM);

            for (int i=0; i < NDIM; i++) {
                q[i] = (2.0*M_PI/path.boxPtr->side[i])*line_data[i];
            }
            qValues.push_back(q);
        }
    }
    if (constants()->wavevectorType() == "file_float") {
        std::cout << "wavevectorType = file-float" << std::endl;

        std::ifstream file(input);
        std::string line;
        // Read one line at a time into the variable line:
        while(std::getline(file, line)) {
            std::vector<int> line_data;
            std::stringstream line_stream(line);
        
            int value;
            // Read an integer at a time from the line
            while(line_stream >> value) {
                // Add the integers from a line to a 1D array (vector)
                line_data.push_back(value);
            }
            PIMC_ASSERT(line_data.size()==NDIM);

            for (int i=0; i < NDIM; i++) {
                q[i] = (2.0*M_PI/path.boxPtr->side[i])*line_data[i];
            }
            qValues.push_back(q);
        }
    }
    if (constants()->wavevectorType() == "help") {
        std::cout << "wavevectorType = help" << std::endl;
        std::cout << std::endl;
        std::cout << "The intermediate scattering function behavior is determined by the isf_input and wavevectorType command line arguments." << std::endl;
        std::cout << "Setting wavevectorType to `help` displays this message." << std::endl;
        std::cout << "Other available options are:" << std::endl;
        std::cout << "    int        - set isf_input to an `N*NDIM` space-separated list of integers `i` where the wavevector components are determined by `i*2*pi/L` for the corresponding simulation cell side `L`" << std::endl;
        std::cout << "    float      - set isf_input to an `N*NDIM` space-separated list of floating point numbers `x`, where sequential values modulo NDIM are the corresponding wavevector components" << std::endl;
        std::cout << "    max-int    - set isf_input to an `NDIM` space-separated list of integers `i` where the wavevector components are determined by all allowable wavevectors between `-i*2*pi/L` to `i*2*pi/L` for the corresponding simulation cell side `L`" << std::endl;
        std::cout << "    max-float  - set isf_input to an `NDIM` space-separated list of floating point numbers `x` where wavevector components are dermined for all allowable wavevectors with magnitudes less than the supplied wavevector" << std::endl;
        std::cout << "    file-int   - set isf_input to the path of a file containing any number of lines with `NDIM` space-separated integers `i` where the wavevector components are determined by `i*2*pi/L` for the corresponding simulation cell side `L`" << std::endl;
        std::cout << "    file-float - set isf_input to the path of a file containing any number of lines `NDIM` space-separated floating point numbers `x` where the wavevector components are determined by the supplied wavevector on each line" << std::endl;
        std::cout << std::endl;

        throw "Set argstring_type to: < int | float | max-int | max-float | file-int | file-float >";
    }
    if (constants()->wavevectorType() == "") {
        std::cout << "wavevectorType not set" << std::endl;
        throw "argstring_type not set (set to: < int | float | max-int | max-float | file-int | file-float | help >)";
    }
    delete[] cstr;

    // Write qValues to disk FIXME should be handled by communicator
    /* std::ofstream outFile((format("qValues-ssf-%s.dat") % constants()->id()).str()); */
    /* for (const auto &e : qValues){ */
    /*    outFile << e << "\n"; */
    /* } */
    /* outFile.flush(); */
    /* outFile.close(); */
    
    numq = qValues.size();
    qValues_dVec.resize(numq);
    for (int nq = 0; nq < numq; nq++) {
        qValues_dVec(nq) = qValues[nq];
    }

    /* Initialize the accumulator for the intermediate scattering function*/
    isf.resize(numq*(int(numTimeSlices/2) + 1));
    isf = 0.0;

    // Create multiple gpu streams
    stream_array.resize(MAX_GPU_STREAMS);
    for (int i = 0; i < MAX_GPU_STREAMS; i++) {
        #ifndef USE_CUDA
        HIP_ASSERT(hipStreamCreate(&stream_array(i)));
        #endif
        #ifdef USE_CUDA
        CUDA_ASSERT(cudaStreamCreate(&stream_array(i)));
        #endif
    }

    /* This is a diagonal estimator that gets its own file */
    initialize(numq*(int(numTimeSlices/2) + 1));

    /* the q-values */
    //header = str(format("#%15.6E") % qMag(0));
    //for (int n = 1; n < numq; n++)
    //    header.append(str(format("%16.6E") % qMag(n)));
    //header.append("\n");

    /* The imaginary time values */
    header = str(format("#%15d") % 0);
    for (unsigned int n = 1; n < isf.size(); n++) {
        header.append(str(format("%16d") % n));
    }
    /* utilize imaginary time translational symmetry */
    norm = 0.5;

    bytes_beads = NDIM*(1 + constants()->initialNumParticles())*sizeof(double);
    bytes_isf = isf.size()*sizeof(double);
    bytes_qvecs = NDIM*numq*sizeof(double);
    #ifndef USE_CUDA
        HIP_ASSERT(hipMalloc(&d_isf, bytes_isf)); // Allocate memory for isf on GPU
        HIP_ASSERT(hipMalloc(&d_qvecs, bytes_qvecs)); // Allocate memory for qvecs on GPU
        HIP_ASSERT(hipMemcpy( d_qvecs, qValues_dVec.data(), bytes_qvecs, hipMemcpyHostToDevice )); // Copy qvecs data to gpu
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaMalloc(&d_isf, bytes_isf)); // Allocate memory for isf on GPU
        CUDA_ASSERT(cudaMalloc(&d_qvecs, bytes_qvecs)); // Allocate memory for qvecs on GPU
        CUDA_ASSERT(cudaMemcpy( d_qvecs, qValues_dVec.data(), bytes_qvecs, cudaMemcpyHostToDevice )); // Copy qvecs data to gpu
    #endif
}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
IntermediateScatteringFunctionEstimatorGpu::~IntermediateScatteringFunctionEstimatorGpu() { 
    for (int i = 0; i < MAX_GPU_STREAMS; i++) {
        #ifndef USE_CUDA
            HIP_ASSERT(hipStreamDestroy(stream_array(i)));
        #endif
        #ifdef USE_CUDA
            CUDA_ASSERT(cudaStreamDestroy(stream_array(i)));
        #endif
    }
    isf.free();
    qValues_dVec.free();
    stream_array.free();

    // Release device memory
    #ifndef USE_CUDA
        HIP_ASSERT(hipFree(d_beads));
        HIP_ASSERT(hipFree(d_qvecs));
        HIP_ASSERT(hipFree(d_isf));
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaFree(d_beads));
        CUDA_ASSERT(cudaFree(d_qvecs));
        CUDA_ASSERT(cudaFree(d_isf));
    #endif
}

/*************************************************************************//**
 *  measure the intermediate scattering function for each value of the 
 *  imaginary time separation tau.
 *
 *  We only compute this for N > 1 due to the normalization.
******************************************************************************/
void IntermediateScatteringFunctionEstimatorGpu::accumulate() {
    int numParticles = path.getTrueNumParticles();
    int numTimeSlices = constants()->numTimeSlices();
    int beta_over_two_idx = numTimeSlices/2;
    int number_of_beads = numParticles*numTimeSlices;
    int NNM = number_of_beads*numParticles;
    //int number_of_connections = int(number_of_beads*(number_of_beads + 1)/2);

    double _inorm = 1.0/number_of_beads;

    auto beads_extent = path.get_beads_extent();
    int full_number_of_beads = beads_extent[0]*beads_extent[1];
    int full_numTimeSlices = beads_extent[0];
    int full_numParticles = beads_extent[1];

    //Size, in bytes, of beads array
    size_t bytes_beads_new = NDIM*full_number_of_beads*sizeof(double);

    #ifndef USE_CUDA
        if (bytes_beads_new > bytes_beads) {
            bytes_beads = bytes_beads_new;
            HIP_ASSERT(hipFree(d_beads));
            HIP_ASSERT(hipMalloc(&d_beads, bytes_beads)); // Allocate memory for beads on GPU
        }
        HIP_ASSERT(hipMemcpy( d_beads, path.get_beads_data_pointer(), bytes_beads, hipMemcpyHostToDevice )); // Copy beads data to gpu
        HIP_ASSERT(hipMemset(d_isf, 0, bytes_isf)); // Set initial isf data to zero
    #endif
    #ifdef USE_CUDA
        if (bytes_beads_new > bytes_beads) {
            bytes_beads = bytes_beads_new;
            CUDA_ASSERT(cudaFree(d_beads));
            CUDA_ASSERT(cudaMalloc(&d_beads, bytes_beads)); // Allocate memory for beads on GPU
        }
        CUDA_ASSERT(cudaMemcpy( d_beads, path.get_beads_data_pointer(), bytes_beads, cudaMemcpyHostToDevice )); // Copy beads data to gpu
        CUDA_ASSERT(cudaMemset(d_isf, 0, bytes_isf)); // Set initial isf data to zero
    #endif

    int grid_size = (NNM + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;

    int stream_idx;
    for (int nq = 0; nq < numq; nq++) {
        for (int Mi = 0; Mi < numTimeSlices; Mi++) {
            stream_idx = (nq*numTimeSlices + Mi) % MAX_GPU_STREAMS;
            #ifndef USE_CUDA
            hipLaunchKernelGGL(gpu_isf, dim3(grid_size), dim3(GPU_BLOCK_SIZE), 0, stream_array(stream_idx),
                d_isf, d_qvecs, d_beads, nq, Mi, _inorm, numq, numTimeSlices, numParticles,
                number_of_beads, full_numTimeSlices, full_numParticles, full_number_of_beads, NNM, beta_over_two_idx);
            #endif
            #ifdef USE_CUDA
            cuda_wrapper::gpu_isf_wrapper(dim3(grid_size), dim3(GPU_BLOCK_SIZE), stream_array(stream_idx),
                d_isf, d_qvecs, d_beads, nq, Mi, _inorm, numq, numTimeSlices, numParticles,
                number_of_beads, full_numTimeSlices, full_numParticles, full_number_of_beads, NNM, beta_over_two_idx);
            #endif
        }
    }
    #ifndef USE_CUDA
        HIP_ASSERT(hipDeviceSynchronize());
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaDeviceSynchronize());
    #endif

    //// Copy isf data back to host
    #ifndef USE_CUDA
        HIP_ASSERT(hipMemcpy(isf.data(), d_isf, bytes_isf, hipMemcpyDeviceToHost)); //Only copy up to beta/2 back to host
    #endif
    #ifdef USE_CUDA
        CUDA_ASSERT(cudaMemcpy(isf.data(), d_isf, bytes_isf, cudaMemcpyDeviceToHost)); //Only copy up to beta/2 back to host
    #endif

    estimator += isf;

}
#endif

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
    estimator(2) += 1.0*numParticles/(M_PI*maxR*maxR*path.boxPtr->side[NDIM-1]);
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
void CylinderLinearPotentialEstimator::accumulate1() {

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

        /* We only take slice 0 */
        bead2[0] = 0;

        /* Sum over particles at slice 0 */
        for (bead2[1] = 0; bead2[1] < path.numBeadsAtSlice(bead2[0]); bead2[1]++) {

            r2 = path(bead2);
            if (!include(r2,maxR)) {
                sep = r2 - r1;
                path.boxPtr->putInBC(sep);
                totV += actionPtr->interactionPtr->V(sep);
            } // bead2 is inside  maxR
        } // bead2

        /* Add the constant piece from the external potential */
        totV += actionPtr->externalPtr->V(r1);

        estimator(n) += totV;
    } // n
}

/*************************************************************************//**
 *  We determine what the effective potential along the axis of the pore.
 *
 *  This method creates a spatially resolved histogram of the total potential
 *  felt by partices in the inner core, removing any self-interactions inside
 *  the cut-off radius.
******************************************************************************/
void CylinderLinearPotentialEstimator::accumulate() {

    double totV = 0.0;
    dVec r1,r2;         // The two bead positions

    dVec sep;           // The bead separation
    beadLocator bead1,bead2;  // The bead locators

    for (int slice = 0; slice < path.numTimeSlices; slice++) {
        bead1[0] = slice;

        for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(slice); bead1[1]++) {

            /* Location of bead 1 */
            r1 = path(bead1);

            /* If we are inside the cutoff cylinder, accumulate the potential */
            if (include(r1,maxR)) {

                totV = 0.0;

                /* Sum over particles */
                bead2[0] = slice;

                for (bead2[1] = 0; bead2[1] < path.numBeadsAtSlice(bead2[0]); bead2[1]++) {

                    r2 = path(bead2);

                    /* Make sure r2 is not inside the central core */
                    if (!include(r2,maxR)) {
                        sep = r2 - r1;
                        path.boxPtr->putInBC(sep);
                        totV += actionPtr->interactionPtr->V(sep);
                    } // bead2 is not inside the core
                } //bead2

                /* Add the contribution of the external potential energy */
                totV += actionPtr->externalPtr->V(r1);

                /* determine the z-index of bead 1 */
                int k = int((0.5*Lz + r1[NDIM-1])/dz);
                if (k < NRADSEP) {
                    estimator(k) += totV; // /constants()->numTimeSlices();
                    norm(k) += 1.0;
                }

            } // bead1 is inside the core

        } // bead1

    } // slice
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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CYLINDER STATIC STRUCTURE FACTOR ESTIMATOR CLASS --------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*************************************************************************//**
 *  Constructor.
 * 
 *  Compute the static structure factor for a fixed set of q vector magnitude.
******************************************************************************/
CylinderStaticStructureFactorEstimator::CylinderStaticStructureFactorEstimator(
        const Path &_path, ActionBase *_actionPtr, const MTRand &_random, 
        double _maxR, int _frequency, string _label) :
    EstimatorBase(_path,_actionPtr,_random,_maxR,_frequency,_label) 
{
    /* The maximum q-vector magnitude to consider (hard-coded for now) */
    double qMax = 4.0; // 1/Å

    /* We choose dq from the smallest possible q-vector, set by PBC */
    double dq = 2.0*M_PI/path.boxPtr->side[NDIM-1];
    
    /* Get the desired q-vectors */
    int numq = 0;
    q = getQVectors2(dq,qMax,numq,"line");

    /* Determine how many q-vector magnitudes  we have */
    numq = q.size();

    /* Initialize the accumulator intermediate scattering function*/
    sf.resize(numq);
    sf = 0.0;

    /* This is a diagonal estimator that gets its own file */
    initialize(numq);

    /* The magnitude of q */
    header = str(format("#%15.6E") % 0.0);
    for (int nq = 1; nq < numq; nq++)  {
        double qMag = sqrt(blitz::dot(q[nq][0],q[nq][0]));
        header.append(str(format("%16.6E") % (qMag)));
    }

    /* Utilize imaginary time translational symmetry */
    norm = 1.0/constants()->numTimeSlices();

    /* Normalize by the number of q-vecs (except when there are none) */
    for (int nq = 0; nq < numq; nq++) {
        if (q[nq].size() > 0)
            norm(nq) /= q[nq].size();
    }

}

/*************************************************************************//**
 *  Destructor.
******************************************************************************/
CylinderStaticStructureFactorEstimator::~CylinderStaticStructureFactorEstimator() { 
    sf.free();
}

/*************************************************************************//**
 *  Measure the static structure factor for each value for each q-magnitude
******************************************************************************/
void CylinderStaticStructureFactorEstimator::accumulate() {

    int numParticles = num1DParticles(path,maxR);
    int numTimeSlices = constants()->numTimeSlices();

    beadLocator bead1,bead2;  // The bead locator
    sf = 0.0; // initialize

    /* q-magnitudes */
    for (auto [nq,cq] : enumerate(q)) {

        /* q-vectors with a fixed q-magnitude*/
        for (const auto &cqvec : cq) {
    
            /* Average over all time slices */
            for (bead1[0] = 0; bead1[0] < numTimeSlices; bead1[0]++) {

                /* This is an equal imaginary time estimator */
                bead2[0] = bead1[0];

                for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(bead1[0]); bead1[1]++) {

                    if (include(path(bead1),maxR)) {
                        /* The bead1 = bead2 part */
                        sf(nq) += 1.0;

                        for (bead2[1] = bead1[1]+1; bead2[1] < path.numBeadsAtSlice(bead2[0]); bead2[1]++) {

                            if (include(path(bead2),maxR)) {
                                sf(nq) += 2*cos(dot(path.getSeparation(bead1,bead2),cqvec));
                            } //include bead2

                        } // bead2[1]
                    } //include bead1
                } // bead1[1]
            } //bead1[0]

        } // q-vectors
    } // q-magnitudes

    estimator += sf/numParticles; 
}

/**************************************************************************//**
 *  Sample the estimator.
 * 
 *  Here we overload the cylinder static structure factor estimator, as
 *  we only measure when we have some relevant particle separations.
******************************************************************************/
void CylinderStaticStructureFactorEstimator::sample() {
    numSampled++;

    if ( baseSample() && (num1DParticles(path,maxR)> 0) ) {
        totNumAccumulated++;
        numAccumulated++;
        accumulate();
    }
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
        estimator(slice/2) += sum(actionPtr->potential(slice)) + tailV;
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
        if ( beadIndex[1] < constants()->initialNumParticles()) {
            r = path.getSeparation(beadIndex0,beadIndex);
            estimator(beadIndex[1]-1) += dot(r,r);
        }
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
 *  @see https://mathworld.wolfram.com/SpherePointPicking.html
 *  @see http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
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
        double phi   = acos(2*random.rand() - 1.0);
        rVec[0] = r*sin(theta)*cos(phi);
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
