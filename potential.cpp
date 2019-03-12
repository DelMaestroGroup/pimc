/**
 * @file potential.cpp
 * @author Adrian Del Maestro
 *
 * @brief Implementation of all potential types.
 */

#include "potential.h"
#include "path.h"
#include "lookuptable.h"
#include "communicator.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// POTENTIAL BASE CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
PotentialBase::PotentialBase () : tailV(0.0) {

}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
PotentialBase::~PotentialBase () {
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * The default version creates a list of particle positions in an equally
 * spaced grid.
******************************************************************************/
Array<dVec,1> PotentialBase::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    /* Get the linear size per particle, and the number of particles */
    double initSide = pow((1.0*numParticles/boxPtr->volume),-1.0/(1.0*NDIM));

    /* We determine the number of initial grid boxes there are in 
     * in each dimension and compute their size */
    int totNumGridBoxes = 1;
    iVec numNNGrid;
    dVec sizeNNGrid;

    for (int i = 0; i < NDIM; i++) {
        numNNGrid[i] = static_cast<int>(ceil((boxPtr->side[i] / initSide) - EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = boxPtr->side[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    /* Now, we place the particles at the middle of each box */
    PIMC_ASSERT(totNumGridBoxes>=numParticles);
    dVec pos;
    for (int n = 0; n < totNumGridBoxes; n++) {

        iVec gridIndex;
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid[j];
            gridIndex[i] = (n/scale) % numNNGrid[i];
        }

        for (int i = 0; i < NDIM; i++) 
            pos[i] = (gridIndex[i]+0.5)*sizeNNGrid[i] - 0.5*boxPtr->side[i];

        boxPtr->putInside(pos);

        if (n < numParticles)
            initialPos(n) = pos;
        else 
            break;
    }

    return initialPos;
}

/**************************************************************************//**
 * Ouptut the potential.
 *
 * For use during comparison and debugging, we output the potential out to
 * a supplied separation.
 *
 * @param maxSep the maximum separation
******************************************************************************/
void PotentialBase::output(const double maxSep) {
    dVec sep;
    sep = 0.0;
    for (double d = 0; d < maxSep; d+= (maxSep/1000.0)) {
        sep[0] = d;
        communicate()->file("debug")->stream() 
            << format("%10.4E\t%16.8E\n") % d % V(sep);
    }
}

/**************************************************************************//**
* Return the minimum image difference for 1D separations 
******************************************************************************/
double PotentialBase::deltaSeparation(double sep1, double sep2) const {
    
    double delta = sep2-sep1;
    while (delta >= 0.5*constants()->L())
        delta -= constants()->L();
    while (delta < -0.5*constants()->L()) 
        delta += constants()->L();

    return delta;
}

/**************************************************************************//**
 * Initialize getExcLen method.  
 *
 * This is only used for Gasparini potential, could probably be better.
******************************************************************************/
Array<double,1> PotentialBase::getExcLen() {
    /* The particle configuration */
    Array<double,1> excLens(0);
    return excLens;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// TABULATED POTENTIAL CLASS -------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
TabulatedPotential::TabulatedPotential() { 
    extV = 0.0;
    extdVdr = 0.0;
    extd2Vdr2 = 0.0;
}

/**************************************************************************//**
 * Destructor. 
******************************************************************************/
TabulatedPotential::~TabulatedPotential() {
    lookupV.free();
    lookupdVdr.free();
    lookupd2Vdr2.free();
}

/**************************************************************************//**
 *  Given a discretization factor and the system size, create and fill
 *  the lookup tables for the potential and its derivative.
******************************************************************************/
void TabulatedPotential::initLookupTable(const double _dr, const double maxSep) {

    /* We now calculate the lookup tables for the interaction potential and 
     * its first and second derivatives. */
    dr = _dr;
    tableLength = int(maxSep/dr);
    lookupV.resize(tableLength);
    lookupdVdr.resize(tableLength);
    lookupd2Vdr2.resize(tableLength);
    lookupV = 0.0;
    lookupdVdr = 0.0;
    lookupd2Vdr2 = 0.0;

    double r = 0;

    for (int n = 0; n < tableLength; n++) {
        lookupV(n)    = valueV(r);
        lookupdVdr(n) = valuedVdr(r);
        lookupd2Vdr2(n) = valued2Vdr2(r);
        r += dr;
    }

//  double rc = constants()->rc();
//  for (int n = 0; n < tableLength; n++) {
//      r += dr;
//      if (r <= rc) {
//          lookupV(n) = valueV(r) - valueV(rc) - valuedVdr(rc)*(r-rc);
//          lookupdVdr(n) = valuedVdr(r) - valuedVdr(rc);
//      }
//      else {
//          lookupV(n) = 0.0;
//          lookupdVdr(n) = 0.0;
//      }
//      cout << format("%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n") % r % lookupV(n) % valueV(r) % 
//          lookupdVdr(n) % valuedVdr(r) % (lookupV(n) - valueV(r)) % (lookupdVdr(n) - valuedVdr(r));
//  }

}

/**************************************************************************//**
 *  Use the Newton-Gregory forward difference method to do a 2-point lookup
 *  on the potential table.  
 *
 *  @see M.P. Allen and D.J. Tildesley, "Computer Simulation of Liquids" 
 *  (Oxford Press, London, England) p 144 (2004).
******************************************************************************/
double TabulatedPotential::newtonGregory(const Array<double,1> &VTable, 
        const TinyVector<double,2> &extVal, const double r) {

    double rdr = r/dr;
    int k = int(rdr);

    if (k <= 0) 
        return extVal[0];

    if (k >= tableLength)
        return extVal[1];

    double xi = rdr - 1.0*k;
    double vkm1 = VTable(k-1);
    double vk = VTable(k);
    double vkp1 = VTable(k+1);

    double T1 = vkm1 + (vk - vkm1) * xi;
    double T2 = vk + (vkp1 - vk) * (xi - 1.0);

    return (T1 + 0.5 * (T2 - T1) * xi);
}

/**************************************************************************//**
 *  Use a direct lookup for the potential table.
 *
 *  This is faster thant Newton-Gregory and may give similar results for a fine
 *  enough mesh.
******************************************************************************/
double TabulatedPotential::direct(const Array<double,1> &VTable, 
        const TinyVector<double,2> &extVal, const double r) {

    int k = int(r/dr);
    if (k <= 0) 
        return extVal[0];

    if (k >= tableLength)
        return extVal[1];

    return VTable(k);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// FREE POTENTIAL CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
FreePotential::FreePotential() : PotentialBase() {
} 

/**************************************************************************//**
 * Destructor.
******************************************************************************/
FreePotential::~FreePotential() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SINGLE WELL POTENTIAL CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
SingleWellPotential::SingleWellPotential() : PotentialBase() {
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
SingleWellPotential::~SingleWellPotential() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARMONIC POTENTIAL CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
HarmonicPotential::HarmonicPotential(double _omega) : PotentialBase() {
    omega2 = _omega*_omega;
}

HarmonicPotential::HarmonicPotential() : PotentialBase() {
    omega2 = 1.0;
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
HarmonicPotential::~HarmonicPotential() {
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * We create particles at random locations close to the origin.
******************************************************************************/
Array<dVec,1> HarmonicPotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    for (int n = 0; n < numParticles; n++) {
        for (int i = 0; i < NDIM; i++) 
            initialPos(n)[i] = 0.1*boxPtr->side[i]*(-1.0 + 2.0*random.rand());
    }

    return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARMONIC TUBE POTENTIAL CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
 *
 * Using a supplied tube radius, setup the soft harmonic tube potential.
******************************************************************************/
HarmonicCylinderPotential::HarmonicCylinderPotential(const double radius) : PotentialBase()
{
    /* c is a dimensionless constant */
    c = 1.20272;

    /* We have to determine the frequency of the oscillator from it's length.
     * w = \hbar / (m R^2).  It is measured in THz */
    w = 6.35077 / (radius*radius*constants()->m());
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
HarmonicCylinderPotential::~HarmonicCylinderPotential() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// DELTA POTENTIAL CLASS -----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  Setup the delta function strength and normalization constant.
 *  @param _sigma The width of the Gaussian
 *  @param _g The integrated strength of the Gaussian in 1D
******************************************************************************/
DeltaPotential::DeltaPotential(double _sigma, double _g) : PotentialBase() 
{
    /* Define the parameters of the delta function potential. */
    norm = _g/sqrt(2.0*_sigma*_sigma*M_PI);
    inv2sigma2 = 1.0/(2.0*_sigma*_sigma);
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
DeltaPotential::~DeltaPotential() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LORENTZIAN POTENTIAL CLASS ------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  Setup the delta function strength and normalization constant.
 *  @param _a The width^2 of the Lorentzian
 *  @param _c The integrated strength of the Lorentzian
******************************************************************************/
LorentzianPotential::LorentzianPotential(double _a, double _c) : PotentialBase() 
{
    /* Define the parameters of the Lorentzian delta function potential. */
    a = _a;
    c = _c;
    norm = 2.0 * c * a / M_PI;
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
LorentzianPotential::~LorentzianPotential() {
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SUTHERLAND POTENTIAL CLASS ------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  Setup the parameters of the potential.
 *  @param _g The interaction strength.
******************************************************************************/
SutherlandPotential::SutherlandPotential(double _g) : PotentialBase() 
{
    g = 2.0*constants()->lambda() * _g * (_g - 1.0);
    pioL = M_PI / constants()->L(); 
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
SutherlandPotential::~SutherlandPotential() {
    // empty
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// DIPOLE POTENTIAL CLASS ----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
******************************************************************************/
DipolePotential::DipolePotential() : PotentialBase() 
{
    /* if (NDIM==2) { */
    /*     int N = constants()->initialNumParticles(); */
    /*     tailV = M_PI * N * N / (constants()->L() * constants()->L() * constants()->rc()); */
    /* } */
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
DipolePotential::~DipolePotential() {
    // empty
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// FIXED AZIZ POTENTIAL CLASS ------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  We load the positions of fixed but interacting particles from disk and
 *  create a new lookup table which will be used to speed up the computation.
 *  The interactions with the fixed particles are assumed to be Aziz.
******************************************************************************/
FixedAzizPotential::FixedAzizPotential(const Container *_boxPtr) :
    aziz(_boxPtr->side) {

    char state;             // Fixed or updateable?
    dVec pos;               // The loaded position

    /* Initialize the cutoff^2 */
    rc2 = constants()->rc2();

    /* We start with an array of size 500 */
    fixedParticles.resize(500);

    /* Here we load both the number and location of fixed helium atoms from disk. */
    numFixedParticles = 0;
    int n = 0;
    while (!communicate()->file("fixed")->stream().eof()) {
        if (communicate()->file("fixed")->stream().peek() == '#') {
            communicate()->file("fixed")->stream().ignore(512,'\n');
        }
        else {
            communicate()->file("fixed")->stream() >> state;
            for (int i = 0; i < NDIM; i++) 
                communicate()->file("fixed")->stream() >> pos[i];

            /* If the particle is labelled with an 'F' it is fixed and should
             * be included here */
            if (state == 'F') {
                numFixedParticles++;
                if (numFixedParticles >= int(fixedParticles.size()))
                    fixedParticles.resizeAndPreserve(numFixedParticles);

                /* Put the initial position in the container */
                _boxPtr->putInside(pos);
                fixedParticles(n) = pos;
                n++;
            }
            communicate()->file("fixed")->stream().ignore();
        }
    }

    fixedParticles.resizeAndPreserve(numFixedParticles);

    /* Now that we have the particle positions, create a new lookup table pointer
     * and initialize it */
    lookupPtr = new LookupTable(_boxPtr,1,numFixedParticles);
    lookupPtr->updateGrid(fixedParticles);

    /* Resize and initialize our local grid box arrays */
    fixedBeadsInGrid.resize(lookupPtr->getTotNumGridBoxes(),numFixedParticles);
    numFixedBeadsInGrid.resize(lookupPtr->getTotNumGridBoxes());
    fixedBeadsInGrid = XXX;
    numFixedBeadsInGrid = 0;

    /* Create a local copy of all beads in each grid box plus nearest neighbors.
     * This will drastically speed up the computing of potential energies. */
    for (n = 0; n < lookupPtr->getTotNumGridBoxes(); n++) {
        lookupPtr->updateFullInteractionList(n,0);
        numFixedBeadsInGrid(n) = lookupPtr->fullNumBeads;
        for (int m = 0; m < lookupPtr->fullNumBeads; m++) 
            fixedBeadsInGrid(n,m) = lookupPtr->fullBeadList(m)[1];
    }

}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
FixedAzizPotential::~FixedAzizPotential() {
    delete lookupPtr;
    fixedParticles.free();
    fixedBeadsInGrid.free();
    numFixedBeadsInGrid.free();
}

/**************************************************************************//**
 *  The total potential coming from the interaction of a particle with all 
 *  fixed particles.
******************************************************************************/
double FixedAzizPotential::V(const dVec &pos) {

    double totV = 0.0;

    /* We first find the grid box number that the particle resides in */
    int gridNumber = lookupPtr->gridNumber(pos);

    /* We now loop over all fixed particles in this grid box, only computing
     * interactions when the separation is less than the cutoff */
    dVec sep;
    for (int n = 0; n < numFixedBeadsInGrid(gridNumber); n++) {
        sep = fixedParticles(fixedBeadsInGrid(gridNumber,n)) - pos;
        lookupPtr->boxPtr->putInBC(sep);
        if (dot(sep,sep) < rc2)
            totV += aziz.V(sep);
    }

    return totV;
}

/**************************************************************************//**
 *  The gradient of the total potential coming from the interaction of a 
 *  particle with all fixed particles.
******************************************************************************/
dVec FixedAzizPotential::gradV(const dVec &pos) {

    dVec totGradV;
    totGradV = 0.0;

    /* We first find the grid box number that the particle resides in */
    int gridNumber = lookupPtr->gridNumber(pos);

    /* We now loop over all fixed particles in this grid box, only computing
     * the gradient of interactions when the separation is less than the cutoff */
    dVec sep;
    for (int n = 0; n < numFixedBeadsInGrid(gridNumber); n++) {
        sep = fixedParticles(fixedBeadsInGrid(gridNumber,n)) - pos;
        lookupPtr->boxPtr->putInBC(sep);
        if (dot(sep,sep) < rc2)
            totGradV += aziz.gradV(sep);
    }

    return totGradV;
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * We load an initial configuration from disk, which consists of a number of
 * updateable positions.  These positions are stored as NDIM vectors and
 * proceeded by a letter 'U'.
******************************************************************************/
Array<dVec,1> FixedAzizPotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(1);
    initialPos = 0.0;

    int locNumParticles = 0;    // Number of updateable particles
    char state;                 // Update or Fix
    dVec pos;                   // The current position

    /* We go through all lines in the fixed input file, discarding any comments
     * and assign the initial positions of the particles */
    int n = 0;
    while (!communicate()->file("fixed")->stream().eof()) {
        if (communicate()->file("fixed")->stream().peek() == '#') {
            communicate()->file("fixed")->stream().ignore(512,'\n');
        }
        else {
            communicate()->file("fixed")->stream() >> state;
            for (int i = 0; i < NDIM; i++)
                communicate()->file("fixed")->stream() >> pos[i];

            /* If the particle is labelled with an 'U' it is updatable and should
             * be included */
            if (state == 'U') {
                locNumParticles++;
                initialPos.resizeAndPreserve(locNumParticles);

                /* Put the initial position in the box */
                boxPtr->putInside(pos);

                /* Assign the position to all time slices*/
                initialPos(n) = pos;
                n++;
            }
            communicate()->file("fixed")->stream().ignore();
        }
    }

    /* Reset the file pointer */
    communicate()->file("fixed")->stream().seekg(0, ios::beg);

    /* Return the initial Positions */
    return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// FixedPositionLJPotential Class---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
FixedPositionLJPotential::FixedPositionLJPotential (double _sigma, double _epsilon, 
        const Container *_boxPtr) : PotentialBase() {

    boxPtr = _boxPtr;
    sigma = _sigma;
    epsilon = _epsilon;

    Lz = boxPtr->side[NDIM-1];

    /* Fixed positions of FILENAME */
    dVec pos;               // The loaded position

    /* We start with an array of size 500 */
    fixedParticles.resize(500);

    /* Here we load both the number and location of fixed positions from disk. */
    numFixedParticles = 0;
    int n = 0;
    while (!communicate()->file("fixed")->stream().eof()) {
        if (communicate()->file("fixed")->stream().peek() == '#') {
            communicate()->file("fixed")->stream().ignore(512,'\n');
        }
        else {
            for (int i = 0; i < NDIM; i++) 
                communicate()->file("fixed")->stream() >> pos[i];
            numFixedParticles++;
            if (numFixedParticles >= int(fixedParticles.size()))
                fixedParticles.resizeAndPreserve(numFixedParticles);
            fixedParticles(n) = pos;
            n++;
            communicate()->file("fixed")->stream().ignore();
        }
    }
    fixedParticles.resizeAndPreserve(numFixedParticles);

    /* print out the potential to disk */
    /* int numPoints = 200; */
    /* double dx = _boxPtr->side[0]/numPoints; */
    /* double dy = _boxPtr->side[1]/numPoints; */
    /* pos[2] = -0.5*_boxPtr->side[2] + 2.635; */
    /* for (int i = 0; i < numPoints; i++) { */
    /*     pos[0] = -0.5*_boxPtr->side[0] + i*dx; */
    /*     for (int j = 0; j < numPoints; j++) { */
    /*         pos[1] = -0.5*_boxPtr->side[1] + j*dy; */
    /*         communicate()->file("debug")->stream() << format("%24.16e %24.16e %24.16e\n") % pos[0] % pos[1] % V(pos); */
    /*     } */
    /* } */

    /* exit(-1); */
}


/**************************************************************************//**
 * Destructor.
******************************************************************************/
FixedPositionLJPotential::~FixedPositionLJPotential() {
    fixedParticles.free();
}

/**************************************************************************//**
 *  Return the value of the van der Waals' interaction between a graphene sheet
 *  and a helium adatom at a position, r, above the sheet. 

 *  @param r the position of a helium particle
 *  @return the van der Waals' potential for graphene-helium
******************************************************************************/
double FixedPositionLJPotential::V(const dVec &r) {

    /* Notes: for now I hard-code the potential at 1.5 \AA and a LJ-cutoff of
     * 20 \AA */
    
    if (r[NDIM-1] < (-0.5*Lz + 1.5) )
        return 87292.0;

    else if (r[NDIM-1] > 0.0)
        return 0.0;

    double v = 0.0;
    double sor = 0.0;
    double x = 0.0;
    dVec sep;
    for (int i = 0; i < numFixedParticles; i++) { 
        sep[0] = fixedParticles(i)[0] - r[0];
        sep[1] = fixedParticles(i)[1] - r[1];
        boxPtr->putInBC(sep);
        sep[2] = fixedParticles(i)[2] - r[2];
        x = sqrt(dot(sep,sep));
        if (x < 20.0) {
            sor = sigma/x;
            v += pow(sor,12)-pow(sor,6);
        }
    }
    return 4.0*epsilon*v;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARD CYLINDER POTENTIAL CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *  @param radius The radius of the cylinder
******************************************************************************/
HardCylinderPotential::HardCylinderPotential(const double radius) : 
    PotentialBase(),
    R(radius) {
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
HardCylinderPotential::~HardCylinderPotential() {
}

/**************************************************************************//**
 *  Constructor.
 *
 *  We create a nested Lennard-Jones cylinder, which uses a lookup table to hold the
 *  value of the integrated 6-12 potential for helium atoms interacting
 *  with the walls of a nanpore pre-plated with an adsorbed gas.
 *  @param radius The radius of the cylinder
******************************************************************************/
PlatedLJCylinderPotential::PlatedLJCylinderPotential(const double Ro_, const double Rw, 
        const double sigmaPlated_, const double epsilonPlated_, const double densityPlated_) : 
    PotentialBase(),
    TabulatedPotential()
{
    /* The outer radius of the tube */
    Ro = Ro_;

    /* Plating substrates */
    densityPlated = densityPlated_; //0.207; // atoms / angstrom^3
    epsilonPlated = epsilonPlated_; //36.13613;    // Kelvin
    sigmaPlated   = sigmaPlated_; //3.0225;    // angstroms
    Ri = Ro-Rw; //3.5; 

    // Neon Values
    /* densityPlated = 0.207; // atoms / angstrom^3 FIXME */
    /* epsilonPlated = 20.161242;    // Kelvin */
    /* sigmaPlated   = 2.711;    // angstroms */
    /* Ri = Ro-(1.52*2); //FIXME Fix hardcoded R2 */

    /* Substrate parameters for MCM-41 obtained by fitting 
     * @see https://nano.delmaestro.org/index.php?title=Effective_external_potential_(nanopores) 
     */
    density = 1.000;   // atoms / angstrom^3
    epsilon = 1.59;    // Kelvin
    sigma   = 3.44;    // angstroms

    /* We choose a mesh consisting of 10^6 points, and create the lookup table */
    dR = (1.0E-6)*Ri;
    initLookupTable(dR,Ri);
    
    /* Find the minimun of the potential */
    minV = 1.0E5;
    for (int n = 0; n < tableLength; n++) {
        if (lookupV(n) < minV)
            minV = lookupV(n);
    }

    /* The extremal values for the lookup table */
    extV = valueV(0.0),valueV(Ri);
    extdVdr = valuedVdr(0.0),valuedVdr(Ri);
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
PlatedLJCylinderPotential::~PlatedLJCylinderPotential() {
}

/**************************************************************************//**
 *  Return the actual value of the LJ Cylinder potential, for a distance r 
 *  from the surface of the wall.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double PlatedLJCylinderPotential::V_(const double r, const double R, 
                                   const double sig, const double eps, const double dens) 
{
    double x = (r - (r>=R)*EPS) / R;

    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    double f1 = 1.0 / (1.0 - x2);
    double sigoR3 = pow(sig/R,3.0);
    double sigoR9 = sigoR3*sigoR3*sigoR3;

    double Kx2 = boost::math::ellint_1(x);
    double Ex2 = boost::math::ellint_2(x);

    double v9 = (1.0*pow(f1,9.0)/(240.0)) * (
            (1091.0 + 11156*x2 + 16434*x4 + 4052*x6 + 35*x8)*Ex2 - 
            8.0*(1.0 - x2)*(1.0 + 7*x2)*(97.0 + 134*x2 + 25*x4)*Kx2);
    double v3 = 2.0*pow(f1,3.0) * ((7.0 + x2)*Ex2 - 4.0*(1.0-x2)*Kx2);

    return ((M_PI*eps*sig*sig*sig*dens/3.0)*(sigoR9*v9 - sigoR3*v3));
}

/**************************************************************************//**
 *  Return the r-derivative of the LJ Cylinder potential for separation r.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double PlatedLJCylinderPotential::dVdr_(const double r, const double R, 
                                   const double sig, const double eps, const double dens) 
{
    double x = (r - (r>=R)*EPS) / R;

    /* dV/dr */
    if (x < EPS)
        return (1.28121E8/pow(R,11.0) - 102245.0/pow(R,5.0))*x;
    else {
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x2*x4;
        double x8 = x4*x4;
        double f1 = 1.0 / (1.0 - x2);
        double sigoR3 = pow(sig/R,3.0);
        double sigoR9 = sigoR3*sigoR3*sigoR3;

        double Kx2 = boost::math::ellint_1(x);
        double Ex2 = boost::math::ellint_2(x);

        double dv9dx =(3.0*pow(f1,10.0)/(80.0*x)) *
            ( (1.0 + x2)*(35.0 + 5108*x2 + 22482*x4 + 5108*x6 + 35*x8)*Ex2 -
              (1.0 - x2)*(35.0 + 3428*x2 + 15234*x4 + 12356*x6 +1715*x8)*Kx2 );
        double dv3dx = (6.0*pow(f1,4.0)/x) *
            ( (1.0 + 14*x2 + x4)*Ex2 - (1.0 + 6*x2 - 7*x4)*Kx2 );
        return ((M_PI*eps*sig*sig*sig*dens/(3.0*R))*(sigoR9*dv9dx - sigoR3*dv3dx));
    }
}

/**************************************************************************//**
 *  Return the actual value of the pre-plated cylinder potential, for a 
 *  distance r from the surface of the wall.
******************************************************************************/
double PlatedLJCylinderPotential::valueV(const double r) {

    return V_(r,Ri,sigmaPlated,epsilonPlated,densityPlated) - 
           V_(r,Ro,sigmaPlated,epsilonPlated,densityPlated) +
           V_(r,Ro,sigma,epsilon,density);
}

/**************************************************************************//**
 *  Return the r-derivative of the LJ Cylinder potential for separation r.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double PlatedLJCylinderPotential::valuedVdr(const double r) {
    return dVdr_(r,Ri,sigmaPlated,epsilonPlated,densityPlated) - 
           dVdr_(r,Ro,sigmaPlated,epsilonPlated,densityPlated) +
           dVdr_(r,Ro,sigma,epsilon,density);
}

/**************************************************************************//**
 * N.B. Need to add this if we want to use the virial estimator.
 *
******************************************************************************/
double PlatedLJCylinderPotential::valued2Vdr2(const double r) {
    return 0.0;
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * Return a set of initial positions inside the cylinder.
******************************************************************************/
Array<dVec,1> PlatedLJCylinderPotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    /* We shift the radius inward to account for the excluded volume from the
     * hard wall.  This represents the largest prism that can be put
     * inside a cylinder. */
    dVec lside;
    lside[0] = lside[1] = sqrt(2.0)*(Ri-sigmaPlated);
    lside[2] = boxPtr->side[NDIM-1];

    /* Get the linear size per particle */
    double initSide = pow((1.0*numParticles/product(lside)),-1.0/(1.0*NDIM));

    /* We determine the number of initial grid boxes there are in 
     * in each dimension and compute their size */
    int totNumGridBoxes = 1;
    iVec numNNGrid;
    dVec sizeNNGrid;

    for (int i = 0; i < NDIM; i++) {
        numNNGrid[i] = static_cast<int>(ceil((lside[i] / initSide) - EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = lside[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    /* Now, we place the particles at the middle of each box */
    PIMC_ASSERT(totNumGridBoxes>=numParticles);
    dVec pos;
    for (int n = 0; n < totNumGridBoxes; n++) {

        iVec gridIndex;
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid[j];
            gridIndex[i] = (n/scale) % numNNGrid[i];
        }

        for (int i = 0; i < NDIM; i++) 
            pos[i] = (gridIndex[i]+0.5)*sizeNNGrid[i] - 0.5*lside[i];

        boxPtr->putInside(pos);

        if (n < numParticles)
            initialPos(n) = pos;
        else 
            break;
    }

    return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LJ CYLINDER POTENTIAL CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *
 *  We create a Lennard-Jones cylinder, which uses a lookup table to hold the
 *  value of the integrated 6-12 potential for helium atoms interacting
 *  with a silicon nitride cylinder.
 *  @see C. Chakravarty J. Phys. Chem. B,  101, 1878 (1997).
 *  @param radius The radius of the cylinder
******************************************************************************/
LJCylinderPotential::LJCylinderPotential(const double radius) : 
    PotentialBase(),
    TabulatedPotential()
{

    /* The radius of the tube */
    R = radius;

    /* The density of nitrogen in silicon nitride */
    density = 0.078; // atoms / angstrom^3
//  density = 0.008; // atoms / angstrom^3

    /* We define the values of epsilon and sigma for N and He */ 
//  double epsilonHe = 10.216;  // Kelvin
//  double sigmaHe   = 2.556;   // angstroms
//  double sigmaN    = 3.299;   // angstroms
//  double epsilonN  = 36.2;    // Kelvin
//  epsilon = sqrt(epsilonHe*epsilonN);
//  sigma = 0.5*(sigmaHe + sigmaN);

    /* We operate under the assumption that the silicon can be neglected in the 
     * silicon-nitride, and thus only consider the Nitrogen.  We use a
     * Kiselov type model to extract the actual parameters.  We assume that
     * silicate and silicon-nitride are roughly equivalent. */
    epsilon = 10.22;    // Kelvin
    sigma   = 2.628;    // angstroms

//  epsilon = 32;   // Kelvin
//  sigma   = 3.08; // angstroms

    /* We choose a mesh consisting of 10^6 points, and create the lookup table */
    dR = (1.0E-6)*R;

    initLookupTable(dR,R);
    
    /* Find the minimun of the potential */
    minV = 1.0E5;
    for (int n = 0; n < tableLength; n++) {
        if (lookupV(n) < minV)
            minV = lookupV(n);
    }

    /* The extremal values for the lookup table */
    extV = valueV(0.0),valueV(R);
    extdVdr = valuedVdr(0.0),valuedVdr(R);
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
LJCylinderPotential::~LJCylinderPotential() {
}

/**************************************************************************//**
 *  Return the actual value of the LJ Cylinder potential, for a distance r 
 *  from the surface of the wall.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double LJCylinderPotential::valueV(const double r) {
    double x = r / R;

    if (x >= 1.0)
        x = 1.0 - EPS;
    
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    double f1 = 1.0 / (1.0 - x2);
    double sigoR3 = pow(sigma/R,3.0);
    double sigoR9 = sigoR3*sigoR3*sigoR3;

    double Kx2 = boost::math::ellint_1(x);
    double Ex2 = boost::math::ellint_2(x);

    double v9 = (1.0*pow(f1,9.0)/(240.0)) * (
            (1091.0 + 11156*x2 + 16434*x4 + 4052*x6 + 35*x8)*Ex2 - 
            8.0*(1.0 - x2)*(1.0 + 7*x2)*(97.0 + 134*x2 + 25*x4)*Kx2);
    double v3 = 2.0*pow(f1,3.0) * ((7.0 + x2)*Ex2 - 4.0*(1.0-x2)*Kx2);

    return ((M_PI*epsilon*sigma*sigma*sigma*density/3.0)*(sigoR9*v9 - sigoR3*v3));
}

/**************************************************************************//**
 *  Return the r-derivative of the LJ Cylinder potential for separation r.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double LJCylinderPotential::valuedVdr(const double r) {

    double x = r / R; 

    if (x >= 1.0) 
        x = 1.0 - EPS;

    /* dV/dr */
    if (x < EPS)
        return (1.28121E8/pow(R,11.0) - 102245.0/pow(R,5.0))*x;
    else {
        double x2 = x*x;
        double x4 = x2*x2;
        double x6 = x2*x4;
        double x8 = x4*x4;
        double f1 = 1.0 / (1.0 - x2);
        double sigoR3 = pow(sigma/R,3.0);
        double sigoR9 = sigoR3*sigoR3*sigoR3;

        double Kx2 = boost::math::ellint_1(x);
        double Ex2 = boost::math::ellint_2(x);

        double dv9dx =(3.0*pow(f1,10.0)/(80.0*x)) *
            ( (1.0 + x2)*(35.0 + 5108*x2 + 22482*x4 + 5108*x6 + 35*x8)*Ex2 -
              (1.0 - x2)*(35.0 + 3428*x2 + 15234*x4 + 12356*x6 +1715*x8)*Kx2 );
        double dv3dx = (6.0*pow(f1,4.0)/x) *
            ( (1.0 + 14*x2 + x4)*Ex2 - (1.0 + 6*x2 - 7*x4)*Kx2 );
        return ((M_PI*epsilon*sigma*sigma*sigma*density/(3.0*R))*(sigoR9*dv9dx - sigoR3*dv3dx));
    }
}

/**************************************************************************//**
 *  Return the second r-derivative of the LJ Cylinder potential for separation r.
 *
 * This has been checked with Mathematica --MTG.
******************************************************************************/
double LJCylinderPotential::valued2Vdr2(const double r) {
 
    double x = r / R; 

    if (x >= 1.0) 
        x = 1.0 - EPS;

    /* d2V/dr2 */
    /*if (x < EPS){
    // related to hard core limit, this will likely need to be implemented. 
    return (1.28121E8/pow(R,11.0) - 102245.0/pow(R,5.0))*x;
    }
    else {*/
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    double x10 = x8*x2;
    double f1 = 1.0 / (1.0 - x2);
    double sigoR3 = pow(sigma/R,3.0);
    double sigoR9 = sigoR3*sigoR3*sigoR3;

    double Kx2 = boost::math::ellint_1(x);
    double Ex2 = boost::math::ellint_2(x);

    double d2v9dx2 = (2.0/240)*pow(f1,11)*((1925.0*x10 + 319451.0*x8 + 2079074.0*x6 + 2711942.0*x4
                + 764873.0*x2 + 20975.0)*Ex2
            - (729400.0*x10 + 430024.0*x8 + 344752.0*x6 + 767248.0*x4 + 386200.0*x2 + 12712.0)*Kx2); // checked --MTG
    double d2v3dx2 = 8.0*pow(f1,5)*((11.0 + 80.0*x2  + 5.0*x4)*Ex2 + 4.0*(5.0*x4 - 4.0*x2 - 1.0)*Kx2); // checked --MTG

    return ((M_PI*epsilon*sigma*sigma*sigma*density/(3.0*R*R))*(sigoR9*d2v9dx2 - sigoR3*d2v3dx2));
    //}
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * Return a set of initial positions inside the cylinder.
******************************************************************************/
Array<dVec,1> LJCylinderPotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    /* We shift the radius inward to account for the excluded volume from the
     * hard wall.  This represents the largest prism that can be put
     * inside a cylinder. */
    dVec lside;
    lside[0] = lside[1] = sqrt(2.0)*(R-sigma);
    lside[2] = boxPtr->side[NDIM-1];

    /* Get the linear size per particle */
    double initSide = pow((1.0*numParticles/product(lside)),-1.0/(1.0*NDIM));

    /* We determine the number of initial grid boxes there are in 
     * in each dimension and compute their size */
    int totNumGridBoxes = 1;
    iVec numNNGrid;
    dVec sizeNNGrid;

    for (int i = 0; i < NDIM; i++) {
        numNNGrid[i] = static_cast<int>(ceil((lside[i] / initSide) - EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = lside[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    /* Now, we place the particles at the middle of each box */
    PIMC_ASSERT(totNumGridBoxes>=numParticles);
    dVec pos;
    for (int n = 0; n < totNumGridBoxes; n++) {

        iVec gridIndex;
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid[j];
            gridIndex[i] = (n/scale) % numNNGrid[i];
        }

        for (int i = 0; i < NDIM; i++) 
            pos[i] = (gridIndex[i]+0.5)*sizeNNGrid[i] - 0.5*lside[i];

        boxPtr->putInside(pos);

        if (n < numParticles)
            initialPos(n) = pos;
        else 
            break;
    }

    return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LJ HOUR GLASS POTENTIAL CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *
 *  We create a Lennard-Jones hourglass, which for now, directly evaluates
 *  the potential based on breaking up the pore into finite size pieces, 
 *  each with their own radius.  
 *  @see C. Chakravarty J. Phys. Chem. B,  101, 1878 (1997).
 *
 *  @param radius The radius of the cylinder
 *  @param deltaRadius R(z=0) = radius - deltaRadius
 *  @param deltaWidth The width of the hourglass constriction
******************************************************************************/
LJHourGlassPotential::LJHourGlassPotential(const double radius, 
        const double deltaRadius, const double deltaWidth) : PotentialBase()
{
    /* The radius of the tube */
    R = radius;

    /* The modification in radius: R(z=0) = R - dR */
    dR = deltaRadius;

    /* the length of the pore */
    L = constants()->L();

    /* The inverse length scale over which the radius changes */
    invd = 1.0/deltaWidth;
    R0 = 1.0/tanh(0.5*L*invd);

    /* The density of nitrogen in silicon nitride */
    density = 0.078; // atoms / angstrom^3

    /* These are the values we have historically used for amorphous silicon
     * nitride. */ 
    epsilon = 10.22;    // Kelvin
    sigma   = 2.628;    // angstroms
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
LJHourGlassPotential::~LJHourGlassPotential() {
}

/**************************************************************************//**
 *  Return the actual value of the LJ Cylinder potential, for a distance r 
 *  from the surface of the wall.
 *
 *  Checked and working with Mathematica (hourglass_potential_test.nb) on 
 *  2014-09-04.
 *
 *  @param r the position of a particle
 *  @return the confinement potential
******************************************************************************/
double LJHourGlassPotential::V(const dVec &r) {

    double x = 0.0;
    for (int i = 0; i < NDIM-1; i++)
        x += r[i]*r[i];

    double Rz = Rtanh(r[NDIM-1]);

    x = sqrt(x)/Rz;

    if (x >= 1.0)
        x = 1.0 - EPS;
    
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    double f1 = 1.0 / (1.0 - x2);
    double sigoR3 = pow(sigma/Rz,3.0);
    double sigoR9 = sigoR3*sigoR3*sigoR3;

    double Kx2 = boost::math::ellint_1(x);
    double Ex2 = boost::math::ellint_2(x);

    double v9 = (1.0*pow(f1,9.0)/(240.0)) * (
            (1091.0 + 11156*x2 + 16434*x4 + 4052*x6 + 35*x8)*Ex2 - 
            8.0*(1.0 - x2)*(1.0 + 7*x2)*(97.0 + 134*x2 + 25*x4)*Kx2);
    double v3 = 2.0*pow(f1,3.0) * ((7.0 + x2)*Ex2 - 4.0*(1.0-x2)*Kx2);

    return ((M_PI*epsilon*sigma*sigma*sigma*density/3.0)*(sigoR9*v9 - sigoR3*v3));
}

/**************************************************************************//**
 * Return a set of initial positions inside the hourglass.
 *
 * @param boxPtr A pointer to the simulation cell
 * @param random The random number generator
 * @param numParticles The initial number of particles
 * @return An array of classical particle positions
******************************************************************************/
Array<dVec,1> LJHourGlassPotential::initialConfig1(const Container *boxPtr, 
        MTRand &random, const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    /* We shift the radius inward to account for the excluded volume from the
     * hard wall.  This represents the largest prism that can be put
     * inside a cylinder. */
    dVec lside;
    lside[0] = lside[1] = sqrt(2.0)*(R-dR-sigma);
    lside[2] = boxPtr->side[NDIM-1];

    /* Make sure our radius isn't too small */
    PIMC_ASSERT(lside[0]>0.0);

    /* Get the linear size per particle */
    double initSide = pow((1.0*numParticles/product(lside)),-1.0/(1.0*NDIM));

    /* We determine the number of initial grid boxes there are in 
     * in each dimension and compute their size */
    int totNumGridBoxes = 1;
    iVec numNNGrid;
    dVec sizeNNGrid;

    for (int i = 0; i < NDIM; i++) {
        numNNGrid[i] = static_cast<int>(ceil((lside[i] / initSide) - EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = lside[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    /* Now, we place the particles at the middle of each box */
    PIMC_ASSERT(totNumGridBoxes>=numParticles);
    dVec pos;
    for (int n = 0; n < totNumGridBoxes; n++) {

        iVec gridIndex;
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid[j];
            gridIndex[i] = (n/scale) % numNNGrid[i];
        }

        for (int i = 0; i < NDIM; i++) 
            pos[i] = (gridIndex[i]+0.5)*sizeNNGrid[i] - 0.5*lside[i];

        boxPtr->putInside(pos);

        if (n < numParticles)
            initialPos(n) = pos;
        else 
            break;
    }

    return initialPos;
}

/**************************************************************************//**
 * Return a set of initial positions inside the hourglass.
 *
 * @param boxPtr A pointer to the simulation cell
 * @param random The random number generator
 * @param numParticles The initial number of particles
 * @return An array of classical particle positions
******************************************************************************/
Array<dVec,1> LJHourGlassPotential::initialConfig(const Container *boxPtr, 
        MTRand &random, const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;

    dVec pos;
    pos = 0.0;

    /* We randomly place the particles inside the cylinder taking acount of the 
     * pinched radius. 
     * @see http://mathworld.wolfram.com/DiskPointPicking.html
     */
	for (int n = 0; n < numParticles; n++) {

            /* Uniform position along the pore */
            pos[NDIM-1] = L*(-0.5 + random.rand());

            /* Uniform position in a disk of z-dependent radius*/
            double theta = 2.0*M_PI*random.rand();
            double r = (Rtanh(pos[NDIM-1]) - sigma)*sqrt(random.rand());

            pos[0] = r*cos(theta);
            pos[1] = r*sin(theta);

            boxPtr->putInside(pos);

            initialPos(n) = pos;
    }

	return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// AZIZ POTENTIAL CLASS ------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 * 
 *  Create the Aziz interaction potential.  We use the standard 1979 values.
 *  @see R.A. Aziz et al. J. Chem. Phys. 70, 4330 (1979).
******************************************************************************/
AzizPotential::AzizPotential(const dVec &side) : PotentialBase(), TabulatedPotential()
{
    /* Define all variables for the Aziz potential */
    /* R.A. Aziz et al. J. Chem. Phys. 70, 4330 (1979) */
    rm      = 2.9673;   // A
    A       = 0.5449E6; 
    epsilon = 10.8;     // K
    alpha   = 13.353; 
    D       = 1.2413; 
    C6      = 1.3732;
    C8      = 0.42538;
    C10     = 0.1781;

    /* The extremal values are all zero here */
    extV = 0.0;
    extdVdr = 0.0;
    extd2Vdr2 = 0.0;

    /* We take the maximum possible separation */
    double L = max(side);

    /* Create the potential lookup tables */
    // initLookupTable(0.00005*rm,L);
    initLookupTable((1.0E-6)*rm,L);

    /* Now we compute the tail correction */
    double rmoL = rm / L;
    double rm3 = rm*rm*rm;
    double t1 = A*exp(-alpha*L/(2.0*rm))*rm*(8.0*rm*rm + 4.0*L*rm * alpha + L*L*alpha*alpha)
        / (4.0*alpha*alpha*alpha);
    double t2 = 8.0*C6*pow(rmoL,3.0)/3.0;
    double t3 = 32.0*C8*pow(rmoL,5.0)/5.0;
    double t4 = 128.0*C10*pow(rmoL,7.0)/7.0;
    
    tailV = 2.0*M_PI*epsilon*(t1 - rm3*(t2+t3+t4));
}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
AzizPotential::~AzizPotential() {
}

/**************************************************************************//**
 *  Return the actual value of the Aziz potential, used to construct a lookup
 *  table.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double AzizPotential::valueV(const double r) {
    double x = r / rm;

    double Urep = A * exp(-alpha*x);

    /* No self interactions */
    if (x < EPS) 
        return 0.0;
    /* Hard core limit */
    else if (x < 0.01)
        return (epsilon * Urep);
    else {
        double ix2 = 1.0 / (x * x);
        double ix6 = ix2 * ix2 * ix2;
        double ix8 = ix6 * ix2;
        double ix10 = ix8 * ix2;
        double Uatt = -( C6*ix6 + C8*ix8 + C10*ix10 ) * F(x);
        return ( epsilon * (Urep + Uatt) );
    }
}

/**************************************************************************//**
 *  Return the r-derivative of the Aziz potential for separation r.
 *
 *  Checked and working with Mathematica.
******************************************************************************/
double AzizPotential::valuedVdr(const double r) {
    double x = r / rm;

    double T1 = -A * alpha * exp(-alpha*x);
    
    /* dV/dR */
    /* No self interactions */
    if (x < EPS) 
        return 0.0;
    /* Hard core limit */
    else if (x < 0.01)
        return ( ( epsilon / rm ) * T1 );
    else {
        /* The various inverse powers of x */
        double ix = 1.0 / x;
        double ix2 = ix*ix;
        double ix6 = ix2 * ix2 * ix2;
        double ix7 = ix6 * ix;
        double ix8 = ix6 * ix2;
        double ix9 = ix8 * ix;
        double ix10 = ix8 * ix2;
        double ix11 = ix10 * ix;
        double T2 = ( 6.0*C6*ix7 + 8.0*C8*ix9 + 10.0*C10*ix11 ) * F(x);
        double T3 = -( C6*ix6 + C8*ix8 + C10*ix10 ) * dF(x);
        return ( ( epsilon / rm ) * (T1 + T2 + T3) );
    }
}

/**************************************************************************//**
 *  Return the second r-derivative of the Aziz potential for separation r.
 *
 *  Double checked and working with Mathematica. -MTG
******************************************************************************/
double AzizPotential::valued2Vdr2(const double r) {
    double x = r / rm;

    double T1 = A * alpha * alpha * exp(-alpha*x);
    
    /* d^2V/dR^2 */
    /* No self interactions */
    if (x < EPS) 
        return 0.0;
    /* Hard core limit */
    else if (x < 0.01)
        return ( ( epsilon / rm ) * T1 );
    else {
        /* The various inverse powers of x */
        double ix = 1.0 / x;
        double ix2 = ix*ix;
        double ix6 = ix2 * ix2 * ix2;
        double ix7 = ix6 * ix;
        double ix8 = ix6 * ix2;
        double ix9 = ix8 * ix;
        double ix10 = ix8 * ix2;
        double ix11 = ix10 * ix;
        double ix12 = ix11 * ix;
        double T2 = - ( 42.0*C6*ix8 + 72.0*C8*ix10 + 110.0*C10*ix12 ) * F(x);
        double T3 = 2.0*( 6.0*C6*ix7 + 8.0*C8*ix9 + 10.0*C10*ix11 ) * dF(x);
        double T4 = - ( C6*ix6 + C8*ix8 + C10*ix10 ) * d2F(x);
        return ( ( epsilon / (rm*rm) ) * (T1 + T2 + T3 + T4) );
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// EXCLUDED VOLUME CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/**************************************************************************//**
 * Constructor.
******************************************************************************/
Gasparini_1_Potential::Gasparini_1_Potential(double _az, double _ay, const Container *_boxPtr) :
    PotentialBase(),
    excZ(0.5*_boxPtr->side[2]-_az),
    excY(0.5*_boxPtr->side[1]-_ay),
    V0(1.0E6)
{
    // Empty Constructor
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
Gasparini_1_Potential::~Gasparini_1_Potential() {
    // Empty Destructor
}

/**************************************************************************//**
 * Return initial particle configuration.
 *
 * Return initial positions to exclude volume in the simulation cell.
******************************************************************************/
Array<dVec,1> Gasparini_1_Potential::initialConfig(const Container *boxPtr,
        MTRand &random, const int numParticles) {

    /* label the lengths of the sides of the simulation cell */
    dVec lside;
    lside[0] = boxPtr->side[0];
    lside[1] = boxPtr->side[1];
    lside[2] = boxPtr->side[NDIM-1];

    /* calculate actual volume */
    double volTot = product(lside);

    /* calculate density */
    double density = (1.0*numParticles/volTot);

    /* calculate excluded volume */
    double volExc = lside[0]*(2.0*excY)*(2.0*excZ);

    /* calculate actual number of particles */
    int correctNum = int(numParticles-density*volExc);

    /* The particle configuration */
    Array<dVec,1> initialPos(correctNum);

    /* get linear size per particle  */
    double initSide = pow((1.0*correctNum/(volTot-volExc)),-1.0/(1.0*NDIM));

    /* For accessible volume, determine the number of 
     * initial grid boxes there are in each dimension and compute
     * their size. */
    int totNumGridBoxes1 = 1;
    int totNumGridBoxes2 = 1;
    iVec numNNGrid1;
    iVec numNNGrid2;
    dVec sizeNNGrid1;
    dVec sizeNNGrid2;

    /* divide space into two regions, insert particles appropriately */
    double V1 = (lside[1]-2.0*excY)*(2.0*excZ)*lside[0];
    double V2 = (lside[2]-2.0*excZ)*lside[1]*lside[0];

    double fracV1 = V1/(V1+V2);

    int numIn1 = int(correctNum*fracV1);
    int numIn2 = (correctNum-numIn1);
    
    /* grid space in volume 1 */
    /* x */
    numNNGrid1[0] = static_cast<int>(ceil((1.0*lside[0]/initSide)-EPS));
    if (numNNGrid1[0] < 1)
        numNNGrid1[0] = 1;
    sizeNNGrid1[0] = 1.0*lside[0]/(1.0*numNNGrid1[0]);
    totNumGridBoxes1 *= numNNGrid1[0];
    /* y */
    numNNGrid1[1] = static_cast<int>(ceil(((lside[1]-2.0*excY)/initSide)-EPS));
    if (numNNGrid1[1] < 1)
        numNNGrid1[1] = 1;
    sizeNNGrid1[1] = (lside[1]-2.0*excY)/(1.0*numNNGrid1[1]);
    totNumGridBoxes1 *= numNNGrid1[1];
    /* z */
    numNNGrid1[2] = static_cast<int>(ceil(((2.0*excZ)/initSide)-EPS));
    if (numNNGrid1[2] < 1)
        numNNGrid1[2] = 1;
    sizeNNGrid1[2] = (2.0*excZ)/(1.0*numNNGrid1[2]);
    totNumGridBoxes1 *= numNNGrid1[2];
    
    /* grid space in volume 2 */
    /* x */
    numNNGrid2[0] = static_cast<int>(ceil((1.0*lside[0]/initSide)-EPS));
    if (numNNGrid2[0] < 1)
        numNNGrid2[0] = 1;
    sizeNNGrid2[0] = 1.0*lside[0]/(1.0*numNNGrid2[0]);
    totNumGridBoxes2 *= numNNGrid2[0];
    /* y */
    numNNGrid2[1] = static_cast<int>(ceil((1.0*lside[1]/initSide)-EPS));
    if (numNNGrid2[1] < 1)
        numNNGrid2[1] = 1;
    sizeNNGrid2[1] = 1.0*lside[1]/(1.0*numNNGrid2[1]);
    totNumGridBoxes2 *= numNNGrid2[1];
    /* z */
    numNNGrid2[2] = static_cast<int>(ceil(((lside[2]-2.0*excZ)/initSide)-EPS));
    if (numNNGrid2[2] < 1)
        numNNGrid2[2] = 1;
    sizeNNGrid2[2] = (lside[2]-2.0*excZ)/(1.0*numNNGrid2[2]);
    totNumGridBoxes2 *= numNNGrid2[2];
    
    /* Place particles in the middle of the boxes -- volume 1 */
    PIMC_ASSERT(totNumGridBoxes1>=numIn1);
    dVec pos1;

    for (int n = 0; n < totNumGridBoxes1; n++) {
        iVec gridIndex1;
        /* update grid index */
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid1[j];
            gridIndex1[i] = (n/scale) % numNNGrid1[i];
        }
        /* place particle in position vector, skipping over excluded volume */
        pos1[0] = (gridIndex1[0]+0.5)*sizeNNGrid1[0] + +0.5*lside[0] + 2.0*EPS;
        pos1[1] = (gridIndex1[1]+0.5)*sizeNNGrid1[1] + excY + 2.0*EPS;
        pos1[2] = (gridIndex1[2]+0.5)*sizeNNGrid1[2] - excZ + 2.0*EPS;

        if ((pos1[1]<-excY || pos1[1]>excY) || (pos1[2]<-excZ || pos1[2]>excZ))
            boxPtr->putInside(pos1);

        if (n < numIn1){
            initialPos(n) = pos1;
        }
        else 
            break;
    }

    /* Place particles in the middle of the boxes -- volume 2 */
    PIMC_ASSERT(totNumGridBoxes2>=numIn2);
    dVec pos2;

    for (int n = 0; n < totNumGridBoxes2; n++) {
        iVec gridIndex2;
        /* update grid index */
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid2[j];
            gridIndex2[i] = (n/scale) % numNNGrid2[i];
        }
        /* place particles in position vectors */
        pos2[0] = (gridIndex2[0]+0.5)*sizeNNGrid2[0] + 0.5*lside[0] + 2.0*EPS;
        pos2[1] = (gridIndex2[1]+0.5)*sizeNNGrid2[1] + 0.5*lside[1] + 2.0*EPS;
        pos2[2] = (gridIndex2[2]+0.5)*sizeNNGrid2[2] + excZ + 2.0*EPS;

        if ((pos2[1]<-excY || pos2[1]>excY) || (pos2[2]<-excZ || pos2[2]>excZ))
            boxPtr->putInside(pos2);
        
        if (n < numIn2){
            initialPos(n+numIn1) = pos2;
        }
        else 
            break;
    }
    /* do we want to output the initial config to disk? */
    bool outToDisk = 1;
    ofstream OF;
    if (outToDisk){
        OF.open("./OUTPUT/initialConfig.dat");
        OF<<"# Cartesian Coordinates of initial Positions (X-Y-Z)"<<endl;
        OF<<"# "<<lside[0]<<"\t"<<lside[1]<<"\t"<<lside[2]<<endl;
        OF<<"# "<< excY <<"\t"<< excZ <<endl;
        for (int i=0; i< int(initialPos.size()); i++)
            OF<<initialPos(i)(0)<< "\t"<<initialPos(i)(1)<<"\t"<<initialPos(i)(2)<<endl;
        OF.close();
    }
    
    return initialPos;
}
/*************************************************************************//**
 *  Returns the exclusion lengths ay and az
******************************************************************************/
Array<double,1> Gasparini_1_Potential::getExcLen(){

    Array<double, 1> excLens(2);
    excLens(0) = excY; // was ay
    excLens(1) = excZ; // was az
    
    return (excLens);
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARD SPHERE POTENTIAL CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *  @param _a The radius of the hard sphere (also the scattering length)
******************************************************************************/
HardSpherePotential::HardSpherePotential(double _a) : 
    PotentialBase(),
    a(_a) {

}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
HardSpherePotential::~HardSpherePotential() {
// empty deconstructor
}

/**************************************************************************//**
 *  The effective potential.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-12.
 *
 *  @param sep1 The first separation
 *  @param sep2 The second separation
 *  @return the two-body effective pair potential
******************************************************************************/
double HardSpherePotential::V(const dVec &sep1, const dVec &sep2) 
{

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if ((r1 <= a ) || (r2 <= a)) 
        return LBIG;

    double cosTheta = dot(sep1,sep2)/(r1*r2);

    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*constants()->lambda()*constants()->tau());

    double t2 = (a*(r1+r2) - a*a)/(r1*r2);
    double t3 = 1.0 - t2*exp(t1);

    return -log(t3);
}

/**************************************************************************//**
 *  The derivative of the effective potential with respect to lambda.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-12.
 *
 *  @param sep1 the first separation
 *  @param sep2 the second separation
 *  @return the derivative of the effective potential with respect to lambda
******************************************************************************/
double HardSpherePotential::dVdlambda(const dVec &sep1, const dVec &sep2) 
{

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    double cosTheta = dot(sep1,sep2)/(r1*r2);

    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*constants()->lambda()*constants()->tau());

    double t2 = (a*(r1+r2) - a*a)/(r1*r2);
    double t3 = 1.0 - t2*exp(t1);

    double t4 = (1.0/t3)*(t2*exp(t1))*(t1/constants()->lambda());

    return -t4;
}

/**************************************************************************//**
 *  The derivative of the effective potential with respect to tau.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-12.
 *
 *  @param sep1 the first separation
 *  @param sep2 the second separation
 *  @param lambda \lambda = \hbar^2/2m
 *  @param tau the imaginary timestep tau
 *  @return the derivative of the effective potential with respect to tau
******************************************************************************/
double HardSpherePotential::dVdtau(const dVec &sep1, const dVec &sep2) 
{

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    double cosTheta = dot(sep1,sep2)/(r1*r2);

    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*constants()->lambda()*constants()->tau());

    double t2 = (a*(r1+r2) - a*a)/(r1*r2);
    double t3 = 1.0 - t2*exp(t1);

    double t4 = (1.0/t3)*(t2*exp(t1))*(t1/constants()->tau());

    return -t4;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// 1D Delta Potential Class -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
*  Constructor.
*  @param g The strength of delta interaction
******************************************************************************/
Delta1DPotential::Delta1DPotential(double _g) :
    PotentialBase(),
    g(_g) 
{
    
    erfCO = 7.0;
    l0 = 2.0*sqrt(constants()->lambda()*constants()->tau());
    li = 4.0*constants()->lambda()/g;
    xi = l0/li;
    
    xiSqOver2 = (0.5)*xi*xi;
    xiSqrtPIOver2 = sqrt(M_PI/2.0)*xi;
}

/**************************************************************************//**
*  Destructor.
******************************************************************************/
Delta1DPotential::~Delta1DPotential() {
    // empty deconstructor
}


/**************************************************************************//**
 *  Interaction weight of pair-propagator.
 *  Uses asymptotic form when xi+yt> erfCO
 *  @param yt = |\tilde{x}| + |\tilde{x}'|
 *  @param dxt = \tilde{x}' - \tilde{x}
 ******************************************************************************/
double Delta1DPotential::Wint(double yt, double dxt) {
    
    double W,erfVal,expVal;
    
    if(xi + yt < erfCO) 
    {
        erfVal = erfc( (xi+yt)/sqrt(2.0) );
        expVal = exp( xiSqOver2 + (0.5)*dxt*dxt + xi*yt );
        W = 1.0 - xiSqrtPIOver2*expVal*erfVal;
    } 
    else 
    {
        expVal = exp((-0.5)*(yt*yt-dxt*dxt));
        W = 1.0 - (xi/(xi+yt))*expVal;
    }
    
    return W;
}

/**************************************************************************//**
*  xi derivative of Wint.
*  Uses asymptotic form when xi+yt> erfCO
*  @param yt = |\tilde{x}| + |\tilde{x}'|
*  @param dxt = \tilde{x}' - \tilde{x}
******************************************************************************/
double Delta1DPotential::dWdxi(double yt, double dxt) {
    
    double dW,erfVal,expVal,expVal2;
    
    expVal = exp( (0.5)*(dxt*dxt-yt*yt) );
    
    if(xi + yt < erfCO)
    {
        erfVal = erfc( (xi+yt)/sqrt(2.0) );
        expVal2 = exp( (0.5)*(xi+yt)*(xi+yt) );
        dW = (-1.0)*xi*expVal*( sqrt(0.5*M_PI)*(xi +yt + (1.0/xi) )*expVal2*erfVal - 1.0 );
    }
    else
         dW = (-1.0)*xi*expVal*( (xi + yt + (1.0/xi) )/( xi + yt ) - 1.0 );
    
    return dW;
}

/**************************************************************************//**
*  \tilde{y} derivative of Wint.
*  Uses asymptotic form when xi+yt> erfCO
*  @param yt = |\tilde{x}| + |\tilde{x}'|
*  @param dxt = \tilde{x}' - \tilde{x}
******************************************************************************/
double Delta1DPotential::dWdyt(double yt, double dxt) {
    
    double dW,erfVal,expVal,expVal2;
    
    expVal = exp( (0.5)*(dxt*dxt-yt*yt) );
    
    if(xi + yt < erfCO)
    {
        erfVal = erfc( (xi+yt)/sqrt(2.0) );
        expVal2 = exp( (0.5)*(xi+yt)*(xi+yt) );
        dW = xi*expVal*( 1 - sqrt(0.5*M_PI)*xi*expVal2*erfVal );
    }
    else
        dW = xi*expVal*( 1 - xi/(xi+yt) );
    
    return dW;
}

/**************************************************************************//**
*  \tidle{dx} derivative of Wint.
*  Uses asymptotic form when xi+yt> erfCO
*  @param yt = |\tilde{x}| + |\tilde{x}'|
*  @param dxt = \tilde{x}' - \tilde{x}
******************************************************************************/
double Delta1DPotential::dWddxt(double yt, double dxt) {
    
    double dW,erfVal,expVal;
    
    if (xi + yt < erfCO)
    {
        erfVal = erfc( (xi+yt)/sqrt(2.0) );
        expVal = exp( (0.5)*( dxt*dxt + xi*( xi + 2.0*yt)) );
        dW = (-1.0)*sqrt(0.5*M_PI)*xi*dxt*expVal*erfVal;
    } 
    else
    {
        expVal = exp( (0.5)*(dxt*dxt-yt*yt) );
        dW = (-1.0)*dxt*expVal*( xi/( xi + yt ) );
    }
    
    return dW;
}

/**************************************************************************//**
 *  The effective potential.
 *
 *  Computes the non-local two-body effective pair potential.
 **
 *  @param sep1 The first separation
 *  @param sep2 The second separation
 *  @return the two-body effective pair potential
 ******************************************************************************/
double Delta1DPotential::V(const dVec &sep1, const dVec &sep2)
{
    
    double dxt = deltaSeparation(sep1[0], sep2[0])/l0;
    double yt = (abs(sep1[0])+abs(sep2[0]))/l0;
    
    double W = Wint(yt,dxt);
    
    return (-1.0)*log(W);
}

/**************************************************************************//**
*  The derivative of the effective potential with respect to lambda.
*
*  Computes the non-local two-body effective pair potential.
*
*  @param sep1 the first separation
*  @param sep2 the second separation
*  @return the derivative of the effective potential with respect to lambda
******************************************************************************/
double Delta1DPotential::dVdlambda(const dVec &sep1, const dVec &sep2)
{
    
    double dxt = deltaSeparation(sep1[0], sep2[0])/l0;
    double yt = (abs(sep1[0])+abs(sep2[0]))/l0;
    
    double W = Wint(yt,dxt);
    double dWdy = dWdyt(yt,dxt);
    double dWddx = dWddxt(yt,dxt);
    double dWdx = dWdxi(yt,dxt);

    double dWdl = ((-1.0)/(2.0*constants()->lambda()))*(yt*dWdy + dxt*dWddx+ xi*dWdx );
    
    return ((-1.0)/W)*dWdl;
}

/**************************************************************************//**
 *  The derivative of the effective potential with respect to tau.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  @param sep1 the first separation
 *  @param sep2 the second separation
 *  @return the derivative of the effective potential with respect to lambda
 ******************************************************************************/
double Delta1DPotential::dVdtau(const dVec &sep1, const dVec &sep2)
{
    
    double dxt = deltaSeparation(sep1[0], sep2[0])/l0;
    double yt = (abs(sep1[0])+abs(sep2[0]))/l0;
    
    double W = Wint(yt,dxt);
    double dWdt = ((1.0)/(2.0*constants()->tau()))  
        * ( (-1.0)*yt*dWdyt(yt,dxt) - dxt*dWddxt(yt,dxt) + xi*dWdxi(yt,dxt) );

    return ((-1.0)/W)*dWdt;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARD ROD POTENTIAL CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *  @param _a The radius of the hard rod (also the scattering length)
******************************************************************************/
HardRodPotential::HardRodPotential(double _a) : 
    PotentialBase(),
    a(_a) {

}

/**************************************************************************//**
 *  Destructor.
******************************************************************************/
HardRodPotential::~HardRodPotential() {
// empty deconstructor
}

/**************************************************************************//**
 *  The effective potential.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-17.
 *
 *  @param sep1 The first separation
 *  @param sep2 The second separation
 *  @return the two-body effective pair potential
******************************************************************************/
double HardRodPotential::V(const dVec &sep1, const dVec &sep2) 
{
    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    /* We need to enforce the distinguishable particle constraint at short
     * imaginary times */
    /* if ( (sep1[0]*sep2[0] < 0.0) || (r1 <= a ) || (r2 <= a) ) */ 
    /*     return LBIG; */

    if ( (r1 <= a ) || (r2 <= a) ) 
        return LBIG;

    double d1 = deltaSeparation(r1,a);
    double d2 = deltaSeparation(r2,a);

    double t1 = -d1*d2/(2.0*constants()->lambda()*constants()->tau());

    /* communicate()->file("debug")->stream() << sep1[0] << "\t" << sep2[0] << "\t" */ 
    /*     << -log(1.0-exp(t1)) << endl; */

    return (-log(1.0 - exp(t1)));
}

/**************************************************************************//**
 *  The derivative of the effective potential with respect to lambda.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-17.
 *
 *  @param sep1 the first separation
 *  @param sep2 the second separation
 *  @return the derivative of the effective potential with respect to lambda
******************************************************************************/
double HardRodPotential::dVdlambda(const dVec &sep1, const dVec &sep2) 
{

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));
    double d1 = deltaSeparation(r1,a);
    double d2 = deltaSeparation(r2,a);

    double t1 = d1*d2;
    double t2 = t1/(2.0*constants()->lambda()*constants()->lambda()*constants()->tau());

    return ((0.5*t1/(exp(t2)-1.0))/(constants()->lambda()*constants()->lambda()*constants()->tau()));
}

/**************************************************************************//**
 *  The derivative of the effective potential with respect to tau.
 *
 *  Computes the non-local two-body effective pair potential.
 *
 *  Tested and working with Mathematica on 2013-06-17.
 *
 *  @param sep1 the first separation
 *  @param sep2 the second separation
 *  @return the derivative of the effective potential with respect to tau
******************************************************************************/
double HardRodPotential::dVdtau(const dVec &sep1, const dVec &sep2) 
{

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));
    double d1 = deltaSeparation(r1,a);
    double d2 = deltaSeparation(r2,a);

    double t1 = d1*d2;
    double t2 = t1/(2.0*constants()->lambda()*constants()->tau());

    return ((0.5*t1/(exp(t2)-1.0))/(constants()->lambda()*constants()->tau()*constants()->tau()));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GraphenePotential Class----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
#include <boost/math/special_functions/bessel.hpp>
GraphenePotential::GraphenePotential (double _strain, double _poisson, double _a0, 
        double _sigma, double _epsilon) : PotentialBase() {
    double strain = _strain;
    double poisson = _poisson;
    double a0 = _a0;

    sigma = _sigma;
    epsilon = _epsilon;

    /* Lattice vectors */
    a1x = (sqrt(3.)*a0/8.)*(4.+strain-(3.*strain*poisson));
    a1y = (3.*a0/8.)*(4.+(3.*strain)-(strain*poisson));
    a2x = -a1x;
    a2y = a1y;

    /* reciprocal lattice vectors */
    g1x = 8.*M_PI*sqrt(3.)*(4. + (3.*strain) - (strain*poisson))/(3.*a0*(4. + strain - (3.*strain*poisson))*(4.+(3*strain)-(strain*poisson)));
    g1y = 8.*M_PI*(4. + strain - (3.*strain*poisson))/(3.*a0*(4. + strain - (3.*strain*poisson))*(4.+(3*strain)-(strain*poisson)));
    g2x = -g1x;
    g2y = g1y;

    /* basis vectors */
    b1x = a2x;
    b1y = -0.5 * (1. + strain) * a0;
    b2x = 0.;
    b2y = -.5*(1. + strain) * a0-(a1x/sqrt(3.));

    /* area of unit cell */
    A = fabs((a1x*a2y) - (a1y*a2x));

}


/**************************************************************************//**
 * Destructor.
******************************************************************************/
GraphenePotential::~GraphenePotential() {

}

/**************************************************************************//**
 *  Return the value of the van der Waals' interaction between a graphene sheet
 *  and a helium adatom at a position, r, above the sheet. 

 *  @param r the position of a helium particle
 *  @return the van der Waals' potential for graphene-helium
******************************************************************************/
double GraphenePotential::V(const dVec &r) {

    double x = r[0];
    double y = r[1];
    double z = r[2]+(constants()->L()/2);
    double g = 0.;
    double gdotb1 = 0.;
    double gdotb2 = 0.;
    double k5term = 0.;
    double k2term = 0.;
    double prefactor = 0.;
    double v_g = 0.;
    double v = (4.*M_PI/A)*epsilon*sigma*sigma*( ((2./5.)*pow((sigma/z),10)) - pow((sigma/z),4) );
    if (z < 1.5) {
        return 10000.;
    }
    
    else {
        for (double m = -1; m < 1+1; m++) {
            for (double n = -1; n < 1+1; n++) {
                if ((m != 0) || (n != 0)){
                g = sqrt(pow((m*g1x + n*g2x),2.) + pow((m*g1y + n*g2y),2.));
                gdotb1 = ((m*g1x + n*g2x)*(b1x+x)) + ((m*g1y + n*g2y)*(b1y+y));
                gdotb2 = ((m*g1x + n*g2x)*(b2x+x)) + ((m*g1y + n*g2y)*(b2y+y));
                k5term = pow((g*sigma*sigma/2./z),5.)*boost::math::cyl_bessel_k(5., g*z)/30.;
                k2term = 2.*pow((g*sigma*sigma/2./z),2.)*boost::math::cyl_bessel_k(2., g*z);
                prefactor = epsilon*sigma*sigma*2.*M_PI/A;

                v_g = prefactor*(k5term-k2term);
                v += (cos(gdotb1)+cos(gdotb2))*v_g;
                }
            }
        }
        

        return v;

    }
}
/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * We create particles at random locations above the graphene sheet.
******************************************************************************/
Array<dVec,1> GraphenePotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {

    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;
    double initSideCube = 1.0*numParticles;
    
    for (int i = 0; i < NDIM - 1; i++) {
        initSideCube /= boxPtr->side[i];
    }

    initSideCube /= ((boxPtr->side[NDIM-1] / 2.0) - 6.0);

    /* Get the linear size per particle, and the number of particles */
    double initSide = pow((initSideCube),-1.0/(1.0*NDIM));

    /* We determine the number of initial grid boxes there are in 
     * in each dimension and compute their size */
    int totNumGridBoxes = 1;
    iVec numNNGrid;
    dVec sizeNNGrid;

    for (int i = 0; i < NDIM - 1; i++) {
        numNNGrid[i] = static_cast<int>(ceil((boxPtr->side[i] / initSide) - EPS));

        /* Make sure we have at least one grid box */
        if (numNNGrid[i] < 1)
            numNNGrid[i] = 1;

        /* Compute the actual size of the grid */
        sizeNNGrid[i] = boxPtr->side[i] / (1.0 * numNNGrid[i]);

        /* Determine the total number of grid boxes */
        totNumGridBoxes *= numNNGrid[i];
    }

    numNNGrid[NDIM-1] = static_cast<int>(ceil(((boxPtr->side[NDIM-1] - 12.0) / (2.0 * initSide)) - EPS));

    /* Make sure we have at least one grid box */
    if (numNNGrid[NDIM-1] < 1)
        numNNGrid[NDIM-1] = 1;

    /* Compute the actual size of the grid */
    sizeNNGrid[NDIM-1] = (boxPtr->side[NDIM-1] - 12.0) / (2.0 * numNNGrid[NDIM-1]);

    /* Determine the total number of grid boxes */
    totNumGridBoxes *= numNNGrid[NDIM-1];

    /* Now, we place the particles at the middle of each box */
    PIMC_ASSERT(totNumGridBoxes>=numParticles);
    dVec pos;
    for (int n = 0; n < totNumGridBoxes; n++) {

        iVec gridIndex;
        for (int i = 0; i < NDIM; i++) {
            int scale = 1;
            for (int j = i+1; j < NDIM; j++) 
                scale *= numNNGrid[j];
            gridIndex[i] = (n/scale) % numNNGrid[i];
        }

        for (int i = 0; i < NDIM - 1; i++) 
            pos[i] = (gridIndex[i]+0.5)*sizeNNGrid[i] - 0.5*boxPtr->side[i];
        
        pos[NDIM - 1] = (gridIndex[NDIM - 1]+0.5)*sizeNNGrid[NDIM - 1] + 3.0;

        boxPtr->putInside(pos);

        if (n < numParticles)
            initialPos(n) = pos;
        else 
            break;
    }
    return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GrapheneLUTPotential Class----------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
#include <boost/math/special_functions/bessel.hpp>
GrapheneLUTPotential::GrapheneLUTPotential (double _strain, double _poisson, double _a0, 
        double _sigma, double _epsilon, const Container *_boxPtr) : PotentialBase() {

    double strain = _strain;
    double poisson = _poisson;
    double a0 = _a0;
    sigma = _sigma;
    epsilon = _epsilon;

    /* get a local copy of the system size */
    Lzo2 = 0.5*_boxPtr->side[NDIM-1];

    /* Lookup Tables */
    tableLength = int((zmax - zmin)/dr);
    
    /* gradvg.resize(gtot, tableLength); */

    /* Lattice vectors */
    a1x = (sqrt(3.)*a0/8.)*(4.+strain-(3.*strain*poisson));
    a1y = (3.*a0/8.)*(4.+(3.*strain)-(strain*poisson));
    a2x = -a1x;
    a2y = a1y;

    /* reciprocal lattice vectors */
    g1x = 8.*M_PI*sqrt(3.)*(4. + (3.*strain) - (strain*poisson))/(3.*a0*(4. + strain - (3.*strain*poisson))*(4.+(3*strain)-(strain*poisson)));
    g1y = 8.*M_PI*(4. + strain - (3.*strain*poisson))/(3.*a0*(4. + strain - (3.*strain*poisson))*(4.+(3*strain)-(strain*poisson)));
    g2x = -g1x;
    g2y = g1y;

    /* basis vectors */
    b1x = a2x;
    b1y = -0.5 * (1. + strain) * a0;
    b2x = 0.;
    b2y = -.5*(1. + strain) * a0-(a1x/sqrt(3.));

    /* area of unit cell */
    A = fabs((a1x*a2y) - (a1y*a2x));
    
    double g = 0.0;
    double prefactor = epsilon*sigma*sigma*2.*M_PI/A;
    double k5term = 0.0;
    double k2term = 0.0;
    /* double dk1term = 0.0; */
    /* double dk2term = 0.0; */
    /* double dk3term = 0.0; */
    /* double dk4term = 0.0; */
    /* double dk5term = 0.0; */
    /* double dk6term = 0.0; */
    double z = 0.0;

    /* FULL LOOKUP TABLE OFFERS NO SPEED BENEFITS */
    /* vg.resize(2*(gnum+1)*gnum+1,tableLength); */

    /* for (int iz = 0; iz < tableLength; iz++) { */
    /*     z = zmin + iz*dr; */
    /*     vg(0,iz) = q*prefactor*( ((2./5.)*pow((sigma/z),10)) - pow((sigma/z),4) ); */
    /* } */

    /* int ig = 1; */
    /* for (int m = -gnum; m <= gnum; m++) { */
    /*     for (int n = 1; n <= gnum; n++) { */
    /*         g = sqrt(pow((m*g1x + n*g2x),2) + pow((m*g1y + n*g2y),2)); */
    /*         for (int iz = 0; iz < tableLength; iz++) { */
    /*             z = zmin + iz*dr; */
    /*             k5term = pow((g*sigma*sigma/2./z),5)*boost::math::cyl_bessel_k(5, g*z)/30.; */
    /*             k2term = 2.*pow((g*sigma*sigma/2./z),2)*boost::math::cyl_bessel_k(2, g*z); */
    /*             vg(ig,iz) = prefactor*(k5term-k2term); */
    /*         } */
    /*         ig++; */
    /*     } */
    /* } */
    /* for (int m=1; m <= gnum; m++) { */
    /*     g = sqrt(pow(m*g1x,2) + pow(m*g1y,2)); */
    /*     for (int iz = 0; iz < tableLength; iz++) { */
    /*         z = zmin + iz*dr; */
    /*         k5term = pow((g*sigma*sigma/2./z),5)*boost::math::cyl_bessel_k(5, g*z)/30.; */
    /*         k2term = 2.*pow((g*sigma*sigma/2./z),2)*boost::math::cyl_bessel_k(2, g*z); */
    /*         vg(ig,iz) = prefactor*(k5term-k2term); */
    /*     } */
    /*     ig++; */
    /* } */

    /* Create a unique list of all g-vector magnitudes */
    set<double> uniquegMag;

    /* This gives us the upper half-plane */
    uniquegMag.insert(0.0);
    for (int m = -gnum; m <= gnum; m++) {
        for (int n = 1; n <= gnum; n++) {
            g = sqrt(pow((m*g1x + n*g2x),2) + pow((m*g1y + n*g2y),2));
            uniquegMag.insert(g);
        }
    }
    for (int m=1; m <= gnum; m++) {
        g = sqrt(pow(m*g1x,2) + pow(m*g1y,2));
        uniquegMag.insert(g);
    }

    /* Convert the set to a vector and sort */
    std::vector<double> gMag(uniquegMag.begin(), uniquegMag.end());
    sort(gMag.begin(),gMag.end());

    /* Create the mapping from g-vectors to g-magnitudes */
    gMagID.resize(2*(gnum+1)*gnum+1);
    gMagID(0) = 0;

    int ig = 1;
    for (int m = -gnum; m <= gnum; m++) {
        for (int n = 1; n <= gnum; n++) {
            g = sqrt(pow((m*g1x + n*g2x),2) + pow((m*g1y + n*g2y),2));
            gMagID(ig) = distance(gMag.begin(), find(gMag.begin(), gMag.end(),g));
            ig++;
        }
    }
    for (int m=1; m <= gnum; m++){
        g = sqrt(pow(m*g1x,2) + pow(m*g1y,2));
        gMagID(ig) = distance(gMag.begin(), find(gMag.begin(), gMag.end(),g));
        ig++;
    }

    /* Fill up the lookup table based on g-magnitudes only */
    vg.resize(gMag.size(), tableLength);

    ig = 0;
    for (auto cg : gMag) {
        for (int iz = 0; iz < tableLength; iz++) {
            z = zmin + iz*dr;
            if (ig == 0)
                vg(ig,iz) = q*prefactor*( ((2./5.)*pow((sigma/z),10)) - pow((sigma/z),4) );
            else {
                k5term = pow((cg*sigma*sigma/2./z),5)*boost::math::cyl_bessel_k(5, cg*z)/30.;
                k2term = 2.*pow((cg*sigma*sigma/2./z),2)*boost::math::cyl_bessel_k(2, cg*z);
                vg(ig,iz) = prefactor*(k5term-k2term);
            }
        }
        ig++;
    }

    /* print out the potential to disk */
    /* int numPoints = 500; */
    /* double dx = _boxPtr->side[0]/numPoints; */
    /* double dy = _boxPtr->side[1]/numPoints; */
    /* dVec pos; */
    /* pos[2] = -0.5*_boxPtr->side[2] + 2.635; */
    /* for (int i = 0; i < numPoints; i++) { */
    /*     pos[0] = -0.5*_boxPtr->side[0] + i*dx; */
    /*     for (int j = 0; j < numPoints; j++) { */
    /*         pos[1] = -0.5*_boxPtr->side[1] + j*dy; */
    /*         communicate()->file("debug")->stream() << format("%24.16e %24.16e %24.16e\n") % pos[0] % pos[1] % V(pos); */
    /*     } */
    /* } */

    /* exit(-1); */

    /* print out the potential to disk */
    /* int numPoints = 500000; */
    /* double dz = 2*_boxPtr->side[2]/numPoints; */
    /* dVec pos; */
    /* pos[0] = 0.0; */
    /* pos[1] = 0.0; */
    /* for (int i = 0; i < numPoints; i++) { */
    /*     pos[2] = -0.5*_boxPtr->side[2] + 2.0 + i*dz; */
    /*     double cz = pos[2] + 0.5*_boxPtr->side[2]; */
    /*     communicate()->file("debug")->stream() << format("%24.16e %24.16e\n") % cz % V(pos); */
    /* } */

    /* exit(-1); */
}


/**************************************************************************//**
 * Destructor.
******************************************************************************/
GrapheneLUTPotential::~GrapheneLUTPotential() {
    gMagID.free();
    vg.free();
    gradvg.free();
}

/**************************************************************************//**
 *  Return the value of the van der Waals' interaction between a graphene sheet
 *  and a helium adatom at a position, r, above the sheet. 

 *  @param r the position of a helium particle
 *  @return the van der Waals' potential for graphene-helium
******************************************************************************/
double GrapheneLUTPotential::V(const dVec &r) {
    
    double z = r[2]+(Lzo2);
    int zindex = int((z-zmin)/dr);
    if (z < zmin) 
        return 400000.;
    
    if (zindex >= tableLength)
        return q*(epsilon*sigma*sigma*2.*M_PI/A)*( ((2./5.)*pow((sigma/z),10)) - pow((sigma/z),4) );

    double x1 = r[0]+b1x;
    double x2 = r[0]+b2x;
    double y1 = r[1]+b1y;
    double y2 = r[1]+b2y;

    double mx,my;
    double v = vg(0,zindex);

    int ig = 1;
    for (int m = -gnum; m <= gnum; m++) {
        mx = m*g1x;
        my = m*g1y;
        for (int n = 1; n <= gnum; n++) {
            v += 2.0*(cos((mx + n*g2x)*x1 + (my + n*g2y)*y1) + 
                      cos((mx + n*g2x)*x2 + (my + n*g2y)*y2)) * vg(gMagID(ig),zindex); 
            ig++;
        }
    }

    for (int m=1; m <= gnum; m++){
        mx = m*g1x;
        my = m*g1y;
        v += 2.0*(cos(mx*x1 + my*y1) + cos(mx*x2 + my*y2)) * vg(gMagID(ig),zindex);
        ig++;
    }

    return v;
}

/**************************************************************************//**
 *  Return the gradient of the van der Waals' interaction between a graphene sheet
 *  and a helium adatom at a position, r, above the sheet. 

 *  @param r the position of a helium particle
 *  @return the gradient of the van der Waals' potential for graphene-helium
******************************************************************************/
dVec GrapheneLUTPotential::gradV(const dVec &r) {

    return 0.0;

    /* double x = r[0]; */
    /* double y = r[1]; */
    /* double z = r[NDIM-1]+(Lzo2); */
    /* int zindex = int((z-zmin)/dr); */
    /* dVec dv = 0.0; */
    
    /* if (z < zmin) */ 
    /*     return dv; */

    /* if (zindex >= tableLength) { */
    /*     dv[NDIM-1] += q*(epsilon*sigma*sigma*2.*M_PI/A)* 4. * (pow(sigma/z,4) - pow(sigma/z,10)) / z; */
    /*     return dv; */
    /* } */
    
    /* dv[NDIM-1] = gradvg(0,zindex); */
    /* double gdotb1 = 0.; */
    /* double gdotb2 = 0.; */

    /* int k = 0; */
    
    /* /1* NB: this has not been optimized!!! *1/ */
    /* for (int m = -gnum; m < gnum+1; m++) { */
    /*     for (int n = -gnum; n < gnum+1; n++) { */
    /*         k = karr(m+gnum,n+gnum); */
    /*         if((m != 0) or (n != 0)) { */
    /*             gdotb1 = ((m*g1x + n*g2x)*(b1x+x)) + ((m*g1y + n*g2y)*(b1y+y)); */
    /*             gdotb2 = ((m*g1x + n*g2x)*(b2x+x)) + ((m*g1y + n*g2y)*(b2y+y)); */
                
    /*             dv[0] += -(g1x*m+g2x*n) * (sin(gdotb1) + sin(gdotb2))*vg(k,zindex); */
    /*             dv[1] += -(g1y*m+g2y*n) * (sin(gdotb1) + sin(gdotb2))*vg(k,zindex); */
    /*             dv[NDIM-1] += (cos(gdotb1)+cos(gdotb2))*gradvg(k,zindex); */
    /*         } */
    /*     } */
    /* } */
    /* return dv; */
}

/**************************************************************************//**
 * Return an initial particle configuration.
 *
 * We create particles at random locations above the graphene sheet.
******************************************************************************/
Array<dVec,1> GrapheneLUTPotential::initialConfig(const Container *boxPtr, MTRand &random,
        const int numParticles) {
    /* The particle configuration */
    Array<dVec,1> initialPos(numParticles);
    initialPos = 0.0;
    
    int Nlayer = round(boxPtr->side[0]*boxPtr->side[1]/a1x/a1y/4);
    int numX = round(boxPtr->side[0]/a1x/2);
    int numY = round(boxPtr->side[1]/a1y/2);
    double gridX = boxPtr->side[0]/numX;
    double gridY = boxPtr->side[1]/numY;
    int gridSize = 0;
    if (3*Nlayer < numParticles){
        gridSize = numParticles - 3*Nlayer;
    }
    int numInLayers = numParticles - gridSize;
    int layerNum = 0;
    int iShifted = 0;
    dVec pos;
    for (int i = 0; i < numInLayers; i++){
        layerNum = i/Nlayer;
        iShifted = i - (layerNum*Nlayer);
        pos[0] = ((iShifted % numX) + 0.5)*gridX - (boxPtr->side[0]/2);
        pos[1] = ((iShifted / numX) + 0.5)*gridY - (boxPtr->side[1]/2);
        pos[NDIM-1] = ((layerNum+1)*3.0)-(boxPtr->side[NDIM-1]/2);
        boxPtr->putInside (pos);
        initialPos(i) = pos;
    }
    if (gridSize) {
        int initCubePart = ceil(pow(1.0*gridSize,1.0/(1.0*NDIM)));
        dVec initSide;
        for (int i = 0; i < NDIM - 1; i++) {
            initSide[i] = boxPtr->side[i]/initCubePart;
        }

        initSide[NDIM-1] = (boxPtr->side[NDIM-1] - 12.0)/initCubePart;
    
        for (int n = 0; n < gridSize; n++) {
            pos[0] = (((n % initCubePart) + 0.5) * initSide[0]) - (boxPtr->side[0]/2);
            pos[1] = (((n / initCubePart) + 0.5) * initSide[1]) - (boxPtr->side[1]/2);
            pos[NDIM-1] = (((n / initCubePart / initCubePart) + 0.5) * initSide[NDIM-1]) - (boxPtr->side[NDIM-1]/2) + 12.0;
            boxPtr->putInside (pos);
            initialPos(n + numInLayers) = pos;
        }
    }
    return initialPos;
}
