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
}

/**************************************************************************//**
 * Destructor. 
******************************************************************************/
TabulatedPotential::~TabulatedPotential() {
	lookupV.free();
	lookupdVdr.free();
}

/**************************************************************************//**
 *  Given a discretization factor and the system size, create and fill
 *  the lookup tables for the potential and its derivative.
******************************************************************************/
void TabulatedPotential::initLookupTable(const double _dr, const double maxSep) {

	/* We now calculate the lookup tables for the interaction potential and 
	 * its derivative. */
	dr = _dr;
	tableLength = int(maxSep/dr);
	lookupV.resize(tableLength);
	lookupdVdr.resize(tableLength);
	lookupV = 0.0;
	lookupdVdr = 0.0;

	double r = 0;

	for (int n = 0; n < tableLength; n++) {
		lookupV(n)    = valueV(r);
		lookupdVdr(n) = valuedVdr(r);
		r += dr;
	}

//	double rc = constants()->rc();
//	for (int n = 0; n < tableLength; n++) {
//		r += dr;
//		if (r <= rc) {
//			lookupV(n) = valueV(r) - valueV(rc) - valuedVdr(rc)*(r-rc);
//			lookupdVdr(n) = valuedVdr(r) - valuedVdr(rc);
//		}
//		else {
//			lookupV(n) = 0.0;
//			lookupdVdr(n) = 0.0;
//		}
//		cout << format("%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E\n") % r % lookupV(n) % valueV(r) % 
//			lookupdVdr(n) % valuedVdr(r) % (lookupV(n) - valueV(r)) % (lookupdVdr(n) - valuedVdr(r));
//	}

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
HarmonicPotential::HarmonicPotential() : PotentialBase() {
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
 *  @param _a The width^2 of the Gaussian
 *  @param _c The integrated strength of the Gaussian
******************************************************************************/
DeltaPotential::DeltaPotential(double _a, double _c) : PotentialBase() 
{
	/* Define the parameters of the delta function potential. */
	a = _a;
	c = _c;
	norm = c/sqrt(a*M_PI);
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

	char state;				// Fixed or updateable?
	dVec pos;				// The loaded position

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

	int locNumParticles = 0;	// Number of updateable particles
	char state;					// Update or Fix
	dVec pos;					// The current position

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
//	density = 0.008; // atoms / angstrom^3

	/* We define the values of epsilon and sigma for N and He */ 
//	double epsilonHe = 10.216; 	// Kelvin
//	double sigmaHe   = 2.556; 	// angstroms
//	double sigmaN    = 3.299; 	// angstroms
//	double epsilonN  = 36.2;	// Kelvin
//	epsilon = sqrt(epsilonHe*epsilonN);
//	sigma = 0.5*(sigmaHe + sigmaN);

	/* We operate under the assumption that the silicon can be neglected in the 
	 * silicon-nitride, and thus only consider the Nitrogen.  We use a
	 * Kiselov type model to extract the actual parameters.  We assume that
	 * silicate and silicon-nitride are roughly equivalent. */
	epsilon = 10.22; 	// Kelvin
	sigma   = 2.628;	// angstroms

//	epsilon = 32; 	// Kelvin
//	sigma   = 3.08;	// angstroms

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

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
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
	lside[0] = lside[1] = sqrt(2.0)*(R-1.0);
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
	rm      = 2.9673;  	// A
	A       = 0.5449E6; 
	epsilon = 10.8; 	// K
	alpha   = 13.353; 
	D       = 1.2413; 
	C6      = 1.3732;
	C8      = 0.42538;
	C10     = 0.1781;

	/* The extremal values are all zero here */
	extV = 0.0;
	extdVdr = 0.0;

	/* We take the maximum possible separation */
	double L = max(side);

	/* Create the potential lookup tables */
	//initLookupTable(0.00005*rm,L);
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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GASPARINI_1_POTENTIAL CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
/**************************************************************************//**
 * Constructor.
******************************************************************************/
Gasparini_1_Potential::Gasparini_1_Potential(double _az, double _ay, const Container *_boxPtr) :
    PotentialBase(),
    az(0.5*_boxPtr->side[2]*_az),
    ay(0.5*_boxPtr->side[1]*_ay),
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

	/* The particle configuration */
	Array<dVec,1> initialPos(numParticles);

	/* label the lengths of the sides of the simulation cell */
    dVec lside;
    lside[0] = boxPtr->side[0];
    lside[1] = boxPtr->side[1];
    lside[2] = boxPtr->side[NDIM-1];

    /* calculate actual volume */
    double volTot = product(lside);

    /* get linear size per particle for actual volume */
    double initSide = pow((1.0*numParticles/volTot),-1.0/(1.0*NDIM));

	/* For the ENTIRE SIMULATION CELL (even excluded volume) 
     * We determine the number of initial grid boxes there are in 
	 * each dimension and compute their size */
	int totNumGridBoxes = 1;
	iVec numNNGrid;
	dVec sizeNNGrid;

	for (int i = 0; i < NDIM; i++) {
		numNNGrid[i] = static_cast<int>(ceil((boxPtr->side[i] / initSide) - EPS));

		/* Make sure we have at least one grid box */
		if (numNNGrid[i] < 1)
			numNNGrid[i] = 1;

		/* Compute the actual size of the grid */
		sizeNNGrid[i] = lside[i] / (1.0 * numNNGrid[i]);

		/* Determine the total number of grid boxes */
		totNumGridBoxes *= numNNGrid[i];
	}

  	/* Now, we place the particles at the middle of each box
     * UNLESS that middle falls within our excluded volume */
	PIMC_ASSERT(totNumGridBoxes>=numParticles);
	dVec pos;
    int numIn = 0;
    double jump = 0.5;
    
    while (numIn < numParticles){
        for (int n = 0; n < totNumGridBoxes; n++) {

            /* put correct number of particles in box, no more */
            if (numIn >= numParticles)
                break;
            
            iVec gridIndex;
            /* update grid index */
            for (int i = 0; i < NDIM; i++) {
                int scale = 1;
                for (int j = i+1; j < NDIM; j++) 
                    scale *= numNNGrid[j];
                gridIndex[i] = (n/scale) % numNNGrid[i];
            }
            /* place particle in position vector, skipping over excluded volume */
            for (int i = 0; i < NDIM; i++)
                pos[i] = (gridIndex[i]+jump)*sizeNNGrid[i] - 0.5*lside[i];

            if (((pos[2] < -az) || (pos[2] > az)) && ((pos[1] < -ay) || (pos[1] > ay))){
                initialPos(numIn) = pos;
                numIn++;
            }
        }
        jump += 0.1;
    }

	return initialPos;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARD SPHERE POTENTIAL CLASS -----------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *  @param a The radius of the hard sphere (also the scattering length)
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
 *  @param lambdaTau The product of \lambda = \hbar^2/2m and \tau
 *  @return the two-body effective pair potential
******************************************************************************/
double HardSpherePotential::V(const dVec &sep1, const dVec &sep2, 
        double lambdaTau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double cosTheta = dot(sep1,sep2)/(r1*r2);


    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*lambdaTau);

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
 *  @param lambda \lambda = \hbar^2/2m
 *  @param tau the imaginary timestep tau
 *  @return the derivative of the effective potential with respect to lambda
******************************************************************************/
double HardSpherePotential::dVdlambda(const dVec &sep1, const dVec &sep2, 
        double lambda, double tau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double cosTheta = dot(sep1,sep2)/(r1*r2);

    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*lambda*tau);

    double t2 = (a*(r1+r2) - a*a)/(r1*r2);
    double t3 = 1.0 - t2*exp(t1);

    double t4 = (1.0/t3)*(t2*exp(t1))*(t1/lambda);

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
double HardSpherePotential::dVdtau(const dVec &sep1, const dVec &sep2, 
        double lambda, double tau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double cosTheta = dot(sep1,sep2)/(r1*r2);

    double t1 = -(r1*r2 + a*a - a*(r1 + r2)) * (1.0 + cosTheta);
    t1 /= (4.0*lambda*tau);

    double t2 = (a*(r1+r2) - a*a)/(r1*r2);
    double t3 = 1.0 - t2*exp(t1);

    double t4 = (1.0/t3)*(t2*exp(t1))*(t1/tau);

    return -t4;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HARD ROD POTENTIAL CLASS --------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Constructor.
 *  @param a The radius of the hard rod (also the scattering length)
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
 *  @param lambdaTau The product of \lambda = \hbar^2/2m and \tau
 *  @return the two-body effective pair potential
******************************************************************************/
double HardRodPotential::V(const dVec &sep1, const dVec &sep2, 
        double lambdaTau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double t1 = sep1[0]*sep2[0] - r1*r2;
    t1 /= (4.0*lambdaTau);

    double t2 = -(r1-a)*(r2-a)/(2.0*lambdaTau);
    double t3 = -log(1.0 - exp(t2));

    return 0.0*t1 + t3;
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
 *  @param lambda \lambda = \hbar^2/2m
 *  @param tau the imaginary timestep tau
 *  @return the derivative of the effective potential with respect to lambda
******************************************************************************/
double HardRodPotential::dVdlambda(const dVec &sep1, const dVec &sep2, 
        double lambda, double tau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double t1 = r1*r2 - sep1[0]*sep2[0];
    double t2 = (r1-a)*(r2-a);
    double t3 = t2/(2.0*lambda*tau);

    return ((0.0*t1 + 2.0*t2/(exp(t3)-1.0))/(4.0*lambda*lambda*tau));
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
 *  @param lambda \lambda = \hbar^2/2m
 *  @param tau the imaginary timestep tau
 *  @return the derivative of the effective potential with respect to tau
******************************************************************************/
double HardRodPotential::dVdtau(const dVec &sep1, const dVec &sep2, 
        double lambda, double tau) {

    double r1 = sqrt(dot(sep1,sep1));
    double r2 = sqrt(dot(sep2,sep2));

    if (r1*r2 < EPS)
        return log(BIG);

    double t1 = r1*r2 - sep1[0]*sep2[0];
    double t2 = (r1-a)*(r2-a);
    double t3 = t2/(2.0*lambda*tau);

    return ((0.0*t1 + 2.0*t2/(exp(t3)-1.0))/(4.0*lambda*tau*tau));
}
