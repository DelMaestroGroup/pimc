/**
 * @file wavefunction.cpp
 * @author Adrian Del Maestro
 * @date 04.09.2013
 *
 * @brief WaveFunction implementation.
 */

#include "wavefunction.h"
#include "path.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// WAVEFUNCTION BASE CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  Setup the path data members for the constant trial wavefunction.
******************************************************************************/
WaveFunctionBase::WaveFunctionBase (const Path &_path, LookupTable &_lookup,
        string _name) :
    name(_name),
	path(_path),
    lookup(_lookup)
{
    // empty constructor
}

/**************************************************************************//**
 *  Empty base constructor.
******************************************************************************/
WaveFunctionBase::~WaveFunctionBase() {
    // empty destructor
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SECH WAVEFUNCTION CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
SechWaveFunction::SechWaveFunction(const Path &_path,LookupTable &_lookup, string _name) :
	WaveFunctionBase(_path,_lookup,_name)
{
	/* Set the parameter to its optimized value */
	a = sqrt(0.5*M_PI);
}

/**************************************************************************//**
 * Destructor.
******************************************************************************/
SechWaveFunction::~SechWaveFunction() {
    // empty destructor
}

/**************************************************************************//**
 * The value of the trial wave function.
******************************************************************************/
double SechWaveFunction::PsiTrial(const int slice) {

    /* The cumulative value */
    double psiT = 1.0;
    beadLocator beadIndex;
    dVec pos;
    double r;

    for (int ptcl = 0; ptcl < path.numBeadsAtSlice(slice); ptcl++) {
        beadIndex = slice,ptcl;
        pos = path(beadIndex);
        r = sqrt(dot(pos,pos));
        psiT *= 1.0/cosh(a*r);
    }

    return psiT;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// JASTROW WAVEFUNCTION CLASS ---------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
******************************************************************************/
JastrowWaveFunction::JastrowWaveFunction(const Path &_path,LookupTable &_lookup, string _name) :
WaveFunctionBase(_path,_lookup,_name)
{
	/* Set the parameter to its optimized value */
	alpha = 19.0;
    beta = 0.12;
    //beta = 3.07;
}

/**************************************************************************//**
* Destructor.
******************************************************************************/
JastrowWaveFunction::~JastrowWaveFunction() {
    // empty destructor
}

/**************************************************************************//**
* The value of the 2-body trial wave function.
******************************************************************************/
double JastrowWaveFunction::PsiTrial(const double r) {
    
    double psiT = exp( -0.5*alpha/(1.0+beta*pow(r,5.0)));
    return psiT;
}

/**************************************************************************//**
* The derivative of psi over psi
******************************************************************************/
double JastrowWaveFunction::delPsiTrial(const double r) {
    
    double delpsiT = 2.5*alpha*beta*pow(r,4.0)/pow((1.0+beta*pow(r,5.0)),2.0);
    return delpsiT;
}

/**************************************************************************//**
* The derivative of psi over psi
******************************************************************************/
double JastrowWaveFunction::delSqPsiTrial(const double r) {
    
    double delSqPsiT = -1.25*alpha*beta*pow(r,3.0)*(-8.0+
                    beta*(4.0-5.0*alpha)*pow(r,5.0)+12.0*pow(beta,2.0)*pow(r,10.0) )/
                        ( pow((1.0+beta*pow(r,5.0)),4.0) );
    return delSqPsiT;
}

/**************************************************************************//**
* The value of the trial wave function.
******************************************************************************/
double JastrowWaveFunction::PsiTrial(const int slice) {
    
    /* The cumulative value */
    double psiT = 1.0;
    int numParticles = path.numBeadsAtSlice(slice);    
    dVec sep;						// The spatial separation between beads.
    double r;                       // Distance between beads
    beadLocator bead1,bead2;
	bead1[0] = bead2[0] = slice;
        
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        /* The loop over all other particles, to find the total interaction
         * potential */
        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);
            r = sqrt(dot(sep,sep));
            psiT *= PsiTrial(r);
        } // bead2
        
	} // bead1
    
    return psiT;
}

/**************************************************************************//**
* The value of the N-body trial wave function.
******************************************************************************/
double JastrowWaveFunction::gradSqPsiTrial(const int slice) {
    
    /* The cumulative value */
    double gradSqPsiT = 1.0;
    double delPsi12;
    int numParticles = path.numBeadsAtSlice(slice);
    dVec sep;						// The spatial separation between beads.
    double r;                       // Distance between beads
    beadLocator bead1,bead2,bead3;
	bead1[0] = bead2[0] = bead3[0] = slice;
    
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        /* The loop over all other particles, to find the total interaction
         * potential */
        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);
            r = sqrt(dot(sep,sep));
            gradSqPsiT += 2.0*delSqPsiTrial(r);
            delPsi12 = delPsiTrial(r);
            for (bead3[1] = bead2[1]+1; bead3[1] < numParticles; bead3[1]++) {
                sep = path.getSeparation(bead3,bead1);
                r = sqrt(dot(sep,sep));
                gradSqPsiT += 2.0*delPsi12*delPsiTrial(r);
            }// bead 3
        } // bead2
        
	} // bead1
    
    return gradSqPsiT;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// LiebLiniger WAVEFUNCTION CLASS --------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
 ******************************************************************************/
LiebLinigerWaveFunction::LiebLinigerWaveFunction(const Path &_path,LookupTable &_lookup, string _name) :
WaveFunctionBase(_path,_lookup,_name)
{
	/* Set the parameter to its optimized value */
    R = constants()->R_LL_wfn();
    k = constants()->k_LL_wfn();
}

/**************************************************************************//**
 * Destructor.
 ******************************************************************************/
LiebLinigerWaveFunction::~LiebLinigerWaveFunction() {
    // empty destructor
}

/**************************************************************************//**
 * The value of the 2-body trial wave function.
 ******************************************************************************/
double LiebLinigerWaveFunction::PsiTrial(const double r) {
    
    double psiT = 1.0;
    if(r<R)
        psiT = cos(k*(abs(r)-R));
    return psiT;
}

/**************************************************************************//**
 * The derivative of psi over psi
 ******************************************************************************/
double LiebLinigerWaveFunction::delPsiTrial(const double r) {
    
    return 0.0;
}

/**************************************************************************//**
 * The derivative of psi over psi
 ******************************************************************************/
double LiebLinigerWaveFunction::delSqPsiTrial(const double r) {
    
    return 0.0;
}

/**************************************************************************//**
* The value of the trial wave function.
******************************************************************************/
double LiebLinigerWaveFunction::PsiTrial(const int slice) {
    
    /* The cumulative value */
    double psiT = 1.0;
    Array <bool,1> doParticles(path.numBeadsAtSlice(slice));
    doParticles = true;

    dVec sep;						// The spatial separation between beads.
    double r;                       // Distance between beads
    beadLocator bead1,bead2;
	bead1[0] = bead2[0] = slice;
    
    /* No cutoff */
//    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
//        /* The loop over all other particles, to find the total interaction
//         * potential */
//        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
//            sep = path.getSeparation(bead2,bead1);
//            r = sqrt(dot(sep,sep));
//            psiT *= PsiTrial(r);
//        } // bead2
//        
//	} // bead1
    
    /* Using cutoff */
	for (bead1[1] = 0; bead1[1] < path.numBeadsAtSlice(slice); bead1[1]++) {
        
		doParticles(bead1[1]) = false;
        
		/* Get the interaction list */
		lookup.updateInteractionList(path,bead1);
        
		/* Sum the interaction potential over all NN beads */
		for (int n = 0; n < lookup.numBeads; n++) {
			bead2 = lookup.beadList(n);
			if (doParticles(bead2[1])) {
				sep = path.getSeparation(bead2,bead1);
                r = sqrt(dot(sep,sep));
                psiT *= PsiTrial(r);
			}
		} // n
        
	} // bead1
    
    return psiT;
}

/**************************************************************************//**
* The weight of the trial wave function for a bead.
******************************************************************************/
double LiebLinigerWaveFunction::PsiTrial(const beadLocator &bead1) {
    
    /* The cumulative value */
    double psiT = 1.0;
    
    dVec sep;						// The spatial separation between beads.
    double r;                       // Distance between beads
    
    /* No cutoff */
    //    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
    //        /* The loop over all other particles, to find the total interaction
    //         * potential */
    //        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
    //            sep = path.getSeparation(bead2,bead1);
    //            r = sqrt(dot(sep,sep));
    //            psiT *= PsiTrial(r);
    //        } // bead2
    //
    //	} // bead1
    
	/* We only continue if bead1 is turned on */
	if (path.worm.beadOn(bead1)) {
        
		/* Fill up th nearest neighbor list */
		lookup.updateInteractionList(path,bead1);
        
		/* Sum the interaction potential over all NN beads */
		for (int n = 0; n < lookup.numBeads; n++) {
            
            sep = path.getSeparation(bead1,lookup.beadList(n));
            r = sqrt(dot(sep,sep));
            psiT *= PsiTrial(r);
		}
	}
    
    return psiT;
}


/**************************************************************************//**
 * The value of the N-body trial wave function.
 ******************************************************************************/
double LiebLinigerWaveFunction::gradSqPsiTrial(const int slice) {
    
    return 0.0;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SUTHERLAND WAVEFUNCTION CLASS ---------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 * Constructor.
 ******************************************************************************/
SutherlandWaveFunction::SutherlandWaveFunction(const Path &_path, LookupTable &_lookup, 
        double _lambda, string _name) :
WaveFunctionBase(_path,_lookup,_name)
{
    // The Sutherland model value of the interaction paramter \lambda
    lambda = _lambda;

    // pi/L
    pioL = M_PI/constants()->L();

    // The wavefunction pre-factor
    // N.B. This only works for canonical simulations 
    int N = constants()->initialNumParticles();
    CN = pow(tgamma(1.0+lambda)/constants()->L(),0.5*N) / sqrt(tgamma(1.0 + lambda*N));
    CN *= pow(2.0,0.5*N*(N-1)*lambda);
}

/**************************************************************************//**
 * Destructor.
 ******************************************************************************/
SutherlandWaveFunction::~SutherlandWaveFunction() {
    // empty destructor
}

/**************************************************************************//**
* The value of the trial wave function.
******************************************************************************/
double SutherlandWaveFunction::PsiTrial(const int slice) {
    
    /* The cumulative value */
    double psiT = 1.0;

    /* The number of particles */
    int numParticles = path.numBeadsAtSlice(slice);    

    dVec sep;						// The spatial separation between beads.
    double r;                       // Distance between beads
    beadLocator bead1,bead2;
	bead1[0] = bead2[0] = slice;
    
    /* No cutoff */
    for (bead1[1] = 0; bead1[1] < numParticles; bead1[1]++) {
        /* The loop over all other particles, to find the total wavefunction */
        for (bead2[1] = bead1[1]+1; bead2[1] < numParticles; bead2[1]++) {
            sep = path.getSeparation(bead2,bead1);
            r = sqrt(dot(sep,sep));
            psiT *= PsiTrial(r);
        } // bead2
    } // bead1

    return CN*psiT;
}
