/** 
 * @file constants.cpp
 * @author Adrian Del Maestro
 *
 * @brief ConstantParameters class implementation.
 */

#include "constants.h"
#include <time.h>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CONSTANT PARAMETERS  CLASS ------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**************************************************************************//**
 *  An empty constructor which simply sets all constants to null.
******************************************************************************/
ConstantParameters::ConstantParameters() : T_(), mu_(), tau_(), lambda_(), m_(),
	dBWavelength_(), rc_(), C0_(), C_(), V_(), L_(), Mbar_(), b_(), numTimeSlices_(), 
	initialNumParticles_(), deltaNumParticles_(), id_(), restart_(), canonical_(), 
	intPotentialType_(), extPotentialType_()
{ }

/**************************************************************************//**
 *  Initialize all constants from passed values.
 *
 *  We initialize all constant parameters used in the simulation.  The ID is
 *  computed as the number of seconds since January 01, 2009.  The value of
 *	lambda = hbar^2/2 m k_B is computed in units where lenghts are measured in 
 *	angstroms and energies in kelvin.
******************************************************************************/
void ConstantParameters::initConstants(bool _canonical, double _T, double _mu, double _m,
		double _rc, double _C0, double _V, double _L, int _initialNumParticles, 
		int _Mbar, int _numTimeSlices, uint32 _id, uint32 _process,
		uint32 _numEqSteps, string _intPotentialType, string _extPotentialType) {

	/* The simulation ID is the number of seconds since January 1 2009 */
	if (_id == 0) {
		time_t seconds = (time(NULL) - 39*365*24*60*60);
		id_ = uint32(seconds + _process);
		restart_ = false;
	}
	else {
		id_ = _id;
		restart_ = true;
	}

	/* Are we working in the grand canonical ensemble? */
	canonical_ = _canonical;

	/* Assigned values */
	b_             = int (ceil(log(1.0*_Mbar) / log(2.0)-EPS));
	Mbar_          = _Mbar;
	T_             = _T;
	mu_            = _mu;
	m_             = _m;
	lambda_        = 24.24 / m_;
	rc_            = _rc;
	rc2_           = rc_*rc_;
	C0_            = _C0;
	numTimeSlices_ = _numTimeSlices;
	tau_           = 1.0/(numTimeSlices_*T_);
	V_	           = _V;
	L_             = _L;
	numEqSteps_    = _numEqSteps;

	initialNumParticles_ = _initialNumParticles;

	/* We arbitrarily set the particle weighting number (for now) */
	deltaNumParticles_ = 10;

	intPotentialType_ = _intPotentialType;
	extPotentialType_ = _extPotentialType;

	PIMC_ASSERT(Mbar_ < numTimeSlices_);
	PIMC_ASSERT(Mbar_ > 1);

	if (Mbar_ >= numTimeSlices_) {
		cout << "Update length is too long" << endl;
		exit(EXIT_FAILURE);
	}

	if (Mbar_ <= 1) {
		cout << "Update length is too short" << endl;
		exit(EXIT_FAILURE);
	}

	/* Computed values */
	dBWavelength_ = 2.0*sqrt(M_PI * lambda_ / T_);
	Delta_        = 0.04*dBWavelength_;
	getC();

	/* Here we set all the attempt probabilities */
    attemptProb_["open"] = 0.4;
    attemptProb_["insert"] = 0.4;
    attemptProb_["close"] = 0.15;
    attemptProb_["advance head"] = 0.075;
    attemptProb_["recede head"] = 0.075;
    attemptProb_["advance tail"] = 0.075;
    attemptProb_["recede tail"] = 0.075;
    attemptProb_["remove"] = 0.15;
    attemptProb_["swap head"] = 0.10;
    attemptProb_["swap tail"] = 0.10;
    attemptProb_["staging"] = 0.00;
    attemptProb_["bisection"] = 0.15;
    attemptProb_["center of mass"] = 0.05;

    double totProb = attemptProb_["close"] + attemptProb_["advance head"] + attemptProb_["recede head"]
        + attemptProb_["advance tail"] + attemptProb_["recede head"] + attemptProb_["remove"]
        + attemptProb_["swap head"] + attemptProb_["swap tail"] + attemptProb_["staging"]
        + attemptProb_["bisection"] + attemptProb_["center of mass"];

	if (abs(totProb - 1.0) > EPS) {
		cout << "Close + AdvanceHead + RecedeHead + AdvanceTail + RecedeTail + Remove + SwapHead " 
			 << "+ Staging + CoM Probability != 1" << endl;
		exit(EXIT_FAILURE);
	}
	PIMC_ASSERT(totProb-1.0 < EPS);

    totProb = attemptProb_["open"] + attemptProb_["insert"] + attemptProb_["staging"]
       + attemptProb_["bisection"] + attemptProb_["center of mass"];
	
	if (abs(totProb - 1.0) > EPS) {
		cout << "Open + Insert + Staging + Bisection + CoM Probability != 1" << endl;
		exit(EXIT_FAILURE);
	}
	PIMC_ASSERT(totProb-1.0 < EPS);
}

/**************************************************************************//**
 *  This public method returns an instance of the constant object,  Only one
 *  can ever exist at a time.
******************************************************************************/
ConstantParameters* ConstantParameters::getInstance ()
{   
	static ConstantParameters inst;
	return &inst;
}
