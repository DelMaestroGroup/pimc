/**
 * @file constants.h
 * @author Adrian Del Maestro
 * @date 03.17.2009	St Patrick's Day
 *
 * @brief ConstantParameters class definition.
 */

#ifndef CONSTANTS_H 
#define CONSTANTS_H

#include "common.h"

// ========================================================================  
// ConstantParameters Class
// ========================================================================  
/** 
 * Constant simulation parameters.
 *
 * Holds a number of parameters that will need to be accessed throughout
 * the simulation and allows their access via the singleton design
 * pattern.
 * @see http://en.wikipedia.org/wiki/Singleton_pattern
 */
class ConstantParameters
{
	public:
		static ConstantParameters* getInstance();
		void initConstants(bool,double,double,double,double,double,double,double,
				int,int,int,uint32,uint32,uint32,string,string);

		/* All the get methods */
		double T() const {return T_;}					///< Get temperature.
		double mu() const {return mu_;} 				///< Get chemical potential.
		double tau() const {return tau_;}				///< Get imaginary time step
		double m() const {return m_;}					///< Get mass.
		double lambda() const {return lambda_;}			///< Get lambda = hbar^2/(2mk_B)
		double rc() const {return rc_;}					///< Get potential cutoff
		double rc2() const {return rc2_;}				///< Get potential cutoff squared
		double C0() const {return C0_;}					///< Get worm factor C0
		double Delta() const {return Delta_;}			///< Get center of mass fhit
		double C() const {return C_;}					///< Get full worm constant
		double V() const {return V_;}					///< Get cell volume
		double L() const {return L_;}					///< Get maximum side length

		/** Get deBroglie wavelength */
		double dBWavelength() const {return dBWavelength_;} ///< Get deBroglie wavelength
		/** Get density matrix normalization factor */
		double rho0Norm(int M, int d) const {return ( pow((4.0*M_PI*lambda_*M*tau_),-0.5*d));}
		/** Get (4lambda/tau)^{-1} */
		double fourLambdaTauInv() const { return (0.25 / (lambda_ * tau_)); }

		/* Get methods for the attempt probabilities */
		double openAttemptProb() const {return openAttemptProb_;}			///< Get open probability
		double closeAttemptProb() const {return closeAttemptProb_;}			///< Get close probability
		double insertAttemptProb() const {return insertAttemptProb_;}		///< Get insert probability
		double removeAttemptProb() const {return removeAttemptProb_;}		///< Get remove probability
		double advanceHeadAttemptProb() const {return advanceHeadAttemptProb_;}	///< Get advance head probability
		double recedeHeadAttemptProb() const {return recedeHeadAttemptProb_;}	///< Get recede head probabililty
		double advanceTailAttemptProb() const {return advanceTailAttemptProb_;}	///< Get advance tail probabililty
		double recedeTailAttemptProb() const {return recedeTailAttemptProb_;}	///< Get recede tail probabililty
		double swapHeadAttemptProb() const {return swapHeadAttemptProb_;}		///< Get swap head probabililty
		double swapTailAttemptProb() const {return swapTailAttemptProb_;}		///< Get swap tail probabililty
		double comAttemptProb() const {return comAttemptProb_;}					///< Get CoM probabililty
		double stagingAttemptProb() const {return stagingAttemptProb_;}			///< Get staging probabililty

		bool restart() const {return restart_;}			///< Get restart state
		bool canonical() const { return canonical_;}	///< Get ensemble

		int Mbar() {return Mbar_;}						///< Get Mbar
		int b() {return b_;}							///< Get bisection level
		int numTimeSlices() {return numTimeSlices_;}	///< Get number of time slices
		int initialNumParticles() { return initialNumParticles_;}	///< Get initial number of particles
		int deltaNumParticles() { return deltaNumParticles_;}	///< Get the canonical fluctuation weight
		uint32 id() {return id_;}						///< Get simulation ID
		uint32 numEqSteps() {return numEqSteps_;}	///< Get the number of equilibration steps


		string intPotentialType() const {return intPotentialType_;}	///< Get interaction potential type
		string extPotentialType() const {return extPotentialType_;}	///< Get external potential type

		/* Set methods */
		void setmu(double _mu) {mu_ = _mu;}				///< Set the value of the chemical potential
		void setDelta(double _Delta) {Delta_ = _Delta;}	///< Set the CoM move size

		/* Increment/decrement methods */
		void shiftC0(double frac) {C0_ += frac*C0_; getC();}		///< Shift the value of C0
		void getC() {C_ = C0_ / (1.0*Mbar_*numTimeSlices_*V_);}		///< Get the value of the worm constant 
		void shiftDelta(double frac) {Delta_ += frac*Delta_; }		///< Shift the CoM move size
		void shiftmu (double frac) { mu_ += frac*mu_; }				///< Shift the chemical potential
        void incid() {++id_;}                                       ///< Increment the PIMCID by 1

	protected:
		ConstantParameters();
		ConstantParameters(const ConstantParameters&);				///< Protected constructor
		ConstantParameters& operator= (const ConstantParameters&);	///< Overload Singleton equals

	private:
		double T_;				// Temperature [K]
		double mu_;				// Chemical Potential [K]
		double tau_;			// Imaginary time step (tau = beta/numTimeSlices)
		double lambda_;			// hbar^2 / (2*mass*k_B) [K A^2]
		double m_;				// The particle mass
		double dBWavelength_;	// Thermal de Broglie wavelength
		double Delta_;			// The size of center of mass and displace moves

		double rc_;				// The potential cutoff length
		double rc2_;			// The potential cutoff length squared
		double C0_;				// Worm probability constant [PRE 74, 036701 (2006)]
		double C_;				// Worm normalization factor
		double V_;				// The volume of the system
		double L_;				// The linear 'length' of the system 

		int Mbar_;				  // Maximum worm algorithm trial length
		int b_;					  // The maximum number of levels
		int numTimeSlices_;		  // Number of imaginary time slices
		int initialNumParticles_; // The initial number of particles
		int deltaNumParticles_;   // The particle number fluctuation weight (for canonical sims)

		uint32 id_;				// The unique simulation ID
		uint32 numEqSteps_;		// Number of equilibration steps

		bool restart_;			// Are we restarting the simulation
		bool canonical_;		// Are we in the canonical ensemble?

		string intPotentialType_; // The type of interaction potential
		string extPotentialType_; // The type of external potential

		/* All the move attempt probabilities */
		double openAttemptProb_;
		double closeAttemptProb_;
		double insertAttemptProb_;
		double removeAttemptProb_;
		double advanceHeadAttemptProb_;
		double recedeHeadAttemptProb_;
		double advanceTailAttemptProb_;
		double recedeTailAttemptProb_;
		double swapHeadAttemptProb_;
		double swapTailAttemptProb_;
		double comAttemptProb_;
		double stagingAttemptProb_;
};

/**************************************************************************//**
 * Global public access to the constants.
******************************************************************************/
inline ConstantParameters* constants() {
	ConstantParameters *temp = ConstantParameters::getInstance();
	return temp;
}

#endif

