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
#include <boost/program_options.hpp>

namespace po = boost::program_options;

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

		void initConstants(po::variables_map &);

		/* All the get methods */
		double T() const {return T_;}					///< Get temperature.
        double imagTimeLength() const { return imagTimeLength_;}    ///< Get the extent in imaginary time.
		double mu() const {return mu_;} 				///< Get chemical potential.
		double tau() const {return tau_;}				///< Get imaginary time step
		double m() const {return m_;}					///< Get mass.
		double lambda() const {return lambda_;}			///< Get lambda = hbar^2/(2mk_B)
		double rc() const {return rc_;}					///< Get potential cutoff
		double rc2() const {return rc2_;}				///< Get potential cutoff squared
		double C0() const {return C0_;}					///< Get worm factor C0
		double C() const {return C_;}					///< Get full worm constant
		double V() const {return V_;}					///< Get cell volume
		double L() const {return L_;}					///< Get maximum side length
		double comDelta() const {return comDelta_;}		///< Get center of mass shift
		double displaceDelta() const {return displaceDelta_;}		///< Get center of mass shift
        int virialWindow() const {return virialWindow_;}    ///< Get centroid virial window size.

		/** Get deBroglie wavelength */
		double dBWavelength() const {return dBWavelength_;} ///< Get deBroglie wavelength
		/** Get density matrix normalization factor */
		double rho0Norm(int M, int d) const {return ( pow((4.0*M_PI*lambda_*M*tau_),-0.5*d));}
		/** Get (4lambda/tau)^{-1} */
		double fourLambdaTauInv() const { return (0.25 / (lambda_ * tau_)); }

        /* Get the move attempt probability */
        double attemptProb(string type) {
            if (attemptProb_.count(type))
                return attemptProb_[type];
            else {
                cerr << "Get: Attempt probability for " << type << " does not exist!" << endl;
                exit(EXIT_FAILURE);
                return 0.0;
            }
        }

        /* Set the move attempt probability */
        void setAttemptProb(string type,double prob) {
            attemptProb_[type] = prob;
            /* if (attemptProb_.count(type)) */
            /*     attemptProb_[type] = prob; */
            /* else { */
            /*     cerr << "Set: Attempt probability for " << type << " does not exist!" << endl; */
            /*     exit(EXIT_FAILURE); */
            /* } */
        }

		bool restart() const {return restart_;}			                    ///< Get restart state
        bool wallClockOn() const {return wallClockOn_;}                     ///< Get wallclockOn
        uint32 wallClock() const {return wallClock_;}                       ///< Get wallclock limit
		bool canonical() const { return canonical_;}                        ///< Get ensemble
        bool window() const { return window_;}                              ///< Get window on/off
        bool startWithState() const { return startWithState_;}               ///< Are we starting from a state?
        int windowWidth() const { return windowWidth_;}                     ///< Get window 1/2 width
        bool gaussianEnsemble() const { return gaussianEnsemble_;}          ///< Get enesemble weight on/off
        double gaussianEnsembleSD() const { return gaussianEnsembleSD_;}    ///< Get enesemble weight standard dev.

		int Mbar() {return Mbar_;}						///< Get Mbar
		int b() {return b_;}							///< Get bisection level
		int numTimeSlices() {return numTimeSlices_;}	///< Get number of time slices
		int initialNumParticles() { return initialNumParticles_;}	///< Get initial number of particles
        int maxWind() { return maxWind_;}               ///< Get the maximum winding number sampled
		uint32 id() {return id_;}						///< Get simulation ID
		uint32 numEqSteps() {return numEqSteps_;}	///< Get the number of equilibration steps
        int numBroken() {return numBroken_;}            //< Get number of broken paths
        double spatialSubregion() {return spatialSubregion_;}          //< Get size of subregion
        bool spatialSubregionOn() {return spatialSubregionOn_;}           //< Get subregion on/off
        int Npaths() {return Npaths_;}                  //< Get number of paths

		string intPotentialType() const {return intPotentialType_;}	///< Get interaction potential type
		string extPotentialType() const {return extPotentialType_;}	///< Get external potential type
        string waveFunctionType() const {return waveFunctionType_;}	///< Get wave function type
        double endFactor() const {return endFactor_;}        ///< Get end factor
        string actionType() const {return actionType_;}	            ///< Get wave action type

        /* Trial wave funciton parameters */
        double R_LL_wfn() const {return R_LL_wfn_;}        ///< Get Lieb-Liniger length scale
        double k_LL_wfn() const {return k_LL_wfn_;}        ///< Get Lieb-Liniger wave number
    
		/* Set methods */
		void setmu(double _mu) {mu_ = _mu;}				///< Set the value of the chemical potential
		void setCoMDelta(double _comDelta) {comDelta_ = _comDelta;}	///< Set the CoM move size
		void setDisplaceDelta(double _displaceDelta) {displaceDelta_ = _displaceDelta;}	///< Set the displace move size

		/* Increment/decrement methods */
		void shiftC0(double frac) {C0_ += frac*C0_; getC();}		        ///< Shift the value of C0
		void getC() {C_ = C0_ / (1.0*Mbar_*numTimeSlices_*V_);}		        ///< Get the value of the worm constant 
		void shiftCoMDelta(double frac) {comDelta_ += frac*comDelta_; }		///< Shift the CoM move size
		void shiftDisplaceDelta(double frac) {displaceDelta_ += frac*displaceDelta_; }		///< Shift the displace move size
		void shiftmu (double frac) { mu_ += frac; }				        ///< Shift the chemical potential
        void incid() {++id_;}                                               ///< Increment the PIMCID by 1

        bool saveStateFiles() { return saveStateFiles_;}                              ///< Are we saving states every MC bin?

	protected:
		ConstantParameters();
		ConstantParameters(const ConstantParameters&);				///< Protected constructor
		ConstantParameters& operator= (const ConstantParameters&);	///< Overload Singleton equals

	private:
		double T_;				// Temperature [K]
		double imagTimeLength_;	// Temperature [K]
		double mu_;				// Chemical Potential [K]
		double tau_;			// Imaginary time step (tau = beta/numTimeSlices)
		double lambda_;			// hbar^2 / (2*mass*k_B) [K A^2]
		double m_;				// The particle mass
		double dBWavelength_;	// Thermal de Broglie wavelength
		double comDelta_;	    // The size of center of mass move
        double displaceDelta_;  // The size of the displace move shift

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
        int numBroken_;           // The number of broken paths
        double spatialSubregion_;      // The limits of the spatial sub region for EE
        bool spatialSubregionOn_;     // True if using a spatial subregion for EE
        int Npaths_;                // Number of paths used

		uint32 id_;				// The unique simulation ID
		uint32 numEqSteps_;		// Number of equilibration steps

		bool restart_;              // Are we restarting the simulation
        uint32 wallClock_;          // The wall clock limit in seconds
        bool wallClockOn_;          // Is the wall clock on?
        bool canonical_;            // Are we in the canonical ensemble?
        bool window_;               // Are we using a particle number window?
        bool startWithState_;       // Are we starting from a supplied state?
        int windowWidth_;           // Half width of particle number window
        bool gaussianEnsemble_;     // Are we using gaussian ensemble weight?
        double gaussianEnsembleSD_; // Standard deviation of ensemble weight

		string intPotentialType_;   // The type of interaction potential
		string extPotentialType_;   // The type of external potential
        string waveFunctionType_;   // The type of trial wave function
        double R_LL_wfn_;           // The length scale of the Lieb-Liniger wave function
        double k_LL_wfn_;           // The wave number of the Lieb-Liniger wave function
        double endFactor_;          // The multiplicative factor of the potential on end beads
        string actionType_;         // The type of action

        int virialWindow_;        // Window size for centroid virial estimator
        int maxWind_;             // The maximum winding number sampled

        bool saveStateFiles_;       // Are we saving a state file every MC bin?
		
        map <string,double> attemptProb_;	// The move attempt probabilities
};

/**************************************************************************//**
 * Global public access to the constants.
******************************************************************************/
inline ConstantParameters* constants() {
	ConstantParameters *temp = ConstantParameters::getInstance();
	return temp;
}

#endif

