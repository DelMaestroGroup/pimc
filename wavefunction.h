/**
 * @file wavefunction.h
 * @author Adrian Del Maestro
 * @date 04.09.2013
 *
 * @brief Action class definitions.
 */

#include "constants.h"

#ifndef WAVEFUNCTION_H 
#define WAVEFUNCTION_H

class Path;

// ========================================================================  
// WaveFunctionBase Class
// ========================================================================  
/** 
 * Holds a base class that all trial wave function classes will be derived from. 
 *
 * For a given set of worldlines, computes the real value of a trial
 * wavefunction.
 */
class WaveFunctionBase {

	public:
		WaveFunctionBase (const Path &, string _name="constant");
		virtual ~WaveFunctionBase();

        /** Name of the trial wave function */
        string name;

		/** The Constant Trial Wave Function*/
		virtual double psiTrial(const double r) {
            return 1.0;
        };
        virtual double PsiTrial(const int slice) {
            return 1.0;
        };
        virtual double delPsiTrial(const double r) {
            return 0.0;
        };
        virtual double delSqPsiTrial(const double r) {
            return 0.0;
        };
        virtual double gradSqPsiTrial(const int slice) {
            return 0.0;
        };

	protected:
		const Path &path;				///< A reference to the paths

};


// ========================================================================  
// SechWaveFunction Class
// ========================================================================  
/**
 * Implementation of the \Psi_T = sech(a*x) trial wave function suitable
 * for the simple harmonic osscilator.
 * @see A. Sarsa, K. E. Schmidt, and W. R. Magro, J. Chem. Phys. 113, 1366 (2000).
 */
class SechWaveFunction: public WaveFunctionBase {

	public:
		SechWaveFunction(const Path &, string _name="SHO sech");
		~SechWaveFunction();

        double psiTrial (const int);

	private:
		double a;			// The parameter of the wave function 

};


// ========================================================================
// JastrowWaveFunction Class
// ========================================================================
/**
 * Implementation of a Jastrow trial wave function suitable for He
 * @see Cuervo, Roy, & Boninsegni,, J. Chem. Phys., 122(11), 114504 (2005).
 */
class JastrowWaveFunction: public WaveFunctionBase {
    
public:
    JastrowWaveFunction(const Path &, string _name="Jastrow");
    ~JastrowWaveFunction();
    
    double psiTrial (const double);
    double PsiTrial (const int);
    
    double delPsiTrial(const double r);
    double delSqPsiTrial(const double r);
    double gradSqPsiTrial(const int);
    
    
private:
    double alpha;			// The parameter of the wave function
    double beta;			// The parameter of the wave function
    
};


#endif
