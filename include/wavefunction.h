/**
 * @file wavefunction.h
 * @author Adrian Del Maestro
 * @date 04.09.2013
 *
 * @brief Action class definitions.
 */

/* #include "common.h" */

#include "constants.h"
#include "lookuptable.h"

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
        WaveFunctionBase (const Path &,LookupTable &_lookup);
        virtual ~WaveFunctionBase();

        /** The Constant Trial Wave Function*/
        virtual double PsiTrial(const int) { return 1.0; };
        virtual double PsiTrial(const double) { return 1.0; };
        virtual double PsiTrial(const beadLocator &) { return 1.0; };
        virtual double delPsiTrial(const double ) { return 0.0; };
        virtual double delSqPsiTrial(const double ) { return 0.0; };
        virtual double gradSqPsiTrial(const int) { return 0.0; };

    protected:
        const Path &path;               ///< A reference to the paths
        LookupTable &lookup;            ///< We need a non-constant reference for updates

};

// ========================================================================  
// ConstantWaveFunction Class
// ========================================================================  
/**
 * Simple constant wavefunction (implements the base class).
 */
class ConstantWaveFunction: public WaveFunctionBase {
    public:
        ConstantWaveFunction(const Path &,LookupTable &_lookup);
        ~ConstantWaveFunction();

        static const std::string name;
        std::string getName() const {return name;}
};


// ========================================================================  
// SechWaveFunction Class
// ========================================================================  
/**
 * Implementation of the Psi_T = sech(a*x) trial wave function suitable
 * for the simple harmonic osscilator.
 * @see A. Sarsa, K. E. Schmidt, and W. R. Magro, J. Chem. Phys. 113, 1366 (2000).
 */
class SechWaveFunction: public WaveFunctionBase {

    public:
        SechWaveFunction(const Path &,LookupTable &_lookup);
        ~SechWaveFunction();

        static const std::string name;
        std::string getName() const {return name;}

        double PsiTrial (const int);

    private:
        double a;           // The parameter of the wave function 

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
        JastrowWaveFunction(const Path &, LookupTable &_lookup);
        ~JastrowWaveFunction();

        double PsiTrial (const double);
        double PsiTrial (const int);
        double delPsiTrial(const double r);
        double delSqPsiTrial(const double r);
        double gradSqPsiTrial(const int);

        double twoBodyPsiTrial (const double);

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        double alpha;           // The parameter of the wave function
        double beta;            // The parameter of the wave function
    
};

// ========================================================================
// LiebLinigerWaveFunction Class
// ========================================================================
/**
 * Implementation of a Jastrow trial wave function suitable for He
 * @see Cuervo, Roy, & Boninsegni,, J. Chem. Phys., 122(11), 114504 (2005).
 */
class LiebLinigerWaveFunction: public WaveFunctionBase {
    
    public:
        LiebLinigerWaveFunction(const Path &, LookupTable &_lookup);
        ~LiebLinigerWaveFunction();

        double PsiTrial (const double);
        double PsiTrial (const beadLocator &bead1);
        double PsiTrial (const int);

        double delPsiTrial(const double r);
        double delSqPsiTrial(const double r);
        double gradSqPsiTrial(const int);

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        double R;           // The parameter length scale of the wave function
        double k;           // The wavevector of the wave function
    
};

// ========================================================================
// SutherlandWaveFunction Class
// ========================================================================
/**
 * Implementation of the Sutherland model exact wavefunction.  
 * @see Eq. (4) in Astrakharchik, G., Gangardt, D., Lozovik, Y. and Sorokin, I. 
 * Off-diagonal correlations of the Calogero-Sutherland model. Phys. Rev. E 74, 021105 (2006).
 */
class SutherlandWaveFunction: public WaveFunctionBase {
    
    public:
        SutherlandWaveFunction(const Path &, LookupTable &_lookup);
        ~SutherlandWaveFunction();

        /** The 2-body trial wavefunction */
        double PsiTrial(const double r) {return pow(2.0*sin(pioL*r),lambda);}
        double PsiTrial (const int);

        static const std::string name;
        std::string getName() const {return name;}
    
    private:
        double lambda;          // Sutherland model \lambda
        double pioL;            // pi / L
};

#endif
